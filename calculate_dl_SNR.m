function [data, attribute] = calculate_dl_SNR(data, attribute, site)

%  [data, attribute] = calculate_dl_SNR(data, attribute, site);
%
%  This function calculates the standard deviation of an individual
%  velocity estimate given the following inputs in the structure data.
%
%  INPUTS (defaults used in brackets if not supplied)
%  
%  Requires beta_raw, range and time   
%
%  bandwidth     (5e7)
%  elevation     (90)
%  lens_diameter (0.06)
%  wavelength    (1.5e-6)
%  focus         (-1)  infinity
%  energy        (1e-5)
%  prf           (15000)
%  num_pulses_m1 (50000)
%
%  OUTPUTS 
%  
%  The following outouts are added to the structure
%
%    v_error: Standard deviation of the velocity estimate, calculated
%             from the SNR variable following Rye and Hardesty (1997) 
%  
%    SNR:     Wideband SNR, calculated from beta_raw and the instrument
%             characteristics noted above
%
%    calculated_noise_floor:  The expected level of the noise at each height
%                             after performing range correction. For a
%                             focussed telescope, the range correction
%                             does not necessarily follow a range-squared
%                             law.  

if nargin < 1 
  help calculate_dl_SNR;
  return;
end

if ~isstruct(data)
  disp('Expecting a structure')
  help calculate_dl_SNR;
  return;
end

if ~isfield(data, 'range') 
  if isfield(data, 'height')
    range = data.height;
  else
    disp('Expecting a range variable')
    return;
  end
else
  range = data.range;
end

% defaults if not present
if ~isfield(data, 'bandwidth')
  data.bandwidth = 5e7;
  attribute.bandwidth.dimensions = {};
  attribute.bandwidth.long_name = 'bandwidth';
end
if ~isfield(data, 'elevation')
  data.elevation = 90;
  attribute.elevation.dimensions = {};
  attribute.elevation.long_name = 'elevation';
end
if ~isfield(data, 'lens_diameter')
  data.lens_diameter = 0.06;
% 0.075 *1/e2 - antenna not fully illuminated
% 0.04 (5cm antenna and 1/e2 = 0.8) for new version
  attribute.lens_diameter = create_attributes({},'diameter of lens','m');
end
if ~isfield(data, 'wavelength')
  data.wavelength = 1.5e-6;
  attribute.wavelength = create_attributes({},'laser wavelength','m');
end
if ~isfield(data, 'focus')
  data.focus = -1;
end
if ~isfield(data, 'energy')
  data.energy = 1e-5;
  attribute.energy = create_attributes({},'laser energy','J');
end
if ~isfield(data, 'prf')
  data.prf = 15000;
  attribute.prf = create_attributes({},'Pulse repetition frequency','Hz');
end
if ~isfield(data, 'num_pulses_m1')
  data.num_pulses_m1 = 50000;
end
if ~isfield(data, 'bandwidth')
    data.bandwidth = 5e+7;
    attribute.bandwidth = create_attributes({},'Bandwidth','Hz');
end

missing_value = -999;
if nargin < 3
  if findstr(lower(attribute.global.source),'chilbolton')
    site = 'chilbolton';
  else
    site = 'repartee';
  end  
end

switch lower(site)
 case 'hornisgrinde'
  % IMK Doppler lidar: WINDTRACER
  theory.PRF = 500; % PRF
  theory.Nyquist = 25;      % Nyquist
  theory.B = theory.Nyquist .* 2;
  theory.M = 48; % points per range gate
  theory.npulses = 100; % if dt = 0.2
  if isfield(data,'ratio_of_signal_to_noise')
    obs_signal = 10.^(data.ratio_of_signal_to_noise./10); 
  elseif isfield(data,'SNR')
    obs_signal = 10.^(data.SNR./10); 
  end
 case {'repartee','malmi'}
  % example for repartee data (and other salford campaigns)
  theory.PRF = 20000; % PRF
  theory.Nyquist = 14;      % Nyquist
  theory.B = theory.Nyquist .* 2;
  theory.M = 6; % points per range gate
  theory.npulses = 20000;
  obs_signal = abs(data.signal-1);
 
 case 'chilbolton'
  theory.PRF = 15000;
  theory.Nyquist = 12;      % Nyquist
  theory.B = theory.Nyquist .* 2;
  theory.M = 12;  
  theory.npulses = 20000 .* 5;
  obs_signal = abs(data.signal-1);
  
 case 'limassol'
  theory.PRF = 15000;
  theory.Nyquist = 19.456;      % Nyquist
  theory.B = theory.Nyquist .* 2;
  theory.M = 16;  
  theory.npulses = 30000;
  obs_signal = abs(data.signal-1);
 otherwise

%data.num_samples_gate = 10;
  switch lower(site)
   case 'juelich'
    corrfactor = 1.75; % system seems to under report number of pulses (says 1 s data, closer to 1.5 s data)
   otherwise
    corrfactor = 1;
  end

  %if these values not found, return
  if ~isfield(data,'num_samples_gate') | ~isfield(data,'prf') | ~isfield(data,'num_pulses_m1')
    disp(['Instrument for site: ' site ' not yet characterised.'])
    return;
  else
   theory.PRF = data.prf;
   theory.Nyquist = 19.456;      % Nyquist
   theory.B = theory.Nyquist .* 2;
   theory.M = data.num_samples_gate; % points per range gate
   theory.npulses = data.num_pulses_m1 .* corrfactor;
   obs_signal = abs(data.signal-1);
  end
end


theory.SNR = [10.^[-5:0.01:2]]';
theory.deltav = [1 1.5 2]; % typical signal spectral width - best value
                           % 1.5, or 2?

theory.ff = theory.deltav./theory.B;
theory.Np = theory.M .* theory.npulses .* theory.SNR;

for ii = 1:length(theory.deltav)
  % alpha = ratio of photon count to speckle count
  % alpha_speckle = mySNR.*6;
  theory.alpha_speckle(:,ii) = theory.SNR./(sqrt(2.*pi) .* theory.ff(ii));
  theory.v_err_approx(:,ii) = sqrt( sqrt(8) .* theory.ff(ii) .^2./(theory.Np .*theory.alpha_speckle(:,ii)) .* ((1 +  1./(sqrt(2.*pi)).*theory.alpha_speckle(:,ii)).^2)).*theory.B;
  theory.direct_detection(:,ii) = theory.deltav(ii)./sqrt(theory.Np);
  theory.v_err_pearson(:,ii) = 2.*sqrt(sqrt(pi)./theory.alpha_speckle(:,ii)).*(theory.deltav(ii)./sqrt(theory.Np)).*(1+(1./(2.*pi)).*theory.alpha_speckle(:,ii));
 
  
  % full treatment of gaussian pulse requires numerical integration
  % The limits should be -0.5 to 0.5. Integrate from 0 to 0.5 and
  % multiply by 2. This provides the Cramer-Rao lower bound (CRLB)
  x = 0:0.001:0.5;
  a = 10.^[-4:0.1:4];  
  g1 = zeros(size(a));
  f = zeros(size(a));
  for ialpha = 1:length(a)
    g1(ialpha) = a(ialpha)./sqrt(2.*pi).*2.*sum((x.^2.*exp(-(x.^2)))./(1 + a(ialpha).*(exp(-(x.^2)./2).^2))).*0.001;
    f(ialpha) = 2.*sum(((x./theory.ff(ii)).^2)./((1 + 1./(a(ialpha).*(exp(-(x.^2)./(2.*theory.ff(ii).^2))))).^2)).*0.001;
  end

  theory.f(:,ii) = interp1(a,f,theory.alpha_speckle(:,ii));
  theory.v_err_crlb(:,ii) = sqrt(theory.ff(ii).^2 ./ theory.npulses ./theory.M./ theory.f(:,ii)).*theory.B;
  
end

data.v_error = interp1(theory.SNR, theory.v_err_approx(:,3), obs_signal);
data.v_error(isnan(data.v_error)|~isfinite(data.v_error)|data.v_error > theory.Nyquist) = theory.Nyquist;
attribute.v_error = create_attributes({'time','range'},'Std of velocity', {'ms-1','m s<sup>-1</sup>'}, 'This variable is the standard deviation of the velocity estimate, calculated from the SNR variable following Rye and Hardesty (1997)');


