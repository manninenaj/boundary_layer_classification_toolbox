function [signal_corr, step_locations, cloud_mask, flags, background] =...
    correctBackground(signal, range, time, varargin)
%CORRECTBACKGROUND function corrects the background signal of the HALO
%Doppler lidar instrument. The background is corrected for step-changes and
%for the shape of the background within the step-changes respectively.
%
% Inputs:
%
% - signal              m-by-n matrix, the uncorrected SNR, 'm' equals the
%                       number of vertical profiles and 'n' equals the
%                       number of range bins (lengths of 'time' and 'range'
%                       inputs, respectively)
%
% - range               row or column vector, distances of each range bin
%                       from the instrument, length equals the number of
%                       range bins
%
% - time                row or column vector, time stamps of the vertical
%                       profiles, length equals the number of profiles
%
%
% Outputs:
%
% - signal_corrected    m-b-n matrix, the corrected SNR, dimensions equal
%                       the input 'signal'
%
% - step_locations      vector, indices of the detected step-changes in the
%                       bakcground signal
%
% - cloud_mask          m-by-n matrix, cloud and aerosol screening result,
%                       dimensions equal the input 'signal'
%
% - flags               vector, flags for the vertical profiles
%                       0 = constrained fit through 1 (SNR)
%                       1 = 1st deg poly fit
%                       2 = 2nd deg poly fit
%                       3 = correction failed, replace with nans
%                       4 = correction failed, replaced with original
%
% - background          n-by-m matrix, the extracted uncorrected background
%                       signal, dimension equal the input 'signal'
%
% Additional arguments follow in the form of property/value pairs.
% Valid properties are: 'win_size', 'n_sub_sect', 'wavelet_level',
% 'correct_remnant'
%
% - 'win_size'          size of the sliding window, which is used in
%                       calculating the 2D variance in variance based cloud
%                       screening. DEFAULT: ('win_size', [33 1])
%
% - 'wavelet_level'     number of iterations for the wavelet decomposition.
%                       DEFAULT: ('wavelet_level',5)
%
% - 'correct_remnant'   determines how to handle the remnant outlier
%                       profiles, which are not presented well by the
%                       averaged  approach in section 3. DEFAULT:
%                       ('correct_remnant', 'original')
%                       'correct'   correct all of the profiles using
%                                   robust linear regression including the
%                                   remnant profiles
%                       'remove'    remove only the remnant outlier
%                                   profiles, replace with nans
%                       'original'  replace the remnant outlier profiles
%                                   with original profiles
%                       'none'      no correction applied
%
% - 'ignore'            By default, ignore (3) of the lowest range gates
%                       due to incontamination by the emitted pulse.
%                       Additionally, with the 'ignore' parameter user can
%                       assign e.g. 30 of the lowest most range gates to
%                       nans when calculating the shape of the background.
%                       Despite the cloud screening, some remnant aerosol
%                       signal might cause errors in calculating the shape.
%
% Algorithm workflow (see more in detailed descriptions below)
% 0. Prepare data
% 1. Cloud screening
%      1.1 Crude cloud screening
%      1.2 Sub 2 km cloud screening
%      1.3 Final cloud screening
%      1.4 Calculation of the background shape
%      1.5 Filling the cloud-screened regions
% 2. Step detection
%      2.1 Multi-level 1D stationary wavelet decomposition
%      2.2 Peak detection i.e. locating step-changes
% 3. Correction for the step-changes and the shape of the background
% 4. OPTIONAL: Removal or correction of the remnant outlier profiles
%
% version 0.9.8
% 28 November 2016
% Antti Manninen
% antti.j.manninen@helsinki.fi

%% SET DEFAULTS
parameters.win_size        = [33 1];
parameters.wavelet_level   = 5;
parameters.correct_remnant = 'original';
parameters.ignore          = 3;
parameters.sizes           = [length(time) length(range)];
%% CHECK THE INPUTS

% Was the 'parameters' struct supplied?
if ~isempty(varargin)
    % Check for overrides of the defaults
    parameters = parsePropertyValuePairs(parameters, varargin);
end

% Check whether the 'parameters' values can be accepted
parameters = checkParameters(parameters);

% Check the dimensions of the first three inputs
if ~ismatrix(signal) && ~isvector(time) && ~isvector(range) &&...
        length(time)  ~= size(signal,1) &&...
        length(range) ~= size(signal,2)
    error 'Check the input dimensions!'
else
    
    %% PREPARE THE DATA
    signal_0 = signal;
    range = range(:);
    time = time(:);
    % By default, ignore three of the lowest range gates due to
    % incontamination by the emitted pulse.
    signal(:,1:4) = nan;
    signal(signal == 0) = nan;
    
    %% CLOUD SCREENING, cloud-aerosol masking and filling for wavelet
    [cloud_mask,signal_cld_scrd_outlr,signal_filled,flag_nofit] = ...
        cloudScreening(signal, parameters, range);
    
    %% MULTI-LEVEL 1D STATIONARY WAVELET DECOMPOSITION
    % The step-changes in the cloud-screened and filled signal are
    % detected by utilizing the stationary 1D wavelet decomposition
    % method. The output of the wavelet decomposition, detail
    % coefficients, describe the step-changes in the signal as peaks.
    % The detail coefficients are summed over each profile respectively
    % in order to make the peaks more pronounced and make the peak
    % detection more robust.
    
    % Use the default value 5 (recommended), OR use the value supplied
    % (suggest manual check)
    w_level = parameters.wavelet_level;
    
    % Initialize;
    detail_coeff_all  = [];
    i_switch_end   = 1; % for the zero padding
    
    % A range bin at a time
    for i_bin = 1:size(signal_filled, 2) * .75
        
        % Initialize
        zero_padded = signal_filled(:, end - i_bin);
        
        % The length of the zeropadded array ('zero_padded') has to be
        % divisible by (wavelet decomposition level)^2
        while round(size(zero_padded, 1) / 2^w_level) ~=...
                size(zero_padded,1 ) / 2^w_level
            
            % Pad signal with zeros to adjust its length
            zero_padded = [0; zero_padded; 0];
            
            % Remove one zero from the end and from the beginning
            % turn-by-turn
            if bitget(i_switch_end,1) % switch: odd or even
                
                % Remove zero from the end
                zero_padded(end) = [];
            else
                
                % Remove zero from the beginning
                zero_padded(1) = [];
            end
            
            % Add for iteration
            i_switch_end = i_switch_end + 1;
        end
        
        % Mark the locations of the added zeros
        cond_zeros = zero_padded == 0;
        
        % Convert zeros to nans for the wavelet function
        zero_padded(zero_padded == 0) = nan;
        
        % Call 'mySWT' function to calculate the stationary wavelet
        % decomposition up to the level indicated by 'w_level',
        % using the haar function for the convolution
        [~, detail_coeff] = mySWT(zero_padded, w_level, 'haar');
        
        % Leave only the highest level coefficients
        detail_coeff = detail_coeff(end,:);
        
        % Transpose to correct orientation
        detail_coeff = transpose(detail_coeff);
        
        % Get rid of nans which are padded zeros
        detail_coeff(cond_zeros) = [];
        
        % Concatenate arrays
        detail_coeff_all = horzcat(detail_coeff_all, detail_coeff);
    end
    
    % Sum up the coefficients along the bins
    detail_coeff_all_sum = nansum(abs(detail_coeff_all),2);
    
    % Convert zeros to nans
    detail_coeff_all_sum(detail_coeff_all_sum == 0) = nan;
    
    %% BACKGROUND STEP-CHANGE DETECTION
    % The step-changes in the signal are detected from inspecting the
    % detail coefficients from the output of the wavelet decomposition;
    % changes appear as peaks in the detail coefficients.
    % 'peakDetection' function detects peaks in signal based on given
    % delta value. The function searches for local minima and local
    % maxima from the signal (see more in detail description
    % 'peakDetection' function help).
    
    % The steps occur as peaks (or outliers) in the detail coeffs.,
    % so look for peaks in a subset of data. Define what is a "peak"
    % i.e. difference between a "valley" and a "peak".
    [max_peaks_final_level, ~] = peakDetection(detail_coeff_all_sum,...
        prctile(detail_coeff_all_sum, 75));
    if ~isempty(max_peaks_final_level)
        % Shift the locations based on the haar wavelet level impulse
        % distances
        step_locations = max_peaks_final_level(:, 1) + ...
            ((2^parameters.wavelet_level)/2) -1;
        step_locations(step_locations >= length(time) - 10) = [];
        
        %% CORRECT THE STEP-CHANGES AND THE SHAPE OF THE BACKGROUND
        % Removes the steps-changes and corrects the shape of the
        % background within the steps respectively. The correction is
        % carried out by fitting a surface to the cloud screened signal
        % within each step. The shape of the fitted surface is determined
        % by goodness-of-fit indicator root-mean-square error.
        
        % Initialize
        i_start_stp = 1;
        signal_shape_corrtd = nan(size(signal_0));
        background = nan(size(signal_0));
        flags = zeros(size(signal_0,1),1);
        
        for i_bkg_step = 1:length(step_locations) + 1
            
            % Determine background step end limit
            if i_bkg_step ~= length(step_locations) + 1
                i_end_stp = step_locations(i_bkg_step);
            else
                i_end_stp = size(signal, 1);
            end
            
            %% DRIFT CORRECTION
            % Median signal in range bins per each time stamp within a step
            med_signal_stp = nanmedian(signal_cld_scrd_outlr(...
                i_start_stp:i_end_stp,:),2);
            
            % Exclude nans
            x_bkg_stp = 1:length(time(i_start_stp:i_end_stp));
            x_bkg_stp(isnan(med_signal_stp)) = [];
            med_signal_stp(isnan(med_signal_stp)) = [];
            
            % Calculate linear fit and evaluate for the whole step
            linear_coeffs = my_robustfit(x_bkg_stp(:), med_signal_stp(:));
            p_coeff_drift = [linear_coeffs(2) linear_coeffs(1)];
            step_fitted = polyval(p_coeff_drift, ...
                1:length(time(i_start_stp:i_end_stp)));
            
            % Correct the drift with a step
            signal_drift_corrtd = signal_0(i_start_stp:i_end_stp,:) -...
                repmat(step_fitted(:),1,length(range)) +...
                p_coeff_drift(end);
            
            %% CORRECT BACKGROUND SHAPE AND STEP CHANGES
            % Determines the shape of the background between the
            % step-changes, and corrects for the shape and for the steps in
            % the background signal.
            
            % Exclude clouds
            signal_drift_corrtd(isnan(signal_cld_scrd_outlr(...
                i_start_stp:i_end_stp,:))) = nan;
            
            % Median signal per range bin within a step
            signal_med_bin = nanmedian(signal_drift_corrtd);
            
            % Determine is there enough data for higher order fit, depends
            % on the number and location of the remaining points. If there
            % are too few data points OR if the data points aren't
            % distributed sparsely enough, set flag --> no higher order
            % fit
            Signal_MedStep_ySel = signal_med_bin(1:length(range)*.5);
            range_cld_scrd_outlr_sel = range(1:length(range)*.5);
            if sum(~isnan(Signal_MedStep_ySel)) < 5 || ...
                    mean(range_cld_scrd_outlr_sel(...
                    ~isnan(Signal_MedStep_ySel))) > prctile(range,45)
                
                % Set flag
                flag_fit_step = 1;
            else
                flag_fit_step = 0;
            end
            
            % Select only non nans
            y_final_valid = signal_med_bin(~isnan(signal_med_bin));
            x_final_valid = range(~isnan(signal_med_bin));
            
            % Calculate 1st, 2nd deg, and constrained polynomial fits
            [B_prof_1deg,stats_1deg] = my_robustfit(x_final_valid(:), ...
                y_final_valid(:));
            [B_prof_2deg,stats_2deg] = my_robustfit([x_final_valid(:) ...
                x_final_valid(:).^2], y_final_valid(:));
            p_1deg_prof = [B_prof_1deg(2) B_prof_1deg(1)];
            p_2deg_prof = [B_prof_2deg(3) B_prof_2deg(2) B_prof_2deg(1)];
            % Evaluate coefficients
            y_fit_prof_1deg = polyval(p_1deg_prof, range);
            y_fit_prof_2deg = polyval(p_2deg_prof, range);

            % If 1st degree polynomial fit is better
            if stats_1deg.ols_s/stats_2deg.ols_s < 1 || ...
                    flag_fit_step == 1
                
                % Correct for step change and background shape
                signal_shape_corrtd(i_start_stp:i_end_stp,:) = ...
                    signal_drift_corrtd - repmat(y_fit_prof_1deg(:)',...
                    length(i_start_stp:i_end_stp), 1) + ...
                    p_1deg_prof(2);
                
                % Set flag to '1' if 1st deg fit was used
                flags(i_start_stp:i_end_stp,1) = ...
                    ones(length(i_start_stp:i_end_stp), 1);
                
                % Combine corrections to form corrected background
                background(i_start_stp:i_end_stp,:) = ...
                    (repmat(step_fitted(:),1,length(range)) - ...
                    p_coeff_drift(end))...
                    + repmat(y_fit_prof_1deg(:)', ...
                    length(i_start_stp:i_end_stp), 1);
                
                % If 2nd degree polynomial fit is better
            else
                % Correct for step change and background shape
                signal_shape_corrtd(i_start_stp:i_end_stp,:) = ...
                    signal_0(i_start_stp:i_end_stp,:) -...
                    repmat(y_fit_prof_2deg(:)', ...
                    length(i_start_stp:i_end_stp), 1) + ...
                    p_2deg_prof(3);
                
                % Set flag to '2' 2nd deg fit was used
                flags(i_start_stp:i_end_stp,1) = ...
                    repmat(2,length(i_start_stp:i_end_stp), 1);
                
                % Combine corrections to form corrected background
                background(i_start_stp:i_end_stp,:) = ...
                    (repmat(step_fitted(:),1,length(range)) - ...
                    p_coeff_drift(end))...
                    + repmat(y_fit_prof_2deg(:)', ...
                    length(i_start_stp:i_end_stp), 1);
            end
            
            % For iteration
            i_start_stp = i_end_stp + 1;
        end
        
        % Correct the step changes and the shape of the background
        signal_shape_corrtd = signal_0 - background + 1;
    else
        % If the step detection failed, skip it
        signal_shape_corrtd = signal_0;
        step_locations = [];
        flags = repmat(4,size(signal_0));
    end
    %% REMOVAL OR CORRECTION OF REMNANT OUTIER PROFILES
    % For instruments that are not operating optimally, the background
    % noise for some profiles may not be very well represented by the
    % the averaged approach. The outlier profiles can then be flagged
    % and rejected, or the user may choose to apply the background
    % noise profile shape detection and correction on a
    % profile-by-profile basis.
    
    % Initialise
    signal_remn = nan(size(signal_0));
    
    switch parameters.correct_remnant
        % 'correct': correct all of the profiles using robust linear
        %            regression including the remnant profiles
        % 'remove':  remove only the remnant outlier profiles
        % 'original':    first correct all profiles using robust linear
        %            regression, and then remove the outlier profiles
        % 'none':    no correction
        case 'correct'
            
            for i_remn = 1:size(signal_0,1)
                if sum(isnan(signal_shape_corrtd(i_remn,:))) ~= ...
                        length(signal_shape_corrtd(i_remn,:))
                    
                    % Select only background signal and apply cloudmask
                    % calculated from the shape corrected signal
                    y_remn = signal_shape_corrtd(i_remn,:);
                    y_remn(cloud_mask(i_remn,:)) = nan;
                    y_remn(1:parameters.ignore) = nan;
                    y_r_val = y_remn(~isnan(y_remn));
                    x_r_val = range(~isnan(y_remn));
                    
                    % Calculate robust bisquare linear fit
                    b_remn = my_robustfit(x_r_val(:),y_r_val(:));
                    p_c_remn = [b_remn(2) b_remn(1)];
                    y_f_r = polyval(p_c_remn,range);
                    signal_remn(i_remn,:) = ...
                        signal_shape_corrtd(i_remn,:) - y_f_r(:)' + 1;
                end
            end
            
            % Final corrected signal
            signal_corr = signal_remn;
            
        case 'remove'
            
            for i_remn = 1:size(signal_0,1)
                if sum(isnan(signal_shape_corrtd(i_remn,:))) ~= ...
                        length(signal_shape_corrtd(i_remn,:))
                    
                    % Select only background signal and apply cloudmask
                    % calculated from the shape corrected signal
                    y_remn = signal_shape_corrtd(i_remn,:);
                    y_remn(cloud_mask(i_remn,:)) = nan;
                    y_remn(1:parameters.ignore) = nan;
                    y_r_val = y_remn(~isnan(y_remn));
                    x_r_val = range(~isnan(y_remn));
                    
                    % Calculate robust bisquare linear fit
                    b_remn = my_robustfit(x_r_val(:),y_r_val(:));
                    p_c_remn = [b_remn(2) b_remn(1)];
                    y_f_r = polyval(p_c_remn,range);
                    signal_remn(i_remn,:) = ...
                        signal_shape_corrtd(i_remn,:) - y_f_r(:)' + 1;
                end
            end
            
            % Initialize with corrected signal
            signal_remn_cld_scrd = signal_remn;
            
            % Remove clouds, take median for each profile
            signal_remn_cld_scrd(isnan(signal_cld_scrd_outlr)) = nan;
            signal_remn_cld_scrd(:,1:parameters.ignore) = nan;
            signal_remn_peaks = nanmedian(signal_remn_cld_scrd,2);
            
            % Find the remant profiles with peak detection from the
            % median signal in profiles by using a 6 x std dev as peak
            % threshold
            [remn_peaks_max, remn_peaks_min] = ...
                peakDetection(signal_remn_peaks, ...
                6 * nanstd(signal_remn_peaks));
            
            % Collect min and max peak locations
            if ~isempty(remn_peaks_max) && ~isempty(remn_peaks_min)
                remn_peak_loc = vertcat(remn_peaks_max(:,1), ...
                    remn_peaks_min(:,1));
                
                % Set outlier profiles to nan
                signal_remn(remn_peak_loc,:) = nan;
                flags(remn_peak_loc) = 3;
            end
            
            % Final corrected signal
            signal_corr = signal_remn;
            
        case 'original'
            
            for i_remn = 1:size(signal_0,1)
                if sum(isnan(signal_shape_corrtd(i_remn,:))) ~= ...
                        length(signal_shape_corrtd(i_remn,:))
                    
                    % Select only background signal and apply cloudmask
                    % calculated from the shape corrected signal
                    y_remn = signal_shape_corrtd(i_remn,:);
                    y_remn(cloud_mask(i_remn,:)) = nan;
                    y_remn(1:parameters.ignore) = nan;
                    y_r_val = y_remn(~isnan(y_remn));
                    x_r_val = range(~isnan(y_remn));
                    
                    % Calculate robust bisquare linear fit
                    b_remn = my_robustfit(x_r_val(:),y_r_val(:));
                    [~, mID_r] = lastwarn; % iteration limit warning off
                    if ~isempty(mID_r), warning('off',mID_r), end
                    p_c_remn = [b_remn(2) b_remn(1)];
                    y_f_r = polyval(p_c_remn,range);
                    signal_remn(i_remn,:) = ...
                        signal_shape_corrtd(i_remn,:) - y_f_r(:)' + 1;
                end
            end
            
            % Initialize with corrected signal
            signal_remn_cld_scrd = signal_remn;
            
            % Remove clouds, take median for each profile
            signal_remn_cld_scrd(isnan(signal_cld_scrd_outlr)) = nan;
            signal_remn_cld_scrd(:,1:parameters.ignore) = nan;
            signal_remn_peaks = nanmedian(signal_remn_cld_scrd,2);
            
            % Find the remant profiles with peak detection from the
            % median signal in profiles by using a 6 x std dev as peak
            % threshold
            [remn_peaks_max, remn_peaks_min] = ...
                peakDetection(signal_remn_peaks, ...
                6 * nanstd(signal_remn_peaks));
            
            % Collect min and max peak locations
            if ~isempty(remn_peaks_max) && ~isempty(remn_peaks_min)
                remn_peak_loc = vertcat(remn_peaks_max(:,1), ...
                    remn_peaks_min(:,1));
                
                % Set outlier profiles to original
                signal_remn(remn_peak_loc,:) = signal(remn_peak_loc,:);
                flags(remn_peak_loc) = 4;
            end
            
            % Final corrected signal
            signal_corr = signal_remn;
            
        otherwise
            % includes option none %
            signal_corr = signal_shape_corrtd;
    end
    background = signal_0 - signal_corr;
end

% Subfunction - parsePropertyValuePairs
%D'Errico, John (2006). Parsing property/value pairs for function input
%(http://www.mathworks.com/matlabcentral/fileexchange/9082-parse-pv-pairs),
%MATLAB Central File Exchange. Retrieved Oct 7, 2015.

    function params = parsePropertyValuePairs(params, pv_pairs)
        %PARSEPROPERTYVALUEPAIRS: parses sets of property value pairs
        % usage: params = parse_pv_pairs(default_parameters, pv_pairs)
        %
        % arguments: (input)
        %  default_parameters - structure, with one field for every
        %   potential property/value pair. Each field will contain the
        %   default value for that property. If no default is supplied for
        %   a given property, then that field must be empty.
        %
        %  pv_array - cell array of property/value pairs.
        %   Case is ignored when comparing properties to the list of field
        %   names. Also, any unambiguous shortening of a field/property
        %   name is allowed.
        %
        % arguments: (output)
        %  params - parameter struct that reflects any updated
        %  property/value pairs in the pv_array.
        %
        % Example usage:
        % First, set default values for the parameters. Assume we
        % have four parameters that we wish to use optionally in
        % the function examplefun.
        %
        %  - 'viscosity', which will have a default value of 1
        %  - 'volume', which will default to 1
        %  - 'pie' - which will have default value 3.141592653589793
        %  - 'description' - a text field, left empty by default
        %
        % The first argument to examplefun is one which will always be
        % supplied.
        %
        %   function examplefun(dummyarg1,varargin)
        %   params.Viscosity = 1;
        %   params.Volume = 1;
        %   params.Pie = 3.141592653589793
        %
        %   params.Description = '';
        %   params=parse_pv_pairs(params,varargin);
        %   params
        %
        % Use examplefun, overriding the defaults for 'pie', 'viscosity'
        % and 'description'. The 'volume' parameter is left at its default.
        %
        %   examplefun(rand(10),'vis',10,'pie',3,'Description',
        %    'Hello world')
        %
        % params =
        %     Viscosity: 10
        %        Volume: 1
        %           Pie: 3
        %   Description: 'Hello world'
        %
        % Note that capitalization was ignored, the property 'viscosity'
        % was truncated as supplied. Also, note that the order the pairs
        % were supplied was arbitrary.
        
        npv = length(pv_pairs);
        n = npv/2;
        
        if n~=floor(n)
            error 'Property/value pairs must come in PAIRS.'
        end
        if n<=0
            % just return the defaults
            return
        end
        
        if ~isstruct(params)
            error 'No structure for defaults was supplied'
        end
        
        % there was at least one pv pair. process any supplied
        propnames = fieldnames(params);
        lpropnames = lower(propnames);
        for i_sf=1:n
            p_i = lower(pv_pairs{2*i_sf-1});
            v_i = pv_pairs{2*i_sf};
            
            ind = strmatch(p_i,lpropnames,'exact');
            if isempty(ind)
                ind = find(strncmp(p_i,lpropnames,length(p_i)));
                if isempty(ind)
                    error(['No matching property found for: ',...
                        pv_pairs{2*i_sf-1}])
                elseif length(ind)>1
                    error(['Ambiguous property name: ',pv_pairs{2*i_sf-1}])
                end
            end
            p_i = propnames{ind};
            
            % override the corresponding default in params
            params = setfield(params,p_i,v_i); %#ok
            
        end
    end

%% Subfunction - checkParameters
    function params = checkParameters(params)
        %CHECKPARAMETERS checks that the input parameters are correct
        
        % amount of range bins to be ignored at the bottom
        if isempty(params.ignore)
            params.ignore = 4;
        elseif not(isnumeric(params.ignore) && ...
                isscalar(params.ignore) && ...
                params.ignore <params.sizes(2) * .5)
            error(['''ignore'' parameter has to be numeric scalar and' ...
                ' cannot be larger than 1/2 times the number of' ...
                ' range gates. Default is %d'],3)
        end
        
        % win_size == 1 by default
        if isempty(params.win_size)
            params.win_size = [33 1];
        else
            if not(isnumeric(params.win_size) && ...
                    isvector(params.win_size) && ...
                    numel(params.win_size) == 2)
                error(['''win_size'' parameter must be a unmeric' ...
                    ' vector with length of two, default [33,1]'])
            end
        end
        
        % correct_remnant - must be one of 3 options
        valid = {'correct', 'remove', 'original', 'none'};
        if isempty(params.correct_remnant)
            % default == 'original'
            params.correct_remnant = 'original';
        end
        ind = find(strncmpi(params.correct_remnant,valid,...
            length(params.correct_remnant)));
        if (length(ind) == 1)
            params.correct_remnant = valid{ind};
        else
            error(['Invalid input for remnant correction parameter.' ...
                ' Valid options are: ''%s'', ''%s'', ''%s''' ...
                '(default), or ''%s'''],valid{1},valid{2},valid{3},...
                valid{4})
        end
        
        % wavelet_level == 5 by default
        if isempty(params.wavelet_level)
            params.wavelet_level = 5;
        elseif not(isnumeric(params.wavelet_level) && ...
                isscalar(params.wavelet_level)) && ...
                (length(params.wavelet_level)>1)
            error ''wavelet_level' must be a numeric scalar. Default is 5.'
        end
        
    end

%-- Subfunctions - CLOUD-AEROSOL MASKING
% Initial cloud-aerosol masking with 2d variance
% Divides the region of the calculated variance array contained in the
% furthest range bins (furthest 20%) into 64 subsections, which are used to
% find a dynamic threshold for the variance-based cloud-aerosol screening.
    function [cloud_mask,signal_cld_scrd_outlr,signal_fill,flag_nofit] =...
            cloudScreening(signal,params,range)
        % Pad with nans to avoid the border effect
        signal_pad_x = horzcat(signal, nan((params.win_size(1)-1)/2,...
            size(signal,1))');
        signal_pad   = horzcat(nan((params.win_size(1)-1)/2,...
            size(signal_pad_x,1))', signal_pad_x);
        
        % 2D running variance with a window having dimensions equal to
        % 'params.win_size'
        signal_var_pad = nan(size(signal_pad));
        for iC = 1+(params.win_size(1)-1)/2:...
                size(signal_pad,2)-(params.win_size(1)-1)/2
            temp_array = signal_pad(:,iC-(params.win_size(1)-1)/2:...
                iC+(params.win_size(1)-1)/2);
            signal_var_pad(:,iC) = nanvar(temp_array,[],2);
        end
        
        % Remove padded nans
        signal_var = signal_var_pad(:,(params.win_size(1)-1)/2+1:...
            end-(params.win_size(1)-1)/2);
        % Select variance in the furthest range bins (furthest 20%)
        signal_var_up20 = signal_var(:,end-round(size(signal,2)*.2)+1:end);
        
        % Find the nearest dimensions which are divisible by the number of
        % provided subsections (default == 8^2 = 64) AND are larger than
        % the original dimensions
        nearest_div_mm = size(signal_var_up20,1);
        while nearest_div_mm / 8 ~= ...
                fix(nearest_div_mm / 8)
            nearest_div_mm = nearest_div_mm + 1;
        end
        nearest_div_nn = size(signal_var_up20,2);
        while nearest_div_nn / 8 ~= ...
                fix(nearest_div_nn / 8)
            nearest_div_nn = nearest_div_nn + 1;
        end
        
        % If needed, pad with zeros so that the variance array can be
        % processed in equal-sized blocks
        if (size(signal_var_up20,2)-nearest_div_nn) ~= 0
            signal_var_up20_pad_range = horzcat(signal_var_up20,...
                zeros(size(signal_var_up20,1),...
                nearest_div_nn-size(signal_var_up20,2)));
        else
            signal_var_up20_pad_range = signal_var_up20;
        end
        if (size(signal_var_up20,1)-nearest_div_mm) ~= 0
            signal_var_up20_pad_all = vertcat(signal_var_up20_pad_range,...
                transpose(zeros(size(signal_var_up20_pad_range,2),...
                nearest_div_mm-size(signal_var_up20_pad_range,1))));
        else
            signal_var_up20_pad_all = signal_var_up20_pad_range;
        end
        
        % Construct indeces for processing the array in blocks, such as:
        %                 1 1 4 4 7 7
        % block_indeces = 2 2 5 5 8 8
        %                 3 3 6 6 9 9
        ind_1st_col = sort(repmat(1:8,...
            1, size(signal_var_up20_pad_all, 2)/8));
        ind_1st_block = repmat(ind_1st_col(:),...
            size(signal_var_up20_pad_all, 1)/8, 1);
        ind_vec_all = ind_1st_block;
        for i_sub_sect = 2:8
            ind_vec_all = vertcat(ind_vec_all, ind_1st_block + ...
                8 * (i_sub_sect - 1));
        end
        block_indeces = transpose(reshape(ind_vec_all,...
            [nearest_div_nn nearest_div_mm]));
        
        % Remove the regions for padded zeros
        if (abs(nearest_div_nn-size(signal_var_up20, 2))) ~= 0
            block_indeces(:, end - (nearest_div_nn -...
                size(signal_var_up20, 2)) + 1:end) = [];
        end
        if (abs(nearest_div_mm-size(signal_var_up20, 1))) ~= 0
            block_indeces(end - (nearest_div_mm -...
                size(signal_var_up20, 1)) + 1:end ,:) = [];
        end
        
        % Calculate median variance in each block
        block_median = nan(1,max(block_indeces(:)));
        for i_med = 1:max(block_indeces(:))
            block_median(i_med) = ...
                nanmedian(signal_var_up20_pad_all(block_indeces == i_med));
        end
        
        % Find the half of the blocks which have lowest median variance,
        % i.e. which are least influenced by clouds and aerosols, and
        % collect the calculated variance from those blocks into a
        % reference array for finding the dynamic threshold
        [~,sorted_block_ind] = sort(block_median);
        icond_refvar = ismember(block_indeces, sorted_block_ind(1:floor(...
            8 * 8/2)));
        signal_var_ref = signal_var_up20(icond_refvar);
        
        % Find dynamic threshold based on the number of screened pixels in
        % the lowest most median variance regions of the signal. The
        % allowed number of screened pixels is 0.5% of the total amount of
        % pixels in the reference region, variable 'ratio'
        th = 50; % initialise
        i = 1;
        signal_var_th = signal_var_ref;
        signal_var_th(signal_var_ref > prctile(signal_var(:),th)) = nan;
        ratio = sum(isnan(signal_var_th(:))) / numel(signal_var_ref);
        while ratio(i) > .005
            i = i + 1;
            % Increase threshold every iteration
            th(i) = th(i-1) + 1;
            signal_var_th = signal_var_ref;
            signal_var_th(signal_var_ref > ...
                prctile(signal_var(:),th(i))) = nan;
            ratio(i) = sum(isnan(signal_var_th(:))) / ...
                numel(signal_var_ref);
            if ratio(i) == ratio(i-1)
                break
            end
        end
        
        % Screen clouds based on dynamic variance threshold
        signal_cloud_scr_var = signal;
        signal_cloud_scr_var(signal_var > ...
            prctile(signal_var(:),th(i))) = nan;
        
        % Final cloud-aerosol screening and filling
        % Final cloud screening handles each vertical profile individually.
        % Each vertical profile is fitted with a 1st degree polynomial.
        % Note that the shape of the background can follow the shape of
        % either 1st or 2nd degree polynomial. However, the magnitude of
        % the change of the background as a function of range is
        % insignificant compared to the magnitude of the outliers caused
        % by the remnant cloud and aerosol signal. Thus, here the profiles
        % are fitted with 1st degree polynomials. The outliers are removed
        % with the 'findOutliers' function, which finds leverage points as 
        % well as points of high influence, and finally the indeces of the
        % outliers based on their Cook's distance. Then, the regions which
        % were cloud-screened have to be filled for the wavelet
        % decomposition. The profiles are fitted with either 1st or 2nd
        % degree polynomials based on the goodness-of-fit indicator,
        % root-mean-square error.
        
        % Initialize
        signal_cld_scrd_outlr = signal_cloud_scr_var;
        signal_fill = nan(size(signal));
        flag_nofit = nan(size(signal,1),1);
        for i_prof = 1:size(signal,1) % 24280 arm-graciosa 4 testing
            if sum(isnan(signal_cld_scrd_outlr(i_prof,:))) ~= ...
                    length(signal_cld_scrd_outlr(i_prof,:))
                % Find outlier indices
                [~, i_outlrs_hi, ~] = findOutliers(range, ...
                    signal_cld_scrd_outlr(i_prof,:), 1);
                
                % Remove outliers
                signal_cld_scrd_outlr(i_prof,i_outlrs_hi) = nan;
                y_outlr_rmvd = signal_cld_scrd_outlr(i_prof,:);
                x_outlr_rmvd = range;
                
                % Ignore additional number of range bins in case some
                % remnant aerosol signal did remain in the lower most range
                % bins. Note that the lower the fit reaches the more
                % accurate it is, tentatively.
                y_outlr_rmvd(1:params.ignore) = nan;
                
                % Select only non nans
                y_outlr_rmvd_valid = y_outlr_rmvd(~isnan(y_outlr_rmvd));
                x_outlr_rmvd_valid = x_outlr_rmvd(~isnan(y_outlr_rmvd));
                
                % Get 1st and 2nd deg polynomial fits
                [B_1deg_final_scrn,stats_final_1deg] = my_robustfit(...
                    x_outlr_rmvd(:), y_outlr_rmvd(:));
                p_1deg_final_scrn = [B_1deg_final_scrn(2) ...
                    B_1deg_final_scrn(1)];
                [B_2deg_final_scrn,stats_final_2deg] = my_robustfit(...
                    [x_outlr_rmvd_valid(:) x_outlr_rmvd_valid(:).^2],...
                    y_outlr_rmvd_valid(:));
                p_2deg_final_scrn = [B_2deg_final_scrn(3) ...
                    B_2deg_final_scrn(2) B_2deg_final_scrn(1)];
                
                % Evaluate for the whole range
                y_fit_final_scrn_1deg = polyval(p_1deg_final_scrn, range);
                y_fit_final_scrn_2deg = polyval(p_2deg_final_scrn, range);
                
                % Root Mean Squared Error (RMSE)
                RMSE_final_scrn_1deg = stats_final_1deg.ols_s;
                RMSE_final_scrn_2deg = stats_final_2deg.ols_s;
                
                % Initialize
                signal_fill(i_prof,:) = signal_cld_scrd_outlr(i_prof,:);
                
                % Ignore additional number of range bins in case some
                % remnant aerosol signal did remain in the lower most range
                % bins. Note that the lower the fit reaches the more
                % accurate it is.
                signal_fill(i_prof,1:params.ignore) = nan;
                
                % Determine if there's enough data for a fit, based on 
                % number of remaining points and location of the remaining 
                % points; if too few data points OR if the data points 
                % aren't distributed sparsely enough, set flag --> 'no fit'
                signal_cld_scrd_outlr_sel = ...
                    signal_cld_scrd_outlr(i_prof, 1:length(range) * .5);
                range_scrd_outlr_sel = range(1:length(range) * .5);
                if sum(~isnan(signal_cld_scrd_outlr_sel)) < 5 || ...
                        mean(range_scrd_outlr_sel(...
                        ~isnan(signal_cld_scrd_outlr_sel))) > ...
                        prctile(range,45)
                    
                    % Do not fit, filled signal remains with nans within
                    signal_fill(i_prof,...
                        isnan(signal_cld_scrd_outlr(i_prof,:))) = nan;
                    
                    % Set flag
                    flag_nofit(i_prof) = 1;
                else
                    
                    % If 1st degree polynomial fit is better
                    if RMSE_final_scrn_1deg/RMSE_final_scrn_2deg < 1.1
                        % Fill with 1st deg fits
                        signal_fill(i_prof,...
                            isnan(signal_fill(i_prof,:))) =...
                            y_fit_final_scrn_1deg(...
                            isnan(signal_fill(i_prof,:)));
                        
                        % If 2nd degree polynomial fit is better
                    else
                        % Fill with 2nd deg fits
                        signal_fill(i_prof,...
                            isnan(signal_fill(i_prof,:))) =...
                            y_fit_final_scrn_2deg(...
                            isnan(signal_fill(i_prof,:)));
                    end
                end
            end
        end
        % Cloud mask
        cloud_mask = isnan(signal_cld_scrd_outlr);
    end

%-- Subfunction - findOutliers
    function [i_outliers, i_outliers_high, i_outliers_low] =...
            findOutliers(x_outlr, y_outlr, adj)
        %OUTLIERS function first carries out a robust linear regression to 
        % the input data. Then it calculates the leverage points and the 
        % indeces of the outliers based on the Cook's distance with 
        % optional adjustment provided by the user.
        %
        % Inputs:
        % - x               x-values (vector), can contain nans
        % - y               y-values (vector), can contain nans
        % - adj             adjustment factor for Cook's distance threshold
        %                   (scalar)
        %
        % Outputs:
        % - iOutliers       indices of the outliers
        % - iOutliersHigh   indices of the outliers, > fitted y-values
        % - iOutliersLow    indices of the outliers, < fitted y-values
        
        % Check number of inputs
        if nargin ~= 3
            error 'Check number of inputs!'
        end
        
        % Check input validity
        if isscalar(adj) ~= 1
            error 'Cook''s distance input has to be scalar!'
        end
        
        % Convert to column vectors
        y_outlr = y_outlr(:); x_outlr = x_outlr(:);
        
        % Select non nan from the original data
        i_sel_org_outlr = find(~isnan(y_outlr)); % to be used later on
        y_val_org_outlr = y_outlr(i_sel_org_outlr);
        x_val_org_outlr = x_outlr(i_sel_org_outlr);
        
        % Calculate the coefficients
        b_outlr = my_robustfit(x_val_org_outlr(:),y_val_org_outlr(:));
        [~, mID] = lastwarn; % iteration limit warning turned off
        if ~isempty(mID)
            warning('off',mID)
        end
        p_coeff_outlr = [b_outlr(2) b_outlr(1)];
        
        % Evaluate polynomial along the whole range
        y_fit_outlr = polyval(p_coeff_outlr,x_outlr);
        
        % Model matrix
        X_outlr = [ones(length(x_val_org_outlr(:)),1), x_val_org_outlr(:)];
        
        % Hat matrix
        HatMatrix = X_outlr * ((X_outlr' * X_outlr) \ X_outlr');
        
        % Leverage points
        Leverage_points = diag(HatMatrix);
        
        % Residuals
        y_fit_valid_outlr = y_fit_outlr(~isnan(y_outlr)); % select non nan
        y_resid_outlr = y_val_org_outlr(:) - y_fit_valid_outlr(:);
        
        % Mean square error (MSE) of the predictor
        MSE_outlr = nan(size(y_val_org_outlr));
        for iMSE_outlr = 1:length(y_val_org_outlr)
            MSE_outlr(iMSE_outlr) = (y_fit_valid_outlr(iMSE_outlr) -...
                y_val_org_outlr(iMSE_outlr))^2;
        end
        MSE_outlr = sum(MSE_outlr) / numel(y_val_org_outlr);
        
        % Cook's distance
        D_Cook = nan(size(y_val_org_outlr));
        for iCook = 1:length(y_val_org_outlr)
            D_Cook(iCook) = y_resid_outlr(iCook)^2 / ...
                (numel(p_coeff_outlr) * MSE_outlr) * ...
                (Leverage_points(iCook)/(1 - Leverage_points(iCook))^2);
        end
        
        % Initialize
        i_outliers = i_sel_org_outlr;
        
        % Indices where Cook's distance is less than threshold,
        % here only the outliers remain while other are deleted
        i_outliers(D_Cook < (4/numel(y_val_org_outlr)) * adj) = [];
        
        % Indices where Cook's distance is less than threshold, and the
        % observed values are larger than the fitted values
        i_outliers_high = i_outliers(y_outlr(i_outliers) > ...
            y_fit_outlr(i_outliers));
        
        % Indices where Cook's distance is less than threshold, and the
        % observed values are lower than the fitted values
        i_outliers_low  = i_outliers(y_outlr(i_outliers) < ...
            y_fit_outlr(i_outliers));
        
    end

%-- Subfunction - mySWT
    function [approx_coeff, detail_coeff] = mySWT(x, w_level, wavelet_type)
        % Stationary wavelet decomposition for correcting background
        %
        % USAGE:
        %
        %   [approx_coeff, detail_coeff] = my_swt(x,w_level,wavelet_type);
        %
        % INPUTS
        %
        %    MY_SWT(X,LEVEL,'wavelet_type') computes the stationary wavelet
        %    decomposition of the signal X up to LEVEL, using wavelet of 
        %    type 'wavelet_type'. LEVEL must be a positive integer and X 
        %    should have a length divisible by 2^LEVEL.
        %
        %    Currently, 'Haar' is the only supported wavelet type.
        %
        %
        % OUTPUTS
        %
        %   Outputs are the approximation and detailed coefficients for 
        %   each level.
        %
        
        % Use row vector.
        x = x(:)';
        lx = length(x);
        % Check that the length of x is divisible by 2^(wavelet 
        % decomposition level)
        if rem(lx,2^w_level)>0
            disp(['SWT input should have length of divisible by ' ...
                '2^(wavelet decomposition level).' ...
                ' Something wrong when zero padding']);
            return
        end
        
        % Get decomposition filters.
        switch wavelet_type
            case 'haar'
                lopass = [ 1./sqrt(2) 1./sqrt(2)];
                hipass = [-1./sqrt(2) 1./sqrt(2)];
            otherwise
                disp('Ask for new version');
                return
        end
        
        % Compute stationary wavelet coefficients.
        approx_coeff = zeros(w_level,lx);
        detail_coeff = zeros(w_level,lx);
        
        for k = 1:w_level
            
            % Extension
            lf = length(lopass); % length of filter
            
            % Discrete Wavelet Transform mode is periodisation
            x = extendPeriodDWT(x,lf/2);
            
            % Decomposition
            detail_coeff(k,:) = extractVector(conv2(x(:)',hipass(:)',...
                'full'),lx,lf+1);
            approx_coeff(k,:) = extractVector(conv2(x(:)',lopass(:)',...
                'full'),lx,lf+1);
            
            % Dyadic upsampling of filters
            tmp = zeros(1,2.*lf);
            tmp(1:2:2 * lf) = lopass;
            lopass = tmp;
            tmp = zeros(1,2.*lf);
            tmp(1:2:2 * lf) = hipass;
            hipass = tmp;
            
            % Update x
            x = approx_coeff(k,:);
        end
    end

%-- Subfunction - extractVector
    function y = extractVector(x, len, start)
        %EXTRACTVECTOR extracts a vector from within a larger vector
        
        y = x;
        finish = start + len - 1;
        y = y(start:finish);
        
    end

%-- Subfunction - extendPeriodDWT
    function x = extendPeriodDWT(x,lf)
        %EXTENDPERIODDWT extends the DWT using periodisation
        
        length_of_x = length(x);
        if rem(length_of_x,2)
            x(length_of_x+1) = x(length_of_x);
            length_of_x = length_of_x+1;
        end
        
        I = [length_of_x-lf+1:length_of_x , 1:length_of_x , 1:lf];
        if length_of_x<lf
            I = mod(I,length_of_x);
            I(I==0) = length_of_x;
        end
        
        x = x(I);
        
    end

%-- Subfunction - peakDetection
% Modified after Billauer, Eli (2012). Peak detection using MATLAB
% (http://www.billauer.co.il/peakdet.html). Last accessed 7 Oct 2015

    function [max_tab_peakd, min_tab_peakd]= peakDetection(vect_peakd,...
            delta_peakd, ex_peakd)
        %PEAKDETECTION Detects peaks in a vector
        % [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local maxima and 
        % minima ("peaks") in the vector V. MAXTAB and MINTAB consists of 
        % two columns. Column 1 contains indices in V, and column 2 the 
        % found values.
        %
        % With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices in 
        % MAXTAB and MINTAB are replaced with the corresponding X-values.
        %
        % A point is considered a maximum peak if it has the maximal value,
        % and was preceded (to the left) by a value lower by DELTA.
        
        max_tab_peakd = [];
        min_tab_peakd = [];
        
        vect_peakd = vect_peakd(:); % in case this wasn't a proper vector
        
        % Check inputs
        if nargin < 3
            ex_peakd = (1:length(vect_peakd))';
        else
            ex_peakd = ex_peakd(:);
            if length(vect_peakd)~= length(ex_peakd)
                error(['Input vectors ''vect_peakd'' and ''ex_peakd'' ' ...
                    'must have the same length']);
            end
        end
        
        if (length(delta_peakd(:)))>1
            error('Input argument ''delta_peakd'' must be a scalar');
        end
        
        if delta_peakd <= 0
            error('Input argument ''delta_peakd'' must be positive');
        end
        
        % Prepare data
        mn_peakd = Inf;
        mx_peakd = -Inf;
        mnpos_peakd = nan;
        mxpos_peakd = nan;
        
        lookformax_peakd = 1;
        
        for i_peakd = 1:length(vect_peakd)
            this_peakd = vect_peakd(i_peakd);
            if this_peakd > mx_peakd
                mx_peakd = this_peakd;
                mxpos_peakd = ex_peakd(i_peakd);
            end
            if this_peakd < mn_peakd
                mn_peakd = this_peakd;
                mnpos_peakd = ex_peakd(i_peakd);
            end
            
            if lookformax_peakd
                if this_peakd < mx_peakd-delta_peakd
                    max_tab_peakd = [max_tab_peakd; mxpos_peakd mx_peakd];
                    mn_peakd = this_peakd; mnpos_peakd = ex_peakd(i_peakd);
                    lookformax_peakd = 0;
                end
            else
                if this_peakd > mn_peakd+delta_peakd
                    min_tab_peakd = [min_tab_peakd; mnpos_peakd mn_peakd];
                    mx_peakd = this_peakd; mxpos_peakd = ex_peakd(i_peakd);
                    lookformax_peakd = 1;
                end
            end
        end
    end
end
