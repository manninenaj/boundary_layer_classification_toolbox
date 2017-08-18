function [data,att,dim] = preprocessHalo(site,d_type,daten,pol_ch,save_bkg_corr,cut)
%preprocessHalo loada HALO lidar vertically pointing data and if needed,
%corrects it for the backrgound artifacts
%
%Inputs:
% site           'hyytiala' / 'juelich' / ...
% d_type         'ucorrected' / 'corrected' / 'calibrated'
% daten          date in matlab datenum format
% pol_ch         'co' / 'cr'
% save_bkg_corr  save or not the corrected data 
%
%Ouputs:
% data          data in struct format
% att           attributes of the date
% dim           dimensions of the data
%
%2017-07-03
%Antti Manninen
%antti.j.manninen@helsinki.fi
%University of Helsinki, Finland

%% Load HAlO data
[data,att,dim] = loadHaloVert(site,daten,pol_ch);
if isempty(data)
    return
else
    ind = find(data.range > cut,1,'first');
    data.cut = ind-1;
    att.cut.dimensions = {};
    att.cut.long_name = 'cut height';
    %% Correct bkg if needed
    if strcmp(d_type,'uncorrected')
        
        % Calculate background from the bkg files
        snr_corr_1 = correctHaloRipples(data.signal,data.time,...
            site,daten);
        
        % Correct shape
        [snr_corr_2, ~, ~, ~, ~] =...
            correctBackground(snr_corr_1, data.range, ...
            data.time, data.cut, 'correct_remnant','correct','ignore',60);
        %% Save old values for now
        switch site
            case 'hyytiala'
                fields_old = {'signal','beta','beta_raw','beta_error','v_error'};
            case 'juelich'
                fields_old = {'signal','beta','beta_raw','beta_error','v_error'};
                
            otherwise
                fields_old = {'signal','beta','beta_raw','beta_error','v_error'};
        end
        fields_new_co = {'signal0','beta0','beta_raw0','beta_error0','v_error0'};
        
        % Assign data and attributes to the new to be corrected fields
        for ifield = 1:length(fields_old)
            data.(fields_new_co{ifield}) = data.(fields_old{ifield});
            att.(fields_new_co{ifield})  = att.(fields_old{ifield});
        end
        
        % Reassign
        data.signal0 = data.signal; % new uncorrected
        data.signal  = snr_corr_2; % new corrected
        
        att.signal.comments = 'Background corrected signal in arbitrary units.';
        att.signal.long_name = 'Signal (bkg corrected)';
        
        %% Correct forcus and re-calculate uncertainties if system's specified
        if ismember(site,{'hyytiala','sodankyla','kuopio','kumpula',...
                'kenttarova','juelich','limassol','arm-oliktok'})
            
            % Re-apply focus correction to get new PR2 (beta)
            [data.beta_raw,data.beta_raw0] = correct_focus(data.focus, data);
            data.beta = data.beta_raw;
            
            % Calculate new uncertainties from new signal
            [data, att] = calculate_dl_SNR(data, att, site);
            data.beta_error = (1 / sqrt(data.num_pulses_m1)) .* ...
                (1 + ( 1 ./ abs(data.signal-1)));
           
        end
    end
    %% Use snr threshold
    switch site
        case 'hyytiala'
            snr_th = 1.01;
        case 'kenttarova'
            snr_th = 1.01;
        case 'limassol'
            snr_th = 1.01;
        case 'juelich'
            snr_th = 1.015;
        case 'arm-olitok'
            snr_th = 1.02;
        otherwise
            snr_th = 1.02;
    end
    data.beta(data.signal<snr_th | isnan(data.signal)) = nan;
    data.beta_error(data.signal<snr_th | isnan(data.signal)) = nan;
    data.v(data.signal<snr_th | isnan(data.signal)) = nan;
    data.v_error(data.signal<snr_th | isnan(data.signal)) = nan;
    
    %% Save corrected data
    if save_bkg_corr == 1
        
        data = rmfield(data,{'signal0','beta0','beta_raw0','beta_error0','v_error0'});
        att = rmfield(att,{'signal0','beta0','beta_raw0','beta_error0','v_error0'});
        
        f = fieldnames(att);
        for i=1:size(f,1)
            switch f{i}
             case 'global'
             case 'time'
              att.time = create_attributes({'time'},'Decimal hours UTC', 'hours since 00:00:00') ;
              att.time.standard_name = 'time';
              att.time.axis = 'T';
             case 'range'
              att.range.dimensions = {'range'};
             case {'azimuth'}
              if length(data.azimuth) == length(data.time)
                att.(f{i}).dimensions = {'time'};
                att.(f{i}).missing_value = -999;
                att.(f{i}).FillValue_ = -999;
              end
             case {'beta' 'beta_raw' 'v' 'v_raw' 'signal' 'v_error' 'beta_error'}
              att.(f{i}).dimensions = {'time' 'range'};
              att.(f{i}).missing_value = -999;
              att.(f{i}).FillValue_ = -999;
             otherwise
              att.(f{i}).dimensions = {};
            end
          if isfield(att.(f{i}), 'dimensions') && isempty(att.(f{i}).dimensions)
            if isfield(att.(f{i}),'missing_value')
              att.(f{i}) = rmfield(att.(f{i}),{'missing_value' 'FillValue_'}); 
            end
          end  
        end

        dimension.time = length(data.time);
        dimension.range = length(data.range);
        
        f_path = '/data/hatpro/jue/cloudnet';
        % Write corrected data into a netcdf file
        write_nc_struct([f_path '/' site '/calibrated/halo-doppler-lidar/' datestr(daten,'yyyy') '/' ...
            datestr(daten,'yyyymmdd') '_juelich_halo-doppler-lidar.nc'],dimension,data,att)
    end
end
end

