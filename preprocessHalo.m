function [data,att,dim] = preprocessHalo(site,d_type,daten,pol_ch,save_bkg_corr)
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
[data,att,dim] = loadHaloVert(site,daten,d_type,pol_ch);
if isempty(data)
    return
else
    %% Correct bkg if needed
    if strcmp(d_type,'uncorrected')
        
        % Calculate background from the bkg files
        snr_corr_1 = correctHaloRipples(data.signal,data.time,...
            site,daten);
        
        % Correct shape
        [snr_corr_2, ~, ~, ~, ~] =...
            correctBackground(snr_corr_1, data.range, ...
            data.time,'correct_remnant','correct','ignore',60);
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
                'kenttarova','juelich','limassol'})
            
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
        otherwise
            snr_th = 1.02;
    end
    data.beta(data.signal<snr_th | isnan(data.signal)) = nan;
    data.beta_error(data.signal<snr_th | isnan(data.signal)) = nan;
    data.v(data.signal<snr_th | isnan(data.signal)) = nan;
    data.v_error(data.signal<snr_th | isnan(data.signal)) = nan;
    
    %% Save corrected data
    if save_bkg_corr == 1
        
        % Pull the instrument id from the global attributes, if not found 'xx'
        fnames_att = fieldnames(att.global);
        im1 = not(cellfun('isempty',regexpi(fnames_att,'id')));
        im2 = not(cellfun('isempty',regexpi(fnames_att,'system')));
        
        % if field found with both 'id' & 'system' in name (case insensitive)
        if any(im1&im2)
            if ischar(att.global.(fnames_att{im1&im2}))
                % if systemID value is string and incase there are spaces
                halo_id = num2str(str2double(att.global.(fnames_att{im1&im2})));
            else
                % if systemID value is a number
                halo_id = num2str(att.global.(fnames_att{im1&im2}));
            end
        else
            % else assign to unknown
            halo_id = 'xx';
        end
        
        
        f_path = '/.../.../';
        % Write corrected data into a netcdf file
        write_nc_silent([f_path '/' site '/corrected/' datestr(daten,'yyyy') '/' ...
            pol_ch '/' datestr(daten,'yyyymmdd') '_halo-doppler-lidar-' halo_id ...
            '-' pol_ch '_corr.nc'],dim,data,att)
    end
end
end

