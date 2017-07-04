function [data,att,dim] = loadHaloVert(site,daten,d_type,pol_ch)
%loadHaloVert loads HALO Doppler lidar data from vertically pointing
%measurements
%
%Inputs:
% site           site name in lower case, e.g. site = 'hyytiala'
% daten          date in matlab format, e.g. daten = datenum(2017,06,01)
% d_type         data type, e.g. dtype = 'uncorrected'; dtype = 'calibrated'
% pol_ch         polarization channel mode, e.g. polch = 'co'
%
%Outputs:
% data          HALO data in struct format
%

% Different sites have different naming practises
switch site
    case {'hyytiala','kuopio','uto','kumpula','sodankyla','ny-alesund'}
        suffix = '*.nc';
    case 'juelich'
        suffix = 'Stare*.nc';
    case 'limassol'
        suffix = '*.nc';
    case 'arm-oliktok'
        suffix = '*.cdf';
    case 'arm-graciosa'
        suffix = '*.cdf';
    otherwise
        suffix = '*.nc';
end

% Construct path to file
path_data = ['/.../.../' site '/' d_type ...
    '/' datestr(daten,'yyyy') '/'  pol_ch];

% Find the date
files_snr = dir([path_data '/' suffix]); % vertical
fnames_snr = {files_snr(~[files_snr.isdir]).name}'; % name list
ifile_snr = strfind(fnames_snr, datestr(daten,'yyyymmdd')); % find dates
ifile_snr = find(not(cellfun('isempty', ifile_snr)));
if isempty(ifile_snr)
    warning('%s %s-pol data doesn''t exist --> skipping',...
        datestr(daten,'yyyymmdd'),pol_ch)
    data = []; att = []; dim = [];
    return
else
    % if more than one file per day
    if length(ifile_snr)>1
        % load first
        [data,att,dim] = load_nc_struct_silent([path_data '/' ...
            fnames_snr{ifile_snr(1)}]);
        % get fieldnames
        fnamesdata = fieldnames(data);
        fnamesdims = fieldnames(dim);
        % transpose so that time is associated with rows
        for i0 = 1:length(fnamesdata)
            % stack fields that are not scalars AND not 'range'
            if ~isscalar(data.(fnamesdata{i0})) ...
                    && ~strcmp(fnamesdata{i0},'range')
                % check which dim is associated with time
                idim = ismember(size(data.(fnamesdata{i0})),...
                    length(data.time));
                if idim(2)
                    % transpose and stack on the first file's fields
                    data.(fnamesdata{i0}) = data.(fnamesdata{i0})';
                end
            end
        end
        % load the rest
        for i1 = 2:length(ifile_snr)
            [tmpdata,~,tmpdim] = load_nc_struct_silent([path_data '/' ...
                fnames_snr{ifile_snr(i1)}]);
            % find which field in dimensions is associated with time
            iftime = find(not(cellfun('isempty',...
                strfind(fnamesdims,'time'))),1);
            % add dimensions to 'time' field
            dim.(fnamesdims{iftime}) = dim.(fnamesdims{iftime}) + ...
                tmpdim.(fnamesdims{iftime});
            for i2 = 1:length(fnamesdata)
                % stack fields that are not scalars AND not 'range'
                if ~isscalar(tmpdata.(fnamesdata{i2})) ...
                        && ~strcmp(fnamesdata{i2},'range')
                    % check which dim is associated with time
                    idim = ismember(size(tmpdata.(fnamesdata{i2})),...
                        length(tmpdata.time));
                    if idim(1)
                        % stack on the first file's fields
                        data.(fnamesdata{i2}) = [data.(fnamesdata{i2}); ...
                            tmpdata.(fnamesdata{i2})];
                    else
                        % transpose and stack on the first file's fields
                        data.(fnamesdata{i2}) = [data.(fnamesdata{i2}); ...
                            tmpdata.(fnamesdata{i2})'];
                    end
                end
            end
        end
        % if only one file per day
    else
        [data,att,dim] = load_nc_struct_silent([path_data '/' ...
            fnames_snr{ifile_snr}]);
    end
    
    switch site
        case 'juelich'
            
            % Assign the missing names of which were different before
            data.azimuth = data.azi;
            data.beta_raw = data.beta;
            data.signal = data.intensity;
            data.v = data.doppler;
            data.v_raw = data.doppler;
            data.num_pulses_m1 = data.pulses_per_ray;
            data.focus = data.focus_range;
            data.num_samples_gate = data.gate_points;
            att.azimuth = att.azi;
            att.beta_raw = att.beta;
            att.signal = att.intensity;
            att.v = att.doppler;
            att.v_raw = att.doppler;
            att.num_pulses_m1 = att.pulses_per_ray;
            att.focus = att.focus_range;
            att.num_samples_gate = att.gate_points;
            data = rmfield(data,{'intensity','doppler','pulses_per_ray',...
                'focus_range','azi','ele','gate_points'});
            att = rmfield(att,{'intensity','doppler','pulses_per_ray',...
                'focus_range','azi','ele','gate_points'});
            
            % Juelich time stamps are in Julian day numbers
            [~, ~, ~, hour, minute, second] = jd2date(data.time);
            data.time = hour + minute/60 + second/60/60;
            data.time(isnan(nansum(data.signal,1))) = [];
            
            % Sometimes last point is from next day
%             data.time(data.time < 0 | data.time > 24) = nan;
            if data.time(end) < data.time(end - 1)
                data.time(end) = data.time(end) + 24;
            end
            
            % Calculate uncertanties for the uncorrected data
            [data, att] = calculate_dl_SNR(data, att, site);
            data.beta_error = (1 / sqrt(data.num_pulses_m1))...
                .* (1 + ( 1 ./ abs(data.signal-1)));
            att.beta_error = att.beta;
            att.beta_error.long_name = [att.beta_error.long_name ...
                ' error estimate'];
            
        case {'arm-oliktok','arm-graciosa'}
            % Assign the missing names of which were different before
            data.beta_raw = data.attenuated_backscatter;
            data.beta = data.beta_raw;
            data.signal = data.intensity;
            data.v = data.radial_velocity;
            data.v_raw = data.radial_velocity;
            data.num_pulses_m1 = str2double(att.global.shots_per_profile);
            data.focus = str2double(att.global.focus_range);
            data.num_samples_gate = str2double(att.global.samples_per_gate);
            % No worries about attributes Ewan'll give better load function
            att.beta_raw = att.attenuated_backscatter;
            att.beta = att.beta_raw;
            att.signal = att.intensity;
            att.v = att.radial_velocity;
            att.v_raw = att.radial_velocity;
            
            data = rmfield(data,{'intensity','radial_velocity',...
                'attenuated_backscatter'});
            att = rmfield(att,{'intensity','radial_velocity',...
                'attenuated_backscatter'});
            
            % arm-oliktok time stamps are in Julian day numbers
            data.time = data.time/3600;
            
            % Sometimes last point is from next day
%             data.time(data.time < 0 | data.time > 24) = nan;
            if data.time(end) < data.time(end - 1)
                data.time(end) = data.time(end) + 24;
            end
            
            % Calculate uncertanties for the uncorrected data
            [data, att] = calculate_dl_SNR(data, att, site);
            data.beta_error = (1 / sqrt(data.num_pulses_m1))...
                .* (1 + ( 1 ./ abs(data.signal-1)));
            att.beta_error = att.beta;
            att.beta_error.long_name = [att.beta_error.long_name ...
                ' error estimate'];
            
        case 'hyytiala'
            data.v = data.v_raw; % bkg noise th. is applied -> undo
        case 'kenttarova'
            data.v = data.v_raw; % bkg noise th. is applied -> undo
        case 'limassol'
            % Calculate uncertainties for the uncorrected signal
            [data, att] = calculate_dl_SNR(data, att, site);
            data.beta_error = (1 / sqrt(data.num_pulses_m1)) .* ...
                (1 + ( 1 ./ abs(data.signal-1)));
    end
    % add write nc struct
end
end


