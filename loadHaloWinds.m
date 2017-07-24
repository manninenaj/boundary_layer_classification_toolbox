function [data,att,dim] = loadHaloWinds(site,daten,s_type)

%loadHaloWinds loads HALO Doppler lidar data from vertically pointing
%measurements
%
%Inputs:
% site           site name in lower case, e.g. 'hyytiala'
% daten          date in matlab datenum format
% s_type         scan type, e.g. 'VAD'
%
%Outputs:
% data          HALO data in struct format
% att
% dim

switch site
    case 'juelich'
        if strcmp(s_type,'dbs')
            scan = 'Wind_*';
        elseif strcmp(s_type,'vad')
            scan = 'wind_vad_*';
        end
        path_data = ['/data/TR32/D2/data/wind_lidar/data/nc/' ...
            datestr(daten,'yyyy') '/' datestr(daten,'mm'),'/',datestr(daten,'dd')];
        
    case 'arm-oliktok'
        if strcmp(s_type,'dbs')
            ifile_winds = []; path_data = []; scan = [];
        elseif strcmp(s_type,'vad')
            path_data = ['/data/hatpro/jue/cloudnet/juelich/calibrated/dopplerlidar/'...
            datestr(daten,'yyyy') '/ftp.cdc.noaa.gov/Public/mmaahn/olidlprofwind4newsM1'];
            scan = 'olidlprofwind4news*';
        end
end

% Find the date
files_winds = dir([path_data '/' scan]);
fnames_winds = {files_winds(~[files_winds.isdir]).name}'; % name list
if strcmp(s_type,'dbs') && strcmp(site,'juelich')
    ifile_profile = strfind(fnames_winds, 'Profile');
    ifile_vad = strfind(fnames_winds, 'vad');
    ifile_winds = find(cellfun('isempty', ifile_profile) & cellfun('isempty', ifile_vad));
else
    ifile_winds = strfind(fnames_winds, datestr(daten,'yyyymmdd')); % find dates
    ifile_winds = find(not(cellfun('isempty', ifile_winds)));
end
if isempty(ifile_winds)
    warning('%s %s data doesn''t exist --> skipping',...
        datestr(daten,'yyyymmdd'),s_type)
    data = []; att = []; dim = [];
    return
else
    % if more than one file per day
    if length(ifile_winds)>1
        % load first
        [data,att,dim] = load_nc_struct([path_data '/' ...
            fnames_winds{ifile_winds(1)}]);
        % get fieldnames
        fnamesdata = fieldnames(data);
        fnamesdims = fieldnames(dim);
        % load the rest
        for i1 = 2:length(ifile_winds)
            [tmpdata,~,tmpdim] = load_nc_struct([path_data '/' ...
                fnames_winds{ifile_winds(i1)}]);
            % find which field in dimensions is associated with time
            iftime = find(not(cellfun('isempty',...
                strfind(fnamesdims,'time'))),1);
            % add dimensions to 'time' field
            dim.(fnamesdims{iftime}) = dim.(fnamesdims{iftime}) + ...
                tmpdim.(fnamesdims{iftime});
            for i2 = 1:length(fnamesdata)
                % stack fields that are not scalars AND not 'range/height'
                if ~isscalar(tmpdata.(fnamesdata{i2})) ...
                        && (~strcmp(fnamesdata{i2},'range') ||...
                        ~strcmp(fnamesdata{i2},'height'))
                    % check which dim is associated with time
                    idim = ismember(size(tmpdata.(fnamesdata{i2})),...
                        length(tmpdata.time));
                    if idim(1)
                        % stack on the first file's fields
                        data.(fnamesdata{i2}) = [data.(fnamesdata{i2}); ...
                            tmpdata.(fnamesdata{i2})];
                    else
                        % stack on the first file's fields
                        data.(fnamesdata{i2}) = [data.(fnamesdata{i2}), ...
                            tmpdata.(fnamesdata{i2})];
                    end
                end
            end
        end
    % if only one file per day
    else
        [data,att,dim] = load_nc_struct([path_data '/' ...
            fnames_winds{ifile_winds(1)}]);
    end
    % modify fields when necessary
    switch site
        case 'hyytiala'
        case 'juelich'
            [~,~,~,Hrs,Mins,Secs] = jd2date(data.time);
            data.time = Hrs + Mins/60 + Secs/60/60;
            data.uwind = abs(data.speed) .* cos(degtorad(data.dir));
            data.vwind = abs(data.speed) .* sin(degtorad(data.dir));
            
            data.signal = data.intensity;
            data.wind_speed = data.speed;            
            data.wind_direction = data.dir;
%             switch s_type
%                 case 'dbs'
            data.range = data.height;
            data = rmfield(data,{'intensity','height','speed',...
                'dir','beta'});
%                 case 'vad'
%                     data = rmfield(data,{'intensity','speed',...
%                         'dir','w_speed','beta'});
%             end

            fn_j = fieldnames(data);
            for i_fn = 1:length(fn_j)
                if ~isvector(data.(fn_j{i_fn}))
                    data.(fn_j{i_fn}) = data.(fn_j{i_fn})';
                end
            end
        case 'limassol'
            data.range = data.height;
            data.uwind = data.u;
            data.uwind(data.uwind>14)=nan;
            data.vwind = data.v;
            data.vwind(data.vwind>7 | data.vwind<-16) = nan;
            data.wind_speed = sqrt(data.uwind.^2 + data.vwind.^2);
            data.wind_direction(data.wind_direction > 360) = nan;
            data.wind_speed(data.wind_direction == 141) = nan;
            data.wind_direction(data.wind_direction == 141) = nan;
            data = rmfield(data,{'u','height','v'});
            
        case 'arm-oliktok'
            data.time = data.time./3600;
            data.signal = data.mean_snr+1;
            data.uwind = data.u;
            data.vwind = data.v;
            data.range = data.height;
        otherwise
            error('Site %s not specified yet!',site)
    end
end
end
