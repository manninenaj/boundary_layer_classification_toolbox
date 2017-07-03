function [data,att,dim] = loadHaloWinds(site,...
    d_type,daten,s_type)
%loadHaloWinds loads HALO Doppler lidar data from vertically pointing
%measurements
%
%Inputs:
% site           site name in lower case, e.g. 'hyytiala'
% d_type         data type, e.g. 'uncorrected'
% daten          date in matlab datenum format
% s_type         scan type, e.g. 'VAD'
%
%Outputs:
% data          HALO data in struct format
% att
% dim

path_data = ['/.../.../' site '/' d_type '/' datestr(daten,'yyyy') '/' s_type];

% Find the date
files_winds = dir([path_data '/' suffix]);
fnames_winds = {files_winds(~[files_winds.isdir]).name}'; % name list
ifile_winds = strfind(fnames_winds, datestr(daten,'yyyymmdd')); % find dates
ifile_winds = find(not(cellfun('isempty', ifile_winds)), 1);
if isempty(ifile_winds)
    warning('%s %s data doesn''t exist --> skipping',...
        datestr(daten,'yyyymmdd'),s_type)
    data = []; att = []; dim = [];
    return
else
    % if more than one file per day
    if length(ifile_winds)>1
        % load first
        [data,att,dim] = load_nc_struct_silent([path_data '/' ...
            fnames_winds{ifile_winds(1)}]);
        % get fieldnames
        fnamesdata = fieldnames(data);
        fnamesdims = fieldnames(dim);
        % load the rest
        for i1 = 2:length(ifile_winds)
            [tmpdata,~,tmpdim] = load_nc_struct_silent([path_data '/' ...
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
        [data,att,dim] = load_nc_struct_silent([path_data '/' ...
            fnames_winds{ifile_winds}]);
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
                'dir','w_speed','beta'});
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
        otherwise
            error('Site %s not specified yet!',site)
    end
end
end


