function [data,att_dbs,att_vad] = combineHaloDBSnVAD(site,d_type,daten)
%combineHaloDBSnVAD combines dbs and vad winds. TBD: supplementary winds

% Load data
[dbs_tday,att_dbs] = loadHaloWinds(site,d_type,daten,'dbs');
[vad_tday,att_vad] = loadHaloWinds(site,d_type,daten,'vad');

% Check what data is there, if any
if isempty(dbs_tday) && ~isempty(vad_tday)      % dbs 0 AND vad 1
    status = 1;
elseif ~isempty(dbs_tday) && isempty(vad_tday)  % dbs 1 AND vad 0
    status = 2;
elseif isempty(dbs_tday) && isempty(vad_tday)   % dbs 0 AND vad 0
    status = 3;
elseif ~isempty(dbs_tday) && ~isempty(vad_tday) % dbs 1 AND vad 1
    status = 4;
end

% for testing
dbs = dbs_tday;
vad = vad_tday;

switch status
    case 1
        % Clean VAD-data
%         vad.wind_speed(:,1:3)     = nan;
%         vad.wind_direction(:,1:3) = nan;
%         vad.uwind(:,1:3)          = nan;
%         vad.vwind(:,1:3)          = nan;
%         vad.range(1:3)           = nan;
        data = vad;        
    case 2
        % Clean DBS data
%         dbs.wind_speed(:,1:4)     = nan;
%         dbs.wind_direction(:,1:4) = nan;
%         dbs.uwind(:,1:4)          = nan;
%         dbs.vwind(:,1:4)          = nan;
%         dbs.wwind(:,1:4)          = nan;
%         dbs.range(1:4)            = nan;
        data = dbs;
    case 3 % load model winds TBD
        data = [];
    case 4
        % Clean DBS data
        dbs.wind_speed(:,1:3)     = nan;
        dbs.wind_direction(:,1:3) = nan;
        dbs.uwind(:,1:3)          = nan;
        dbs.vwind(:,1:3)          = nan;
        dbs.wwind(:,1:3)          = nan;
        
        % Clean VAD-data
        vad.wind_speed(:,1:3)     = nan;
        vad.wind_direction(:,1:3) = nan;
        vad.uwind(:,1:3)          = nan;
        vad.vwind(:,1:3)          = nan;
                
        % Add new data from VAD
        % include small offset (fraction of dz) to avoid double-counting in range
        offset = median(diff(dbs.range)) * 0.05;
        ind = vad.range < (dbs.range(1) - offset);
        dbs.range = [vad.range(ind) ; dbs.range];
        
        % Interpolate VAD wind data to the common time vector.
        % Use nearest neighbour?
        % Using linear interpolation here:
        speed_i  = interp1(vad.time, vad.wind_speed,     dbs.time, 'linear', 'extrap');
        dir_i    = interp1(vad.time, vad.wind_direction, dbs.time, 'nearest');
        uspeed_i = interp1(vad.time, vad.uwind,          dbs.time, 'linear', 'extrap');
        vspeed_i = interp1(vad.time, vad.vwind,          dbs.time, 'linear', 'extrap');
        
        % find index where no direction data but speed data is found
        i_fd = find(~isnan(speed_i) & isnan(dir_i));
        
        % find the direction based on u ans v components
        % conversion factor
        r2d = 45.0/atan(1.0);
        % calculate dir
        dir_i(i_fd) = atan2(uspeed_i(i_fd), vspeed_i(i_fd)) * r2d + 180;
        
        % Combine wind speeds
        % (:) makes sure that wind/dir vector at surface is (117x1), not (1x117)
        windspd_tmp = [speed_i(:,ind)  dbs.wind_speed];
        winddir_tmp = [dir_i(:,ind)    dbs.wind_direction];
        uwind_tmp   = [uspeed_i(:,ind) dbs.uwind];
        vwind_tmp   = [vspeed_i(:,ind) dbs.vwind];
        data.range  = dbs.range;
        data.time   = dbs.time;
        
        % Try infilling some gaps
        windspd_smooth = medianfilter(windspd_tmp);
        winddir_smooth = medianfilter(winddir_tmp);
        uwind_smooth   = medianfilter(uwind_tmp);
        vwind_smooth   = medianfilter(vwind_tmp);
        [~, counts]    = medianfilter(isfinite(windspd_tmp));
        
        % kernel is [3 3], max count is 9
        windspd_smooth(counts < 3) = nan; windspd_smooth(:,1:3) = nan;
        winddir_smooth(counts < 3) = nan; winddir_smooth(:,1:3) = nan;
        uwind_smooth(counts < 3)   = nan; uwind_smooth(:,1:3)   = nan;
        vwind_smooth(counts < 3)   = nan; vwind_smooth(:,1:3)   = nan;
        ind_nm = isnan(windspd_tmp) & isfinite(uwind_smooth);
        data.wind_speed     = windspd_tmp;
        data.wind_direction = winddir_tmp;
        data.uwind          = uwind_tmp;
        data.vwind          = vwind_tmp;
        data.wind_speed(ind_nm(:,1:10))     = windspd_smooth(ind_nm(:,1:10));
        data.wind_direction(ind_nm(:,1:10)) = winddir_smooth(ind_nm(:,1:10));
        data.uwind(ind_nm(:,1:10))          = uwind_smooth(ind_nm(:,1:10));
        data.vwind(ind_nm(:,1:10))          = vwind_smooth(ind_nm(:,1:10));
        
        data.wind_speed(~isfinite(data.wind_speed)) = nan;
        data.wind_direction(~isfinite(data.wind_direction)) = nan;
        data.uwind(~isfinite(data.uwind)) = nan;
        data.vwind(~isfinite(data.vwind)) = nan;
        
end
end

