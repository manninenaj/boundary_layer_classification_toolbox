function [epsilon,epsilon_error,epsilon_time,mean_v,mean_v_error,...
    sigma_v,sigma_v_error,mean_b,mean_b_error,L1,L2,nsamples,signal] = ...
    calcTKE(data_vert,data_wind_tday,data_wind_yday,...
    data_wind_tmrw,dt)
%CALCTKE

% Initialize
epsilon = struct;
epsilon_error = struct;
epsilon_time = struct;
signal = struct;

% Defaults
beamwidth = 0.00015; % approx two-way half power full width in degrees
wind_error = 1.5; % ms-1

% Delta time
dt_raw = median(diff(data_vert.time));

% Different lidars have different resolutions
if isempty(dt)
    if dt_raw > 0.003
        dt = [2 3 5 15 30 60]./60; % selection in mins, output in decimal hrs
    else
        dt = [1 2 3 5 15 30 60]./60; % selection in mins, output in decimal hrs
    end
else
    dt = dt./60;
end

% Calculate sigma_v and mean_v at different temporal resolutions
for ii = 1:length(dt)
    atime = dt(ii)/2:dt(ii):24 - dt(ii)/2;
    mean_v_tmp        = zeros(length(atime), length(data_vert.range));
    mean_v_error_tmp  = zeros(length(atime), length(data_vert.range));
    sigma_v_tmp       = zeros(length(atime), length(data_vert.range));
    sigma_v_error_tmp = zeros(length(atime), length(data_vert.range));
    mean_b_tmp        = zeros(length(atime), length(data_vert.range));
    mean_b_error_tmp  = zeros(length(atime), length(data_vert.range));
    nsamples_tmp      = zeros(length(atime), length(data_vert.range));
    mean_signal_tmp   = zeros(length(atime), length(data_vert.range));
    for kk = 1:length(atime)
        time_box = find(data_vert.time > atime(kk)-dt(ii)./2 ...
            & data_vert.time <= atime(kk)+dt(ii)./2);
        if length(time_box) > 1
            tmp_v = reshape(data_vert.v(time_box,1:length(data_vert.range)), length(time_box), length(data_vert.range));
            tmp_v_error = reshape(data_vert.v_error(time_box,1:length(data_vert.range)), length(time_box), length(data_vert.range));
            tmp_b = reshape(data_vert.beta(time_box,1:length(data_vert.range)), length(time_box), length(data_vert.range));
            tmp_b_error = reshape(data_vert.beta_error(time_box,1:length(data_vert.range)), length(time_box), length(data_vert.range));
            tmp_snr = reshape(data_vert.signal(time_box,1:length(data_vert.range)), length(time_box), length(data_vert.range));
            
            nsamples_tmp(kk,:) = sum(~isnan(tmp_v));
            index = find(nsamples_tmp(kk,:) > min(10,size(tmp_v,1) .* .75));
            if length(index) >= 1
                mean_v_tmp(kk,index)        = nanmean(tmp_v(:,index));
                mean_v_error_tmp(kk,index)  = nanmean(tmp_v_error(:,index));
                sigma_v_tmp(kk,index)       = nanstd(tmp_v(:,index));
                sigma_v_error_tmp(kk,index) = nanstd(tmp_v_error(:,index));
                mean_b_tmp(kk,index)        = nanmean(tmp_b(:,index));
                mean_b_error_tmp(kk,index)  = nanmean(tmp_b_error(:,index));
                mean_signal_tmp(kk,index)   = nanmean(tmp_snr(:,index));
            end
        end
    end
%     atime(sum(mean_v_tmp==0,2) == size(mean_v_tmp,2)) = [];
%     mean_v_tmp(sum(mean_v_tmp==0,2) == size(mean_v_tmp,2),:) = [];
%     mean_v_error_tmp(sum(mean_v_error_tmp==0,2) == size(mean_v_error_tmp,2),:) = [];
%     sigma_v_tmp(sum(sigma_v_tmp==0,2) == size(sigma_v_tmp,2),:) = [];
%     sigma_v_error_tmp(sum(sigma_v_error_tmp==0,2) == size(sigma_v_error_tmp,2),:) = [];
    atime(atime==0) = nan;
    mean_v_tmp(mean_v_tmp==0) = nan;
    mean_v_error_tmp(mean_v_error_tmp==0) = nan;
    sigma_v_tmp(sigma_v_tmp==0) = nan;
    sigma_v_error_tmp(sigma_v_error_tmp==0) = nan;
    mean_b_tmp(mean_b_tmp==0) = nan;
    mean_b_error_tmp(mean_b_error_tmp==0) = nan;
    mean_signal_tmp(mean_signal_tmp==0) = nan;
    
    if dt(ii)==0.5/60
        mean_v.t_0_5min = mean_v_tmp;
        mean_v_error.t_0_5min = mean_v_error_tmp;
        sigma_v.t_0_5min = sigma_v_tmp;
        sigma_v_error.t_0_5min = sigma_v_error_tmp;
        mean_b.t_0_5min = mean_b_tmp;
        mean_b_error.t_0_5min = mean_b_error_tmp;
        nsamples.t_0_5min = nsamples_tmp;
        time_ref.t_0_5min = atime;
        signal.t_0_5min = mean_signal_tmp;
    else
        eval(['mean_v.t_' sprintf('%s',num2str(dt(ii)*60)) 'min = mean_v_tmp;'])
        eval(['mean_v_error.t_' sprintf('%s',num2str(dt(ii)*60)) 'min = mean_v_error_tmp;'])
        eval(['sigma_v.t_' sprintf('%s',num2str(dt(ii)*60)) 'min = sigma_v_tmp;'])
        eval(['sigma_v_error.t_' sprintf('%s',num2str(dt(ii)*60)) 'min = sigma_v_error_tmp;'])
        eval(['mean_b.t_' sprintf('%s',num2str(dt(ii)*60)) 'min = mean_b_tmp;'])
        eval(['mean_b_error.t_' sprintf('%s',num2str(dt(ii)*60)) 'min = mean_b_error_tmp;'])
        eval(['nsamples.t_' sprintf('%s',num2str(dt(ii)*60)) 'min = nsamples_tmp;'])
        eval(['time_ref.t_' sprintf('%s',num2str(dt(ii)*60)) 'min = atime;'])
        eval(['signal.t_' sprintf('%s',num2str(dt(ii)*60)) 'min = mean_signal_tmp;'])
    end
end

% for later use
clearvars sigma_v_tmp sigma_v_error_tmp mean_v_tmp mean_v_error_tmp mean_b_tmp mean_b_error_tmp mean_signal_tmp

% Calculate epsilon for first value of dt only unless high res data available
% if dt_raw > 0.003
%     timeperiod = 1:3;
% else
%     timeperiod = 1:4;
% end
for ii = 1:length(dt)%timeperiod
    
    % require sigma_v, v_error and n_samples
    if ~isfield(data_vert,'v_error') && ~isfield(data_vert,'n_samples')
        disp('The parameters "v_error" and "n_samples" are also require to estimate dissipation rate');
        return;
    end
    
    if dt(ii) == 0.5/60
        sigma_v_ii = sigma_v.t_0_5min;
        mean_v_error_ii = mean_v_error.t_0_5min;
        time_ii = time_ref.t_0_5min;
    else
        eval(['sigma_v_ii = sigma_v.t_' sprintf('%s',num2str(dt(ii)*60)) 'min;'])
        eval(['mean_v_error_ii = mean_v_error.t_' sprintf('%s',num2str(dt(ii)*60)) 'min;'])
        eval(['time_ii = time_ref.t_' sprintf('%s',num2str(dt(ii)*60)) 'min;'])
    end
    n_samples = data_vert.num_samples_gate;
        
    obs_variance = sigma_v_ii.^2;
    noise_variance = mean_v_error_ii.^2;
    true_variance = obs_variance - noise_variance;
    
    [Xr,Yr] = meshgrid(data_vert.range, time_ii);
    if ~isempty(data_wind_tmrw) && ~isempty(data_wind_yday)
        [Xo,Yo] = meshgrid(data_wind_tday.range,...
            [data_wind_yday.time(:)-24;...
            data_wind_tday.time(:);...
            data_wind_tmrw.time(:)+24]);
        uwind_tmp = interp2(Xo,Yo,...
            [data_wind_yday.uwind;...
            data_wind_tday.uwind;...
            data_wind_tmrw.uwind],Xr,Yr);
        vwind_tmp = interp2(Xo,Yo,[data_wind_yday.vwind;...
            data_wind_tday.vwind;...
            data_wind_tmrw.vwind],Xr,Yr);
    else
        [Xo,Yo] = meshgrid(data_wind_tday.range, data_wind_tday.time);
        uwind_tmp = interp2(Xo,Yo,data_wind_tday.uwind,Xr,Yr);
        vwind_tmp = interp2(Xo,Yo,data_wind_tday.vwind,Xr,Yr);
    end
    
    % Try infilling some gaps
    usmooth = medianfilter(uwind_tmp);
    vsmooth = medianfilter(vwind_tmp);
%     vsmooth = reshape(vsmooth,size(vwind_tmp,1),size(vwind_tmp,2));
    [~, count] = medianfilter(isfinite(uwind_tmp));
    
    % kernel is [3 3], max count is 9
    usmooth(count < 3) = NaN;
    ind_nm = isnan(uwind_tmp) & isfinite(usmooth);
    uwind = uwind_tmp; vwind = vwind_tmp; 
    uwind(ind_nm) = usmooth(ind_nm);
    vwind(ind_nm) = vsmooth(ind_nm);
    
    % nans
    uwind(isnan(sigma_v_ii))=nan;
    vwind(isnan(sigma_v_ii))=nan;
    
    wind = sqrt(uwind .* uwind + vwind .* vwind);
    
    beamdiameter = ones(length(time_ii),1) * data_vert.range' .* (beamwidth*pi/180);
    
    L1_ii = ((3600.*dt(ii).*wind) + 2 .* (ones(length(time_ii),1) * data_vert.range') .* sin(beamdiameter./2));
    L2_ii = ((3600.*dt_raw.*wind) +  2 .* (ones(length(time_ii),1) * data_vert.range') .* sin(beamdiameter./2));
    
%     L1(isnan(L1)) = 0; L2(isnan(L2)) = 0;
    true_variance = abs(true_variance);
    true_variance(isnan(true_variance) | true_variance == 0) = nan;
    
    epsilon_tmp = (2/(3.*.55)).^1.5 .* true_variance.^(3/2) .* 2.*pi .* ((L1_ii.^(2/3) - L2_ii.^(2/3)).^(-3/2));
    epsilon_error_tmp = 3.*sqrt(sqrt(4./n_samples.*(mean_v_error_ii.^2)./(sigma_v_ii.^2))) + wind_error ./ 10;
    
    % final clean
    index = find(epsilon_tmp <= 0 | epsilon_error_tmp > 10);
    epsilon_tmp(index) = NaN;
    epsilon_error_tmp(index) = NaN;
    
    if dt(ii)==.5/60
        epsilon.t_0_5min = epsilon_tmp;
        epsilon_error.t_0_5min = epsilon_error_tmp;
        epsilon_time.t_0_5min = time_ii;
        L1.t_0_5min = L1_ii;
        L2.t_0_5min = L2_ii;
    else
        eval(['epsilon.t_' sprintf('%s',num2str(dt(ii)*60)) 'min = epsilon_tmp;'])
        eval(['epsilon_error.t_' sprintf('%s',num2str(dt(ii)*60)) 'min = epsilon_error_tmp;'])
        eval(['epsilon_time.t_' sprintf('%s',num2str(dt(ii)*60)) 'min = time_ii;'])
        eval(['L1.t_' sprintf('%s',num2str(dt(ii)*60)) 'min = L1_ii;'])
        eval(['L2.t_' sprintf('%s',num2str(dt(ii)*60)) 'min = L2_ii;'])
    end
    
end
end
