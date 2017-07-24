function [time_o,beta,skewn,epsilon,epsilon_error,shear_vec,shear_dir,...
    aero_top,velo,velo_error,sigma_w,sigma_w_error,aero_layer_mask,...
    nsamples,speed,direc,signal] = calcWindQuantities(data_vert,data_wind_tday,...
    data_wind_yday,data_wind_tmrw,model_wind_tday,model_wind_yday,...
    model_wind_tmrw,dt,site,cut_h)
%CALCWINDQUANTITIES calculates different quantities from Halo Doppler wind
%lidar data
%
%Antti Manninen
%University of Helsinki
%Department of Physics
%antti.j.manninen@helsinki.fi

%% Calculate TKE, velo, beta
% Load vertical data and wind data from yesterday, today and tomorrow
[epsilon,epsilon_error,t_epsilon,velo,velo_error,sigma_w,sigma_w_error,...
    beta,~,~,nsamples,signal] = calcTKE(data_vert,data_wind_tday,...
    data_wind_yday,data_wind_tmrw,dt);
% Fill missing lidar TKE with model TKE
fnames = fieldnames(beta);
if ~isempty(model_wind_tday)
    eps_model = calcTKE(data_vert,model_wind_tday,...
        model_wind_yday,model_wind_tmrw,dt);
    for ifn = 1:length(fnames)
        cond = isnan(epsilon.(fnames{ifn})) & ~isnan(eps_model.(fnames{ifn}));
        epsilon.(fnames{ifn})(cond) = eps_model.(fnames{ifn})(cond);
    end
end

%% Aerosol layer
%%Calculate top of aerosol layer
for it = 1:length(dt)
    if dt(it)==0.5
        beta_it = beta.t_0_5min;
    else
        eval(['beta_it = ' sprintf('beta.t_%imin;',dt(it))])
    end
    % Look for the top of the aerosol layer in the beta data
    itop_beta_it = ones(size(beta_it));
    aero_top_beta_it = nan(size(beta_it,1),1);
    aero_layer_mask_it = zeros(size(beta_it));
    for ip = 1:size(beta_it,1)
        for jp = cut_h+1:size(beta_it,2)-6
            itop_beta_it(ip,jp) = not((isnan(beta_it(ip,jp)) | beta_it(ip,jp) == 0) &...
                (all(isnan(beta_it(ip,jp:jp+6))) | all(beta_it(ip,jp:jp+6)==0)));
        end
        itop_beta_it(ip,[1:cut_h,find(itop_beta_it(ip,:) == 0,1,'first'):end]) = 0;
        tmp = itop_beta_it(ip,:); tmp(1:cut_h) = 1;
        aero_top_beta_it(ip) = find(tmp == 0,1,'first')-1;
        aero_layer_mask_it(ip,1:aero_top_beta_it(ip)) = 1;
    end
    aero_layer_mask_it(:,1:cut_h) = 0;
    aero_top_beta_it(aero_top_beta_it<cut_h+1) = nan;
    aero_top_beta_it = round(smooth(inpaint_nans(aero_top_beta_it,4),5));
    aero_top_beta_it(all(isnan(beta_it),2)) = nan;
    if dt(it)==.5
        aero_top.t_0_5min = aero_top_beta_it;
        aero_layer_mask.t_0_5min = aero_layer_mask_it;
    else
        eval(['aero_top.t_' sprintf('%imin = aero_top_beta_it;',dt(it))])
        eval(['aero_layer_mask.t_' sprintf('%imin = aero_layer_mask_it;',dt(it))])
    end
end

%% Wind shear
aatimes = struct2cell(t_epsilon);
% dt = [2 3 5];
for kk = 1:length(aatimes)
    [Xr,Yr] = meshgrid(data_vert.range, aatimes{kk});
    % Take day before and after to fill some gaps
    if ~isempty(data_wind_tmrw) && ~isempty(data_wind_yday)
    [Xo,Yo] = meshgrid(data_wind_tday.range,...
        [data_wind_yday.time(:)-24;...
        data_wind_tday.time(:);...
        data_wind_tmrw.time(:)+24]);
    uwind = interp2(Xo,Yo, ...
        [data_wind_yday.uwind; ...
        data_wind_tday.uwind; ...
        data_wind_tmrw.uwind],Xr,Yr);
    vwind = interp2(Xo,Yo, ...
        [data_wind_yday.vwind; ...
        data_wind_tday.vwind; ...
        data_wind_tmrw.vwind],Xr,Yr);
    else
        [Xo,Yo] = meshgrid(data_wind_tday.range,data_wind_tday.time);
        uwind = interp2(Xo,Yo,data_wind_tday.uwind,Xr,Yr);
        vwind = interp2(Xo,Yo,data_wind_tday.vwind,Xr,Yr);
        
    end
    
    speed_tmp = sqrt(uwind.^2+vwind.^2);  
    direc_tmp = rad2deg(atan2(-uwind,-vwind));
    direc_tmp(direc_tmp < 0) = direc_tmp(direc_tmp < 0)+360;
    
    % Calculate wind shear
    win_size = 5;
    shear_dir_tmp = nan(size(uwind));
    shear_vec_tmp = nan(size(uwind));
    for ir = 1:length(aatimes{kk})
        for ic = floor(win_size/2)+1:length(data_wind_tday.range)-floor(win_size/2)
            duwind = uwind(ir,ic+floor(win_size/2))-...
                uwind(ir,ic-floor(win_size/2));
            dvwind = vwind(ir,ic+floor(win_size/2))-...
                vwind(ir,ic-floor(win_size/2));
            dz = data_wind_tday.range(ic + floor(win_size/2))-...
                data_wind_tday.range(ic - floor(win_size/2));
            shear_vec_tmp(ir,ic) = sqrt((duwind).^2 + (dvwind).^2) ./ dz;
            shear_dir_tmp(ir,ic) = 180 * atan2(dvwind,duwind) / pi;
        end
    end
    
    % Try infilling some gaps
    shear_vec_smooth = medianfilter(shear_vec_tmp,[3,5]);
    shear_dir_smooth = medianfilter(shear_dir_tmp,[3,5]);
    [~, counts]      = medianfilter(isfinite(shear_vec_tmp),[3,5]);
    % kernel is [3 3], max count is 9
    shear_vec_smooth(counts < 3) = nan;
    shear_dir_smooth(counts < 3) = nan;
    ind_nm = isnan(shear_vec_tmp) & isfinite(shear_vec_smooth);
    shear_vec_tmp(ind_nm(:,1:10)) = shear_vec_smooth(ind_nm(:,1:10));
    shear_dir_tmp(ind_nm(:,1:10)) = shear_dir_smooth(ind_nm(:,1:10));
    
    shear_vec_tmp(isnan(beta.(fnames{kk}))) = nan;
    shear_dir_tmp(isnan(beta.(fnames{kk}))) = nan;
    if dt(kk)==.5
        shear_vec.t_0_5min = shear_vec_tmp;
        shear_dir.t_0_5min = shear_dir_tmp;
        speed.t_0_5min = speed_tmp;
        direc.t_0_5min = direc_tmp;
    else
        eval(['shear_vec.t_' sprintf('%s',num2str(dt(kk))) 'min = shear_vec_tmp;'])
        eval(['shear_dir.t_' sprintf('%s',num2str(dt(kk))) 'min = shear_dir_tmp;'])
        eval(['speed.t_' sprintf('%s',num2str(dt(kk))) 'min = speed_tmp;'])
        eval(['direc.t_' sprintf('%s',num2str(dt(kk))) 'min = direc_tmp;'])
    end
end

%% Skewness
velo_4skwn = velo.(fnames{1});
% velo_4skwn(real(log10(beta.(fnames{1})))>-5) = nan;
velo_4skwn(velo_4skwn>nanmean(velo_4skwn(:))+6*nanstd(velo_4skwn(:)) | ...
    velo_4skwn < nanmean(velo_4skwn(:))-6*nanstd(velo_4skwn(:))) = nan;
[skewn_tmp,~,~] = windowSlider(velo_4skwn,[120,3],@skewness);
skewn_tmp(isnan(beta.(fnames{1}))) = nan;
t_sk_ref = t_epsilon.(fnames{1});
t_sk_ref(all(isnan(skewn_tmp),2)) = [];
skewn_tmp(all(isnan(skewn_tmp),2),:) = [];

%% Grid data into same reso and clean up
for ifn = 1:length(fnames)
    % skewness
    [Xr,Yr] = meshgrid(data_vert.range,t_epsilon.(fnames{ifn}));
    [Xo,Yo] = meshgrid(data_vert.range,t_sk_ref);
    skewn.(fnames{ifn}) = interp2(Xo,Yo,skewn_tmp,Xr,Yr);
    skewn.(fnames{ifn})(isnan(beta.(fnames{ifn}))) = nan;
    
    % clean cut_h lowest most range gates
    beta.(fnames{ifn})(:,1:cut_h) = nan;
    velo.(fnames{ifn})(:,1:cut_h) = nan;
    epsilon.(fnames{ifn})(:,1:cut_h) = nan;
    skewn.(fnames{ifn})(:,1:cut_h) = nan;
    sigma_w.(fnames{ifn})(:,1:cut_h) = nan;
    shear_vec.(fnames{ifn})(:,1:cut_h) = nan;
    shear_dir.(fnames{ifn})(:,1:cut_h) = nan;
    speed.(fnames{ifn})(:,1:cut_h) = nan;
    direc.(fnames{ifn})(:,1:cut_h) = nan;
end
%% Time out
time_o = t_epsilon;
end
