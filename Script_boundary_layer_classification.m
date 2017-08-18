%%--- --- --- --- --- --- --- --- --- --- --- --- --- --- ---%%
% ADD CORRECT PATHS AND MAKE SURE YOU HAVE THE REQUIRED DATA  %
%%--- --- --- --- --- --- --- --- --- --- --- --- --- --- ---%%

addpath(genpath('/home/tmarke/Dokumente/PhD/programs/matlab/'))
addpath(genpath('/home/tmarke/Dokumente/PhD/programs/matlab/bl_class/function/export_fig/'))
addpath(genpath('/home/hatpro/matlab/cloudnet/'))

DefFontSize = 12;
set(0,'DefaultAxesLineWidth',1);
set(0,'DefaultLineLineWidth',3);
set(0,'DefaultAxesFontSize',DefFontSize);
set(0,'DefaultAxesFontWeight','Normal');
set(0,'DefaultTextFontSize',DefFontSize);
set(0,'DefaultTextFontWeight','Normal');
set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')

clc
clear
close all

%% GIVE INPUTS
d_type = 'uncorrected';
pol_ch = 'co';
site   = 'arm-oliktok'; %'arm-graciosa';
save_bkgco = 1; % save bkg corrected data?
save_blc = 1; % save blc product?
plot_vis = 'off'; % show blc plot?
save_plot = 0; % save blc plot

% unreliable range in m
cut = 100;

%% LOAD DATA
for daten = datenum(yyyy,mm,dd):datenum(yyyy,mm,dd)
    
    % load winds that are needed
    data_wind_tday = combineHaloDBSnVAD(site,daten,cut);
    data_wind_yday = combineHaloDBSnVAD(site,daten-1,cut);
    data_wind_tmrw = combineHaloDBSnVAD(site,daten+1,cut);
    if isempty(data_wind_tday)
        continue
    else
        
        if ~isempty(data_wind_yday) && size(data_wind_tday.range,1) ~= size(data_wind_yday.range,1)
            data_wind_yday = [];
        end
        if ~isempty(data_wind_tmrw) && size(data_wind_tday.range,1) ~= size(data_wind_tmrw.range,1)
            data_wind_tmrw = [];
        end
        
        % Load models winds if exist and speficied for site, otherwise []
        model_wind_tday = loadModel(site,daten,'ecmwf'); %dwd-lmk-9-11
        model_wind_yday = loadModel(site,daten-1,'ecmwf');
        model_wind_tmrw = loadModel(site,daten+1,'ecmwf');
        if ~isempty(model_wind_yday) && ~isempty(model_wind_tmrw) && ~strcmp(site,'arm-oliktok')
            model_wind_yday.time(end)       = [];
            model_wind_tmrw.time(1)         = [];
            model_wind_yday.uwind(end,:)    = [];
            model_wind_tmrw.uwind(1,:)      = [];
            model_wind_yday.vwind(end,:)    = [];
            model_wind_tmrw.vwind(1,:)      = [];
        end
        
        % load Halo vertical data and correct bkg if needed
        [data_vert,att] = preprocessHalo(site,d_type,daten,pol_ch,...
            save_bkgco,cut);
        
        %% Locate low level jets
        [LLJ] = detectLLJ(data_wind_tday,daten);
        
        %% CALCULATE DWL QUANTITIES
        dt = [0.5 1 2 3];
        [time_o,beta,skewn,epsilon,epsilon_error,shear_vector,...
            shear_direction,aero_top,velo,velo_error,sigma_w,sigma_w_error,...
            aero_layer_mask,nsamples,speed,direc,signal] = calcWindQuantities(...
            data_vert,data_wind_tday,data_wind_yday,data_wind_tmrw,...
            model_wind_tday,model_wind_yday,model_wind_tmrw,dt,site);
        fnames = fieldnames(time_o);
        
        %% GET FLUX OR SUNRISE/SET DATA
        
        switch site
            case 'arm-oliktok'
                path_flux = ['/data/hatpro/jue/cloudnet/juelich/calibrated/dopplerlidar/oli-ecor/'...
                    'oli30ecorM1.b1.' datestr(daten,'yyyymmdd') '.000000.cdf'];
            otherwise
                path_flux = [];
        end
            
        if exist(path_flux,'file') == 2
            flux = load_nc_struct(path_flux);
            t_flux = flux.time./3600;
            flux0 = flux.h;
            flux0(flux.qc_h > 0) = NaN;

        % INTERPOLATE TIME STAMPS
            for ifn = 1:length(fnames)
                flux.(fnames{ifn}) = interp1(t_flux,flux0,time_o.(fnames{ifn}),'pchip');
            end

        else

            for ifn = 1:length(fnames)
                flux.(fnames{ifn}) = [];
            end
        end
                    
        switch site
            case 'hyytiala'
                Location.latitude = 61.845225;
                Location.longitude = 24.287037;
                Location.altitude = 163; % m
                UTCoffset = 2;
            case 'juelich'
                Location.latitude = 50.908547;
                Location.longitude = 6.413536;
                Location.altitude = 111; % m
                UTCoffset = 1;
            case 'limassol'
                Location.latitude = 34.1676998;
                Location.longitude = 33.037998;
                Location.altitude = 5; % m
                UTCoffset = 2;
            case 'arm-oliktok'
                Location.latitude = 70.494856;
                Location.longitude = -149.886470;
                Location.altitude = 2; % m
                UTCoffset = -9;
        end
            sun_rise_set = suncycle(Location.latitude,Location.longitude,daten);

%         %%-- OR --%%
%         time_daten = decimal2daten(time_o.t_1min,daten);
%         time_sa = pvl_maketimestruct(time_daten, 0); % time is already UTC
%         [SunAz, SunEl, AppSunEl, SolarTime] = pvl_ephemeris(time_sa,Location);
%         %%-- OR --%%
        
        %% SET THRESHOLDS
        th.heatflux = 0;
        th.heatflux_hi = 10;
        th.epsilon = -4; % log scale
        th.windshear = 0.03;
        th.cloud = -5;
        th.vert_velo = 0;
        
        %% TKE CONNECTED WITH?
        for ifn = 4:4%length(fnames)
            
            cloudmask = zeros(size(beta.(fnames{ifn})));
            cloudmask(real(log10(beta.(fnames{ifn}))) > th.cloud) = 1;
            % Dilate
            shift_m1 = [cloudmask(:,2:end) zeros(size(cloudmask,1),1)];
            shift_p1 = [zeros(size(cloudmask,1),1) cloudmask(:,1:end-1)];
            cloudmask_dilated = cloudmask + shift_m1 + shift_p1;
            cloudmask_dilated(cloudmask_dilated<0) = 1;
            cloudmask_dilated(aero_layer_mask.(fnames{ifn})==0) = 0;
            cloudmask = logical(cloudmask_dilated);
            % Smooth TKE
            window_e = 10/(nanmedian(diff(time_o.(fnames{ifn})))*60); % 10 mins
            Eps_med = medianfilter(epsilon.(fnames{ifn}),[window_e,3]);
            Eps_med(isnan(epsilon.(fnames{ifn}))) = nan;
            fubarfield.(fnames{ifn}) = associateTKEwith(...
                Eps_med,skewn.(fnames{ifn}),cloudmask,th.epsilon,data_vert.cut);
            
            bitfield.(fnames{ifn}) = createBitfield(th,...
                time_o.(fnames{ifn}),beta.(fnames{ifn}),...
                aero_top.(fnames{ifn}),flux.(fnames{ifn}),Eps_med,...
                shear_vector.(fnames{ifn}),sun_rise_set,fubarfield.(fnames{ifn}));
            
            % BL classification
            [product, product_attribute] = ...
                createBLC(bitfield.(fnames{ifn}),...
                fubarfield.(fnames{ifn}),data_vert.cut);
            product.bl_classification(product.bl_classification == 3 & ...
                isnan(epsilon.(fnames{ifn}))) = 0;
            product.TKE_connected(product.bl_classification == 0) = 0;
            product.time = time_o.(fnames{ifn})(:);
            product.height = data_vert.range(:);
            product_attribute.time = att.time;
            product_attribute.time.dimensions = {'time'};
            product_attribute.height = att.range;
            product_attribute.height.units = 'm';
            product_attribute.height.dimensions = {'height'};
            product_attribute.global = att.global;
            dimensions=struct('time',numel(time_o.(fnames{ifn})(:)),'height',...
                numel(data_vert.range),'sun',numel(sun_rise_set));
            
            product.speed = speed.(fnames{ifn});
            product.dir = direc.(fnames{ifn});
            product.beta = real(log10(beta.(fnames{ifn})));
            product.shear = shear_vector.(fnames{ifn});
            product.eps = real(log10(epsilon.(fnames{ifn})));
            product.skew = skewn.(fnames{ifn});
            product.velo = velo.(fnames{ifn});
            product.vel_err = real(log10(velo_error.(fnames{ifn})));
            product.eps_err = epsilon_error.(fnames{ifn});
            product.snr = signal.(fnames{ifn});
            
            if ~isempty(flux.(fnames{ifn}))
                product.flux = flux.(fnames{ifn});
                product_attribute.flux.dimensions = {'time'};
                product_attribute.flux.units = 'W m^{-2}';
                product_attribute.flux.long_name = 'eddy covariance sensible heat flux';
            end
            
            product.sun = sun_rise_set(:);
            product_attribute.sun.dimensions = {'sun'};
            product_attribute.sun.units = 'UTC';
            product_attribute.sun.long_name = 'sun rise and sun set times';

            product_attribute.speed.dimensions = {'time','height'};
            product_attribute.speed.units = 'm s^{-1}';
            product_attribute.speed.long_name = 'horizontal wind speed';

            product_attribute.dir.dimensions = {'time','height'};
            product_attribute.dir.units = 'degree';
            product_attribute.dir.long_name = 'wind direction';

            product_attribute.beta.dimensions = {'time','height'};
            product_attribute.beta.units = 'm^{-1} sr^{-1}';
            product_attribute.beta.long_name = 'attenuated backscatter coefficient background corrected';

            product_attribute.shear.dimensions = {'time','height'};
            product_attribute.shear.units = 's^{-1}';
            product_attribute.shear.long_name = 'vector wind shear over 5 range gates';

            product_attribute.eps.dimensions = {'time','height'};
            product_attribute.eps.units = 'm^{2} s^{-3}';
            product_attribute.eps.long_name = 'TKE dissipation rate';

            product_attribute.skew.dimensions = {'time','height'};
            product_attribute.skew.units = 'm^{-1} s^{-1}';
            product_attribute.skew.long_name = 'skewness averaged over 120 min and 3 range gates';

            product_attribute.velo.dimensions = {'time','height'};
            product_attribute.velo.units = 'm s^{-1}';
            product_attribute.velo.long_name = 'vertical velocity';

            product_attribute.vel_err.dimensions = {'time','height'};
            product_attribute.vel_err.units = 'm s^{-1}';
            product_attribute.vel_err.long_name = 'error in vertical velocity';

            product_attribute.eps_err.dimensions = {'time','height'};
            product_attribute.eps_err.units = '';
            product_attribute.eps_err.long_name = 'fractional error in TKE dissipation rate';

            product_attribute.snr.dimensions = {'time','height'};
            product_attribute.snr.units = '';
            product_attribute.snr.long_name = 'signal to noise ratio';
            
            if save_blc
                write_nc_struct(sprintf('/data/hatpro/jue/cloudnet/juelich/products/bl-classification/%s/%s_bl_classification_%s_%s.nc', datestr(daten,'yyyy'), datestr(daten,'yyyymmdd'), site, fnames{ifn}),...
                    dimensions, product, product_attribute)
            end
%         end
        %%
        %%-- EXAMPLE OF PLOTTING ROUTINE ASSUMING export_fig FUNCTION EXISTS --%%
        fh = figure('Visible',plot_vis); fh.Color = 'w'; fh.Units = 'normalized'; fh.Position = [.3 0 .4 1];
        %-- Beta --%
        sp1 = subplot(4,1,1);
        pcolor(time_o.(fnames{ifn}),data_vert.range,real(log10(beta.(fnames{ifn})))');
        shading flat; colormap(sp1,flipud(cbrewer('div','Spectral',99))); caxis([-7 -4]);
        set(sp1,'Color',[.8 .8 .8],'position',[.1 .8 .62 .16],...
            'XTick',0:3:24,'XTickLabel',{'0:00','3:00','6:00','9:00',...
            '12:00','15:00','18:00','21:00','0:00'},'layer','top',...
            'YTick',0:500:2000)
        cb = colorbar;
        cb.Ticks = -7:-4;
        cb.TickLabels = cellstr(strcat(repmat('10^{',numel(cb.Ticks(:)),1),num2str(cb.Ticks(:)),...
            repmat('}',numel(cb.Ticks(:),1))));
        cb.Label.String = '(m^{-1} sr^{-1})';
        cb.Position(1) = .74; cb.Position(3) = .015;
        axis(sp1,[0 24 0 2000])
        ylabel('Height (m a.g.l.)')
        text(0.21,1800,'a)','Color','k')
        text(0.21,2200,'Attenuated backscatter coefficients','Color','k')
        %-- Skewness --%
        sp2 = subplot(4,1,2);
        pcolor(time_o.(fnames{ifn}),data_vert.range,skewn.(fnames{ifn})');
        shading flat; colormap(sp2,cmocean('balance',99)); caxis([-2 2]);
        set(sp2,'Color',[.8 .8 .8],'position',[.1 .57 .62 .16],...
            'XTick',0:3:24,'XTickLabel',{'0:00','3:00','6:00','9:00',...
            '12:00','15:00','18:00','21:00','0:00'},'layer','top',...
            'YTick',0:500:2000)
        cb = colorbar;
        cb.Ticks = -2:.5:2;
        cb.Label.String = '(m^{-1} s^{-1})';
        cb.Position(1) = .74; cb.Position(3) = .015;
        axis(sp2,[0 24 0 2000])
        ylabel('Height (m a.g.l.)')
        text(0.21,1800,'b)','Color','k')
        text(0.21,2200,'Radial vertical velocity skewness','Color','k')
        %-- Epsilon --%
        sp3 = subplot(4,1,3);
        pcolor(time_o.(fnames{ifn}),data_vert.range,real(log10(epsilon.(fnames{ifn})))');
        shading flat; colormap(sp3,morgenstemning(99)); caxis([-6 -1]);
        set(sp3,'Color',[.8 .8 .8],'position',[.1 .32 .62 .18],...
            'XTick',0:3:24,'XTickLabel',{'0:00','3:00','6:00','9:00',...
            '12:00','15:00','18:00','21:00','0:00'},'layer','top',...
            'YTick',0:500:2000)
        cb = colorbar;
        cb.Ticks = -6:-1;
        cb.TickLabels = cellstr(strcat(repmat('10^{',numel(cb.Ticks(:)),1),num2str(cb.Ticks(:)),...
            repmat('}',numel(cb.Ticks(:),1))));
        cb.Label.String = '(m^{2} s^{-3})';
        cb.Position(1) = .74; cb.Position(3) = .015;
        axis(sp3,[0 24 0 2000])
        ylabel('Height (m a.g.l.)')
        text(0.21,1800,'c)','Color','k')
        text(0.21,2200,'TKE dissipation rate','Color','k')
        %-- Vector wind shear --%
        sp4 = subplot(4,1,4);
        pcolor(time_o.(fnames{ifn}),data_vert.range,shear_vector.(fnames{ifn})');
        shading flat; colormap(sp4,morgenstemning(99)); caxis([0 .06]);
        set(sp4,'Color',[.8 .8 .8],'Position',[.1 .07 .62 .18],...
            'XTick',0:3:24,'XTickLabel',{'0:00','3:00','6:00','9:00',...
            '12:00','15:00','18:00','21:00','0:00'},'layer','top',...
            'YTick',0:500:2000)
        cb = colorbar;
        cb.Ticks = 0:0.01:0.12;
        cb.Label.String = '(s^{-1})';
        cb.Position(1) = .74; cb.Position(3) = .015;
        axis(sp4,[0 24 0 2000])
        ylabel('Height (m a.g.l.)')
        text(0.21,1800,'d)','Color','k')
        text(0.21,2200,'Vector wind shear','Color','k')
        xlabel('Time UTC')
        %-- Save --%
        if save_plot == 1
            pause(2)
            pathout = '/data/hatpro/jue/cloudnet/juelich/products/bl-classification/oli/';
            export_fig('-png','-nocrop','-painters','-m2',...
                sprintf('%s%s_lidar_quantities_%s_%s.png',pathout,...
                datestr(daten,'yyyy-mm-dd'),site,fnames{ifn}))
            pause(2)
        end

        % Prepare for plotting
        BLclass = double(product.bl_classification);
        BLclass(BLclass==0) = nan;
        tmp = BLclass;
        tmp(isnan(beta.(fnames{ifn})))=nan;
        BLclass(:,data_vert.cut+1:end) = tmp(:,data_vert.cut+1:end);
        %
        TKEconnected = double(product.TKE_connected);
        TKEconnected(TKEconnected==0) = nan;
        TKEconnected(isnan(beta.(fnames{ifn}))) = nan;

        fh = figure('Visible',plot_vis); fh.Color = 'w'; fh.Units = 'normalized'; fh.Position = [.3 0 .4 1];
        %-- TKE connceted with --%
        sp1 = subplot(4,1,1);
        cmap_tkecw = [product_attribute.TKE_connected.legend_key_red;...
            product_attribute.TKE_connected.legend_key_green;...
            product_attribute.TKE_connected.legend_key_blue]';
        pcolor(time_o.(fnames{ifn}),data_vert.range,TKEconnected')
        shading flat; colormap(sp1,cmap_tkecw); caxis([1 6]);
        set(sp1,'Color',[.8 .8 .8],'position',[.1 .8 .62 .16],...
            'XTick',0:3:24,'XTickLabel',{'0:00','3:00','6:00','9:00',...
            '12:00','15:00','18:00','21:00','0:00'},'layer','top',...
            'YTick',0:500:2000)
        cb = colorbar;
        cb.Ticks = 1.5:5.5;
        cb.TickLabels = sprintf(product_attribute.TKE_connected.definition); cb.FontSize = 10;
        cb.Position(1) = .74; cb.Position(3) = .015;
        axis(sp1,[0 24 0 2000])
        ylabel('Height (m a.g.l.)')
        text(0.21,1800,'a)','Color','k')
        text(0.21,2200,'TKE connected with','Color','k')
        %-- Boundary layer classification --%
        sp2 = subplot(4,1,2);
        cmap_blc = [product_attribute.bl_classification.legend_key_red(1:8);...
            product_attribute.bl_classification.legend_key_green(1:8);...
            product_attribute.bl_classification.legend_key_blue(1:8)]';
        pcolor(time_o.(fnames{ifn}),data_vert.range,BLclass')
        shading flat; colormap(sp2,cmap_blc); caxis([1 9]);
        set(sp2,'Color',[.8 .8 .8],'position',[.1 .57 .62 .16],...
            'XTick',0:3:24,'XTickLabel',{'0:00','3:00','6:00','9:00',...
            '12:00','15:00','18:00','21:00','0:00'},'layer','top',...
            'YTick',0:500:2000)
        cb = colorbar;
        cb.Ticks = 1.5:8.5;
        cb.TickLabels = sprintf(product_attribute.bl_classification.definition); cb.FontSize = 10;
        cb.Position(1) = .74; cb.Position(3) = .015;
        axis(sp2,[0 24 0 2000])
        ylabel('Height (m a.g.l.)')
        text(0.21,1800,'b)','Color','k')
        text(0.21,2200,'Boundary layer classification','Color','k')
        %-- Wind speed and LLJ --%
        sp3 = subplot(4,1,3); hold on
        time_LLJs = (LLJ.time_n-daten)*24;
        LLJs_h = LLJ.height; %LLJs_h(time_LLJs>sun_rise_set(1) & time_LLJs<sun_rise_set(2)) = nan;
        labels = LLJ.label; labels(isnan(LLJs_h)) = nan;
        uniqs = unique(labels(~isnan(labels)));
        pcolor(data_wind_tday.time,data_wind_tday.range,data_wind_tday.wind_direction');
        for i = 1:numel(uniqs)
            plot(time_LLJs(labels==uniqs(i)),LLJs_h(labels==uniqs(i)),'k-','LineWidth',.5)
            plot(time_LLJs(labels==uniqs(i)),LLJs_h(labels==uniqs(i)),'k.','MarkerSize',10)
        end
        shading flat; colormap(sp3,colorcet('C8')); caxis([0 360]);
        set(sp3,'Color',[.8 .8 .8],'position',[.1 .32 .62 .18],...
            'XTick',0:3:24,'XTickLabel',{'0:00','3:00','6:00','9:00',...
            '12:00','15:00','18:00','21:00','0:00'},'layer','top',...
            'YTick',0:500:2000,'Box','on')
        cb = colorbar;
        cb.Ticks = 0:45:360;
        cb.Label.String = '(deg)';
        cb.Position(1) = .74; cb.Position(3) = .015;
        axis(sp3,[0 24 0 2000])
        ylabel('Height (m a.g.l.)')
        text(0.21,1800,'c)','Color','k')
        text(0.21,2200,'Wind direction and LLJ','Color','k')
        hold off
        %-- Wind speed --%
        sp4 = subplot(4,1,4); hold on
        pcolor(data_wind_tday.time,data_wind_tday.range,data_wind_tday.wind_speed');
        for i = 1:numel(uniqs)
            plot(time_LLJs(labels==uniqs(i)),LLJs_h(labels==uniqs(i)),'k-','LineWidth',.5)
            plot(time_LLJs(labels==uniqs(i)),LLJs_h(labels==uniqs(i)),'k.','MarkerSize',10)
        end
        shading flat; colormap(sp4,cmocean('thermal',99)); caxis([0 30]);
        set(sp4,'Color',[.8 .8 .8],'Position',[.1 .07 .62 .18],...
            'XTick',0:3:24,'XTickLabel',{'0:00','3:00','6:00','9:00',...
            '12:00','15:00','18:00','21:00','0:00'},'layer','top',...
            'YTick',0:500:2000,'Box','on')
        cb = colorbar;
        cb.Ticks = 0:5:30;
        cb.Label.String = '(m s^{-1})';
        cb.Position(1) = .74; cb.Position(3) = .015;
        axis(sp4,[0 24 0 2000])
        ylabel('Height (m a.g.l.)')
        text(0.21,1800,'d)','Color','k')
        text(0.21,2200,'Wind speed and LLJ','Color','k')
        xlabel('Time UTC')
        hold off
        %-- Save --%
        if save_plot == 1
            pause(2)
            pathout = '/data/hatpro/jue/cloudnet/juelich/products/bl-classification/oli/';
            export_fig('-png','-nocrop','-painters','-m2',...
                sprintf('%s%s_BL_classification_and_LLJ_%s_%s.png',pathout,...
                datestr(daten,'yyyy-mm-dd'),site,fnames{ifn}))
            pause(2)
        end
        end
    end
%     clearvars -except d_type pol_ch site save_bkgco save_blc plot_vis save_plot cut_h cut_v 
end

