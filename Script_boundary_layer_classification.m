%%--- --- --- --- --- --- --- --- --- --- --- --- --- --- ---%%
% ADD CORRECT PATHS AND MAKE SURE YOU HAVE THE REQUIRED DATA  %
%%--- --- --- --- --- --- --- --- --- --- --- --- --- --- ---%%

% clc
% clear
% close all

%% GIVE INPUTS
d_type = 'uncorrected';
pol_ch = 'co';
site   = 'hyytiala'; %'arm-graciosa';
save_bkgco = 1; % save bkg corrected data?
save_blc = 1; % save blc product?

%% LOAD DATA
for daten = datenum(YYYY,MM,DD):datenum(YYYY,MM,DD)
    
    % load winds that are needed
    data_wind_tday = combineHaloDBSnVAD(site,d_type,daten);
    data_wind_yday = combineHaloDBSnVAD(site,d_type,daten-1);
    data_wind_tmrw = combineHaloDBSnVAD(site,d_type,daten+1);
    if isempty(data_wind_tday)||isempty(data_wind_yday)||isempty(data_wind_tmrw)
        continue
    else
        
        % Load models winds if exist and speficied for site, otherwise []
        model_wind_tday = loadModel(site,daten,'gdas1');
        model_wind_yday = loadModel(site,daten-1,'gdas1');
        model_wind_tmrw = loadModel(site,daten+1,'gdas1');
        
        % load Halo vertical data and correct bkg if needed
        [data_vert,att] = preprocessHalo(site,d_type,daten,pol_ch,...
            save_bkgco);
        
        %% Locate low level jets
        [LLJ] = detectLLJ(data_wind_tday,daten);
        
        %% CALCULATE DWL QUANTITIES
        dt = [0.5 1 2 3];
        [time_o,beta,skewn,epsilon,epsilon_error,shear_vector,...
            shear_direction,aero_top,velo,velo_error,sigma_w,sigma_w_error,...
            aero_layer_mask,nsamples] = calcWindQuantities(...
            data_vert,data_wind_tday,data_wind_yday,data_wind_tmrw,...
            model_wind_tday,model_wind_yday,model_wind_tmrw,dt);
        fnames = fieldnames(time_o);
        
        %% GET FLUX OR SUNRISE/SET DATA        
%         %%-- TBD --%%
%         path_flux = ['/.../.../' site '/'...
%             datestr(daten,'yyyy') '/23m/' datestr(daten,'yyyy') '_sens_heat_flux.csv'];
%         if exist(path_flux,'file') == 2
%             tmp_flux = importdata(path_flux);
%             t_flux_tmp = datenum(tmp_flux.data(:,1:6));
%             index = t_flux_tmp >= daten & t_flux_tmp < daten+1;
%             flux0 = double(tmp_flux.data(index,7));
%             t_flux_tmp = t_flux_tmp(index);
%             [~,~,~,HH,MM,SS] = datevec(t_flux_tmp);
%             t_flux = HH+MM/60+SS/60/60; % decimal hours
%             sun_rise_set = [];
%             flag_flux = true;
%         else
%         %%-- TBD --%%      
        flux = [];
        time_daten = decimal2daten(time_o.t_1min,daten);
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
        end
        sun_rise_set = suncycle(Location.latitude,Location.longitude,daten);
        
%         %%-- OR --%%
%         time_sa = pvl_maketimestruct(time_daten, 0); % time is already UTC
%         [SunAz, SunEl, AppSunEl, SolarTime] = pvl_ephemeris(time_sa,Location);
%         %%-- OR --%%
        
        flag_flux = false;
        %% SET THRESHOLDS
        th.heatflux = 0;
        th.heatflux_hi = 10;
        th.epsilon = -4; % log scale
        th.windshear = 0.03;
        th.cloud = -5;
        th.vert_velo = 0;
        
        %% TKE CONNECTED WITH?
        for ifn = 1:length(fnames)
            
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
            Eps_log10 = real(log10(epsilon.(fnames{ifn})));
            window_e = 10/(nanmedian(diff(time_o.(fnames{ifn})))*60); % 10 mins
            Eps_log10_med = medianfilter(Eps_log10,[window_e,3]);
            Eps_log10_med(isnan(beta.(fnames{ifn}))) = nan;
            fubarfield.(fnames{ifn}) = associateTKEwith(...
                Eps_log10_med,skewn.(fnames{ifn}),cloudmask,th.epsilon);
            
            %% INTERPOLATE TIME STAMPS
            if flag_flux
                for ifn = 1:length(fnames)
                    flux.(fnames{ifn}) = interp1(t_flux,flux0,time_o.(fnames{ifn}),'pchip');
                end
            else
                for ifn = 1:length(fnames)
                    flux.(fnames{ifn}) = [];
                end
            end
            
            % CREATE BITFIELD
            Eps_med = medianfilter(epsilon.(fnames{ifn}),[window_e,3]);
            Eps_med(isnan(beta.(fnames{ifn}))) = nan;
            
            bitfield.(fnames{ifn}) = createBitfield(th,...
                time_o.(fnames{ifn}),beta.(fnames{ifn}),...
                aero_top.(fnames{ifn}),flux.(fnames{ifn}),Eps_med,...
                shear_vector.(fnames{ifn}),sun_rise_set,fubarfield.(fnames{ifn}));
            
            % BL classification
            [product, product_attribute] = ...
                createBLC(bitfield.(fnames{ifn}),...
                fubarfield.(fnames{ifn}));
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
            dimensions=struct('time',numel(time_o.t_3min),'height',...
                numel(data_vert.range));
            if save_blc
                write_nc_silent(['/.../.../datestr(daten,'yyyymmdd') '_' ...
                    site '_bl_classfication_%s.nc'],dimensions,product,...
                    product_attribute,fnames{ifn})
            end
        end
        
%         %%-- EXAMPLE OF PLOTTING ROUTINE ASSUMING export_fig FUNCTION EXISTS --%%
%
%         fh = figure; fh.Color = 'w'; fh.Units = 'normalized'; fh.Position = [.3 .2 .3 .5];
%         %-- Beta --%
%         sp1 = subplot(4,1,1);
%         pcolor(time_o.(fnames{ifn}),data_vert.range,real(log10(beta.(fnames{ifn})))');
%         shading flat; colormap(sp1,flipud(cbrewer('div','Spectral',99))); caxis([-7 -4]);
%         set(sp1,'Color',[.8 .8 .8],'position',[.1 .8 .62 .16],...
%             'XTick',0:3:24,'XTickLabel',{'0:00','3:00','6:00','9:00',...
%             '12:00','15:00','18:00','21:00','0:00'},'layer','top',...
%             'YTick',0:500:2000)
%         cb = colorbar;
%         cb.Ticks = -7:-4;
%         cb.TickLabels = cellstr(strcat(repmat('10^{',numel(cb.Ticks(:)),1),num2str(cb.Ticks(:)),repmat('}',numel(cb.Ticks(:),1))));
%         cb.Label.String = '(m^{-1} sr^{-1})';
%         cb.Position(1) = .74; cb.Position(3) = .015;
%         axis(sp1,[0 24 0 2000])
%         ylabel('Height (m a.g.l.)')
%         text(0.21,1800,'a)','Color','k')
%         text(0.21,2200,'Attenuated backscatter coefficients','Color','k')
%         %-- Skewness --%
%         sp2 = subplot(4,1,2);
%         pcolor(time_o.(fnames{ifn}),data_vert.range,skewn.(fnames{ifn})');
%         shading flat; colormap(sp2,cmocean('balance',99)); caxis([-2 2]);
%         set(sp2,'Color',[.8 .8 .8],'position',[.1 .57 .62 .16],...
%             'XTick',0:3:24,'XTickLabel',{'0:00','3:00','6:00','9:00',...
%             '12:00','15:00','18:00','21:00','0:00'},'layer','top',...
%             'YTick',0:500:2000)
%         cb = colorbar;
%         cb.Ticks = -2:.5:2;
%         cb.Label.String = '(m^{-1} s^{-1})';
%         cb.Position(1) = .74; cb.Position(3) = .015;
%         axis(sp2,[0 24 0 2000])
%         ylabel('Height (m a.g.l.)')
%         text(0.21,1800,'b)','Color','k')
%         text(0.21,2200,'Radial vertical velocity skewness','Color','k')
%         %-- Epsilon --%
%         sp3 = subplot(4,1,3);
%         pcolor(time_o.(fnames{ifn}),data_vert.range,real(log10(epsilon.(fnames{ifn})))');
%         shading flat; colormap(sp3,morgenstemning(99)); caxis([-6 -1]);
%         set(sp3,'Color',[.8 .8 .8],'position',[.1 .32 .62 .18],...
%             'XTick',0:3:24,'XTickLabel',{'0:00','3:00','6:00','9:00',...
%             '12:00','15:00','18:00','21:00','0:00'},'layer','top',...
%             'YTick',0:500:2000)
%         cb = colorbar;
%         cb.Ticks = -6:-1;
%         cb.TickLabels = cellstr(strcat(repmat('10^{',numel(cb.Ticks(:)),1),num2str(cb.Ticks(:)),repmat('}',numel(cb.Ticks(:),1))));
%         cb.Label.String = '(m^{2} s^{-3})';
%         cb.Position(1) = .74; cb.Position(3) = .015;
%         axis(sp3,[0 24 0 2000])
%         ylabel('Height (m a.g.l.)')
%         text(0.21,1800,'c)','Color','k')
%         text(0.21,2200,'TKE dissipation rate','Color','k')
%         %-- Vector wind shear --%
%         sp4 = subplot(4,1,4);
%         pcolor(time_o.(fnames{ifn}),data_vert.range,shear_vector.(fnames{ifn})');
%         shading flat; colormap(sp4,morgenstemning(99)); caxis([0 .06]);
%         set(sp4,'Color',[.8 .8 .8],'Position',[.1 .07 .62 .18],...
%             'XTick',0:3:24,'XTickLabel',{'0:00','3:00','6:00','9:00',...
%             '12:00','15:00','18:00','21:00','0:00'},'layer','top',...
%             'YTick',0:500:2000)
%         cb = colorbar;
%         cb.Ticks = 0:0.01:0.12;
%         cb.Label.String = '(s^{-1})';
%         cb.Position(1) = .74; cb.Position(3) = .015;
%         axis(sp4,[0 24 0 2000])
%         ylabel('Height (m a.g.l.)')
%         text(0.21,1800,'d)','Color','k')
%         text(0.21,2200,'Vector wind shear','Color','k')
%         xlabel('Time UTC')
%         %-- Save --%
%         pause(2)
%         pathout = '/.../.../';
%         export_fig('-png','-nocrop','-painters','-m2',...
%             sprintf('%s%s_lidar_quantities_%s_%s.png',pathout,...
%             datestr(daten,'yyyy-mm-dd'),site,fnames{ifn}))
%         pause(2)
%         
%         % Prepare for plotting
%         BLclass = double(product.bl_classification);
%         BLclass(BLclass==0) = nan;
%         tmp = BLclass;
%         tmp(isnan(beta.(fnames{ifn})))=nan;
%         BLclass(:,4:end) = tmp(:,4:end);
%         %
%         TKEconnected = double(product.TKE_connected);
%         TKEconnected(TKEconnected==0) = nan;
%         TKEconnected(isnan(beta.(fnames{ifn}))) = nan;
%         
%         fh = figure; fh.Color = 'w'; fh.Units = 'normalized'; fh.Position = [.3 .2 .3 .5];
%         %-- TKE connceted with --%
%         sp1 = subplot(4,1,1);
%         cmap_tkecw = [product_attribute.TKE_connected.legend_key_red;...
%             product_attribute.TKE_connected.legend_key_green;...
%             product_attribute.TKE_connected.legend_key_blue]';
%         pcolor(time_o.(fnames{ifn}),data_vert.range,TKEconnected')
%         shading flat; colormap(sp1,cmap_tkecw); caxis([1 6]);
%         set(sp1,'Color',[.8 .8 .8],'position',[.1 .8 .62 .16],...
%             'XTick',0:3:24,'XTickLabel',{'0:00','3:00','6:00','9:00',...
%             '12:00','15:00','18:00','21:00','0:00'},'layer','top',...
%             'YTick',0:500:2000)
%         cb = colorbar;
%         cb.Ticks = 1.5:5.5;
%         cb.TickLabels = sprintf(product_attribute.TKE_connected.definition); cb.FontSize = 10;
%         cb.Position(1) = .74; cb.Position(3) = .015;
%         axis(sp1,[0 24 0 2000])
%         ylabel('Height (m a.g.l.)')
%         text(0.21,1800,'a)','Color','k')
%         text(0.21,2200,'TKE connected with','Color','k')
%         %-- Boundary layer classification --%
%         sp2 = subplot(4,1,2);
%         cmap_blc = [product_attribute.bl_classification.legend_key_red(1:8);...
%             product_attribute.bl_classification.legend_key_green(1:8);...
%             product_attribute.bl_classification.legend_key_blue(1:8)]';
%         pcolor(time_o.(fnames{ifn}),data_vert.range,BLclass')
%         shading flat; colormap(sp2,cmap_blc); caxis([1 9]);
%         set(sp2,'Color',[.8 .8 .8],'position',[.1 .57 .62 .16],...
%             'XTick',0:3:24,'XTickLabel',{'0:00','3:00','6:00','9:00',...
%             '12:00','15:00','18:00','21:00','0:00'},'layer','top',...
%             'YTick',0:500:2000)
%         cb = colorbar;
%         cb.Ticks = 1.5:8.5;
%         cb.TickLabels = sprintf(product_attribute.bl_classification.definition); cb.FontSize = 10;
%         cb.Position(1) = .74; cb.Position(3) = .015;
%         axis(sp2,[0 24 0 2000])
%         ylabel('Height (m a.g.l.)')
%         text(0.21,1800,'b)','Color','k')
%         text(0.21,2200,'Boundary layer classification','Color','k')
%         %-- Wind speed and LLJ --%
%         sp3 = subplot(4,1,3); hold on
%         time_LLJs = (LLJ.time_n-daten)*24;
%         LLJs_h = LLJ.height; LLJs_h(time_LLJs>sun_rise_set(1) & time_LLJs<sun_rise_set(2)) = nan;
%         labels = LLJ.label; labels(isnan(LLJs_h)) = nan;
%         uniqs = unique(labels(~isnan(labels)));
%         pcolor(data_wind_tday.time,data_wind_tday.range,data_wind_tday.wind_direction');
%         for i = 1:numel(uniqs)
%             plot(time_LLJs(labels==uniqs(i)),LLJs_h(labels==uniqs(i)),'k-','LineWidth',.5)
%         end
%         shading flat; colormap(sp3,colorcet('C8')); caxis([0 360]);
%         set(sp3,'Color',[.8 .8 .8],'position',[.1 .32 .62 .18],...
%             'XTick',0:3:24,'XTickLabel',{'0:00','3:00','6:00','9:00',...
%             '12:00','15:00','18:00','21:00','0:00'},'layer','top',...
%             'YTick',0:500:2000,'Box','on')
%         cb = colorbar;
%         cb.Ticks = 0:45:360;
%         cb.Label.String = '(deg)';
%         cb.Position(1) = .74; cb.Position(3) = .015;
%         axis(sp3,[0 24 0 2000])
%         ylabel('Height (m a.g.l.)')
%         text(0.21,1800,'c)','Color','k')
%         text(0.21,2200,'Wind direction and LLJ','Color','k')
%         hold off
%         %-- Wind speed --%
%         sp4 = subplot(4,1,4); hold on
%         pcolor(data_wind_tday.time,data_wind_tday.range,data_wind_tday.wind_speed');
%         for i = 1:numel(uniqs)
%             plot(time_LLJs(labels==uniqs(i)),LLJs_h(labels==uniqs(i)),'k-','LineWidth',.5)
%             plot(time_LLJs(labels==uniqs(i)),LLJs_h(labels==uniqs(i)),'k.','MarkerSize',10)
%         end
%         shading flat; colormap(sp4,cmocean('thermal',99)); caxis([0 20]);
%         set(sp4,'Color',[.8 .8 .8],'Position',[.1 .07 .62 .18],...
%             'XTick',0:3:24,'XTickLabel',{'0:00','3:00','6:00','9:00',...
%             '12:00','15:00','18:00','21:00','0:00'},'layer','top',...
%             'YTick',0:500:2000,'Box','on')
%         cb = colorbar;
%         cb.Ticks = 0:2:20;
%         cb.Label.String = '(m s^{-1})';
%         cb.Position(1) = .74; cb.Position(3) = .015;
%         axis(sp4,[0 24 0 2000])
%         ylabel('Height (m a.g.l.)')
%         text(0.21,1800,'d)','Color','k')
%         text(0.21,2200,'Wind speed and LLJ','Color','k')
%         xlabel('Time UTC')
%         hold off
%         %-- Save --%
%         pause(2)
%         pathout = '/.../.../';
%         export_fig('-png','-nocrop','-painters','-m2',...
%             sprintf('%s%s_BL_classification_and_LLJ_%s_%s.png',pathout,...
%             datestr(daten,'yyyy-mm-dd'),site,fnames{ifn}))
%         pause(2)

    end
end

