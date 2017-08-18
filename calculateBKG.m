function [bkg_out, fit_out, bkg_times] = calculateBKG(bkg_path,daten,n_range_gates)

%% find co and cross background files

dates=datestr(daten,'ddmmyy');

filess=dir([bkg_path 'Background_' dates '*.txt*']);

bkg_times=nan(length(filess),1); % col1: time
for i=1:length(filess)
    b_daten=datenum([filess(i).name(12:17) filess(i).name(19:24)],'ddmmyyHHMMSS');
    bkg_times(i,:)=b_daten;
end

%% read in backgrounds
bkg=nan(length(filess),n_range_gates);
for i=1:length(filess)
    
    if ~isempty(strfind(filess(i).name,'.gz'))
        gunzip([bkg_path filess(i).name],['/data/hatpro/jue/cloudnet/juelich/calibrated/halo-doppler-lidar/'...
            datestr(daten,'yyyy') '/'])
        files=dir(['/data/hatpro/jue/cloudnet/juelich/calibrated/halo-doppler-lidar/'...
            datestr(daten,'yyyy') '/' 'Background_' dates '*.txt']);
        fn=['/data/hatpro/jue/cloudnet/juelich/calibrated/halo-doppler-lidar/'...
            datestr(daten,'yyyy') '/' files.name];
        fid=fopen(fn,'r');
        bk=fscanf(fid,'%s');
        fclose(fid);
        delete(['/data/hatpro/jue/cloudnet/juelich/calibrated/halo-doppler-lidar/'...
            datestr(daten,'yyyy') '/' files.name])
    else
    
        fn=[bkg_path filess(i).name];
        fid=fopen(fn,'r');
        bk=fscanf(fid,'%s');
        fclose(fid);
    
    end
    
    dot_i=find(bk=='.');
    end_i=[1;(dot_i+7)'];
%     disp([fn ' number of points: ' num2str(sum(bk=='.')) ])
    bi=1;
    %         for ii=2:length(end_i)
    for ii=2:n_range_gates+1
        bkg(i,bi)=str2num(bk(end_i(ii-1):end_i(ii)-1));
        bi=bi+1;
    end
    
    %         field_width=median(diff(find(bk=='.')));
    %         disp([fn ' number of points: ' num2str(sum(bk=='.')) ' field width: ' num2str(field_width)])
    %         bi=1;
    %         for ii=1:field_width:(field_width*n_range_gates)
    %             temppi1(1,bi)=str2num(bk(ii:min((ii+field_width-1),length(bk))));
    %             bi=bi+1;
    %         end
    % return
end

%% gapfilling

if isempty(bkg) % no data for a day
    bkg=nan(24,n_range_gates+1);
    bkg(:,1)=(0:23)/24+daten;
else
    
    bkg=[bkg_times bkg];
    if str2num(datestr(bkg(1,1),'HH'))~=0
        bkg=[[floor(bkg_times(1));bkg(:,1)] [bkg(1,2:end)*nan;bkg(:,2:end)]];
    end
    if str2num(datestr(bkg_times(end),'HH'))~=23
        bkg=[[bkg(:,1); daten+23/24] [bkg(:,2:end);bkg(1,2:end)*nan]];
    end
    
    dd=(diff(bkg(:,1)-daten))*24;
    for i=1:length(dd)
        if dd(i)>1.5
            new_t=((bkg(i,1)+1/24):1/24:(bkg(i,1)+floor(dd(i))/24))';
            bkg=[[bkg(:,1); new_t] [bkg(:,2:end);repmat(bkg(1,2:end)*nan,length(new_t),1)]];
        end
    end
    
    if any(dd>1.5)
        bkg=sortrows(bkg,1);
    end
end
bkg_raw=bkg;
clear bkg;
%% to SNR
% bkg_snr=bkg_raw;
% bkg_snr(:,2:end)=nan;
fits_1=bkg_raw(:,1:3)*nan; % p(1) p(2) rmse
fits_2=bkg_raw(:,1:4)*nan; % p(1) p(2) p(3) rmse

fit_out = nan(length(bkg_raw(:,1)),n_range_gates);
bkg_out = nan(length(bkg_raw(:,1)),n_range_gates);
for i=1:length(bkg_raw(:,1));
    b_temp=bkg_raw(i,2:end);
    if ~isnan(b_temp(1)) % fit 1 and 2 order polynomial
        %             b_temp(b_temp<5e6)=b_temp(b_temp<5e6)+1e7;
        fitti_1=polyfit((4:n_range_gates),b_temp(4:n_range_gates),1);
        bkg_fitted_1=(1:n_range_gates)*fitti_1(1)+fitti_1(2);
        fitti_2=polyfit((4:n_range_gates),b_temp(4:n_range_gates),2);
        bkg_fitted_2=((1:n_range_gates).^2)*fitti_2(1)+(1:n_range_gates)*fitti_2(2)+fitti_2(3);
        
        rmse_1=sqrt(mean((b_temp(4:n_range_gates)-bkg_fitted_1(4:n_range_gates)).^2,2));
        rmse_2=sqrt(mean((b_temp(4:n_range_gates)-bkg_fitted_2(4:n_range_gates)).^2,2));
        
%         fits_1(i,:)=[fitti_1 rmse_1];
%         fits_2(i,:)=[fitti_2 rmse_2];
        
        % with plotting
        %             figure(1)
        %             plot(4:n_range_gates,b_temp(4:n_range_gates),'.b',4:n_range_gates,bkg_fitted_1(4:n_range_gates),'-r',4:n_range_gates,bkg_fitted_2(4:n_range_gates),'-g')
                    if rmse_2<(0.9*rmse_1)
                        fit_out(i,:) = bkg_fitted_2;
        %                 title([num2str(rmse_1) ' ' num2str(rmse_2) ' 2nd order poly chosen!'])
%                         bkg_snr(i,2:end)=(b_temp./bkg_fitted_2 -1); % SNR, not SNR+1!
        %                 disp([datestr(bkg_raw(i,1)) ' 2nd order poly chosen!'])
                    else
                        fit_out(i,:) = bkg_fitted_1;
        %                 title([num2str(rmse_1) ' ' num2str(rmse_2)])
%                         bkg_snr(i,2:end)=(b_temp./bkg_fitted_1 -1); % SNR, not SNR+1!
                    end
        %             drawnow
        %             pause(0.5)
        
        bkg_out(i,:) = b_temp;
    end
end
%% save background for one day
%     save([out_path 'Background_' datestr(daten,'yyyymmdd') '_2.mat' ],'bkg_raw','bkg_snr','fits_1','fits_2');
% end

end




