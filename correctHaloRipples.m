function [snr1] = correctHaloRipples(snr0,t_snr,site,daten)
%correct_HALO_bkg_ripples correct the ripples in the HALO bakcground, see
%Vakkari et al. (201?)

% Calculate background from the bkg files
bkg_path = ['/.../.../' site ...
    '/uncorrected/' datestr(daten,'yyyy') '/background/'];

files_bkg = dir(bkg_path); % vertical
fnames_bkg = {files_bkg(~[files_bkg.isdir]).name}'; % name list
ifile_bkg = strfind(fnames_bkg, datestr(daten,'ddmmyy')); % find dates
ifile_bkg = find(not(cellfun('isempty', ifile_bkg)), 1);
if ~isempty(ifile_bkg)
% if exist(path_bkg,'dir') == 7 % if bkg files exist
    [b_file, b_fit, bkg_times] = calculateBKG(bkg_path,daten,size(snr0,2));
    b_file(:,1:3) = nan; b_fit(:,1:3) = nan;
    b_file(all(isnan(b_file),2),:) = []; 
    b_fit(all(isnan(b_fit),2),:) = [];
    
    % Find steps
    [time_snr_dnum] = decimal2daten(t_snr,daten);
    step_locations = nan(1,length(bkg_times)-1);
    for i = 2:length(bkg_times)
        if ~isempty(find(time_snr_dnum>bkg_times(i),1,'first'))
            step_locations(i-1) = find(time_snr_dnum>bkg_times(i),1,'first');
        end
    end
    
    % Correct ripples
    istep = 1;
    snr1 = nan(size(snr0));
    for i = 1:size(snr0,1)
        if istep < size(b_file,1)
            if i >= step_locations(istep)
                istep = istep + 1;
            end
        end
        b_snr = b_file(istep,:)./b_fit(istep,:);
        snr1(i,:) = snr0(i,:) + b_snr - 1;
    end
else 
    snr1 = snr0;
end
end

