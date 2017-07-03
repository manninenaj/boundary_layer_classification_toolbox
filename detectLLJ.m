function [LLJout] = detectLLJ(data,daten)
%LLJ detection finds low level jets in the data, see Tuononen et al. (2017)

%% Combine winds


%% Find min and max
%take only data below 1510m
ind = find(data.range<1510);
LLJ.wind_speed = data.wind_speed(:,ind);
LLJ.wind_direction = data.wind_direction(:,ind);
LLJ.uwind = data.uwind(:,ind);
LLJ.vwind = data.vwind(:,ind);
LLJ.range = data.range(ind);
LLJ.time = data.time;
%LLJ.nval  = data.nval(:,ind);

%how many min and max values
nvalues = 10;

height_min = zeros(length(LLJ.time),nvalues) * NaN;
height_max = height_min;
ws_max = height_min;
ws_min = height_min;
dir_max = height_min;
dir_min = height_min;
ind_max = height_min;
ind_min = height_min;


%loop over all time steps
for iray = 1:length(LLJ.time)
    
    index = find(isfinite(LLJ.wind_speed(iray,:)));
    
    %find minimums and maximums only if at least 3 finite values in array
    if size(index,2)>3
        
        %-----
        
        % find maxima
        [P_max,L_max]= findpeaks(LLJ.wind_speed(iray,index),LLJ.range(index),'MinPeakProminence',1,'MinPeakWidth',50); %60
        
        height_max(iray, 1:length(L_max)) = L_max;
        ws_max(iray, 1:length(P_max)) = P_max;
        
        %find indeces of the maximums
        for i = 1:nvalues
            if ~isnan(height_max(iray,i))
                ind_max(iray, i) = find(LLJ.range == height_max(iray,i));
                dir_max(iray, i) = LLJ.wind_direction(iray,ind_max(iray,i));
            end
        end
        
        %-----
        
        %find bottom minimum
        ind = find(LLJ.range == height_max(iray,1));
        [ws_min_bottom ws_min_ind] = nanmin(LLJ.wind_speed(iray,1:ind));
        
        %find the height of a minimum
        height_min_bottom = LLJ.range(ws_min_ind);
        
        if  ~isempty(ws_min_bottom)
            ws_min(iray,1) = ws_min_bottom;
            height_min(iray,1) = height_min_bottom;
        end
        
        %-----
        
        %find minimumn between the maximums if there is more than one local maximums
        if sum(isnan(ws_max(iray,:))) < nvalues - 1
            for j = 1:nvalues-1
                if ~isnan(ind_max(iray,j+1))
                    %find nan between maximums
                    [ws_min_mid, ws_min_ind_mid] = nanmin(LLJ.wind_speed(iray,ind_max(iray,j):ind_max(iray,j+1)));
                    %ws_min_ind_mid:iin pit???? lis??t?? ekan maksimin indeksi ja v??hent???? yksi!
                    height_min_mid = LLJ.range(ind_max(iray,j)+ws_min_ind_mid-1);
                    
                    if ~isempty(ws_min_mid)
                        ws_min(iray,j+1) = ws_min_mid;
                        height_min(iray,j+1) = height_min_mid;
                    end
                end
            end
        end
        
        %-----
        
        %find top minimum
        %index where the top max are in ws_max array
        top_max_ind = max(find(~isnan(ws_max(iray,:))));
        
        if ~isempty(ws_max(iray,top_max_ind))
            real_top_max_ind = find(LLJ.wind_speed(iray,:) == ws_max(iray,top_max_ind)); %this is causing problems!!!
            if size(real_top_max_ind,2)>1
                real_top_max_ind = real_top_max_ind(end);
            end
            
            top_min_ind = max(find(~isnan(ws_min(iray,:))));
            
            %find the index where first NaN
            nan_boundary_max = top_max_ind+1;
            nan_boundary_min = top_min_ind+1;
            
            %find minimum above top_max_ind
            [ws_min_top ind_ws_min_top] = min(LLJ.wind_speed(iray,real_top_max_ind:end));
            ind_ws_min_top = ind_ws_min_top-1 + real_top_max_ind;
            
            if size(ind_ws_min_top,2)>1
                ind_ws_min_top = ind_ws_min_top(1);
            end
            
            %find the height of a top min
            height_min_top = LLJ.range(ind_ws_min_top);
            
            if ~isempty(ws_min_top) && (height_min_top > LLJ.range(real_top_max_ind))
                ws_min(iray,:) = [ws_min(iray,1:nan_boundary_min-1) ws_min_top ws_min(iray, nan_boundary_min+1:nvalues)];
                height_min(iray,:) = [height_min(iray,1:nan_boundary_min-1) height_min_top height_min(iray, nan_boundary_min+1:nvalues)];
            end
            
        end
        
        %-----
        
        %this end is for if size(index,2)>3
    end
    
    %find indeces for minimums and find the ws_min_directions
    for i = 1:nvalues
        if ~isnan(height_min(iray,i))
            ind_min(iray, i) = find(LLJ.range == height_min(iray,i));
            dir_min(iray, i) = LLJ.wind_direction(iray,ind_min(iray,i));
        end
    end
    
end

LLJ.height_max = height_max;
LLJ.height_min = height_min;
LLJ.ws_max = ws_max;
LLJ.ws_min = ws_min;
LLJ.dir_max = dir_max;
LLJ.dir_min = dir_min;
LLJ.ind_max = ind_max;
LLJ.ind_min = ind_min;

%% LLJ criteria

%number of allowed LLJs per profile num-1
%only num-1 max are tested
num=4;

lowleveljet(1:size(LLJ.time,1),1:num) = 0;
lowleveljet_speed(1:size(LLJ.time,1),1:num) = NaN;
lowleveljet_height(1:size(LLJ.time,1),1:num) = NaN;
lowleveljet_dir(1:size(LLJ.time,1),1:num) = NaN;
lowleveljet_speed_min(1:size(LLJ.time,1),1:num+1) = NaN;
lowleveljet_height_min(1:size(LLJ.time,1),1:num+1) = NaN;
lowleveljet_dir_min(1:size(LLJ.time,1),1:num+1) = NaN;

abs_limit = 2;
rel_limit = 1.25;

for iray = 1: size(LLJ.time,1)
    
    for j=1:num-1
        %varmistetaan etta maksimi ja ala ja yla-puolinen minimi ei nanneja
        if ~isnan(LLJ.ws_max(iray,j)) && ~isnan(LLJ.ws_min(iray,j)) && ~isnan(LLJ.ws_min(iray,j+1))
            if LLJ.ws_max(iray,j) > LLJ.ws_min(iray,j) + abs_limit && ...
                    LLJ.ws_max(iray,j) > LLJ.ws_min(iray,j+1) + abs_limit && ...
                    LLJ.ws_max(iray,j) > LLJ.ws_min(iray,j) * rel_limit && ...
                    LLJ.ws_max(iray,j) > LLJ.ws_min(iray,j+1) * rel_limit
                
                lowleveljet(iray,j) = 1;
                %LLJ max
                lowleveljet_speed(iray,j)  = LLJ.ws_max(iray,j);
                lowleveljet_height(iray,j) = LLJ.height_max(iray,j);
                lowleveljet_dir(iray,j)    = LLJ.dir_max(iray,j);
                %LLJ minima
                lowleveljet_speed_min(iray,j:j+1)  = [LLJ.ws_min(iray,j) LLJ.ws_min(iray,j+1)];
                lowleveljet_height_min(iray,j:j+1) = [LLJ.height_min(iray,j) LLJ.height_min(iray,j+1)];
                lowleveljet_dir_min(iray,j:j+1)    = [LLJ.dir_min(iray,j) LLJ.dir_min(iray,j+1)];
            end
        end
    end
    
end

LLJ.lowleveljet = lowleveljet;
LLJ.speed = lowleveljet_speed;
LLJ.height = lowleveljet_height;
LLJ.dir = lowleveljet_dir;
LLJ.speed_min = lowleveljet_speed_min;
LLJ.height_min = lowleveljet_height_min;
LLJ.dir_min = lowleveljet_dir_min;
%LLJ.nval = LLJ.nval;

%% Additional criteria
% elapsed time (use matlab time)

LLJ.time_n  = daten + LLJ.time./24;

% Find jet duration
% need to keep track of jet - give it a label!
% Jet A, B, C
% Keep speed, direction, height

% thresholds for determining whether jet is a new one,
% or a continuation of a previously found jet

% time is in matlab time
thresh_time = 1/24;     % 1 hour --> %kesto v??h. puoli tuntia jos thres_time/2
thresh_height = 135;    % 135 m, four gates in DBS LLJ
thresh_speed = 0.3;     % percentage value 20 % --> if earlier jet more than 20% smaller than next jet --> not good..
thresh_dir = 45;        % 45 degrees, this is a conservative number to remove isolated maxima
thresh_duration = 1/24; %duration threshold ==> same as time thresh

% Keep track of active jets, as we loop over these..

active_jets = [];
next_jet = 0;

% output LLJ
LLJout.time_n         = [];
LLJout.height     = [];
LLJout.speed      = [];
LLJout.dir        = [];
LLJout.height_min = [];
LLJout.speed_min  = [];
LLJout.dir_min    = [];
LLJout.label      = [];

% storage for individual active jets
aj.time_n         = [];
aj.height     = [];
aj.speed      = [];
aj.dir        = [];
aj.height_min = [];
aj.speed_min  = [];
aj.dir_min    = [];
aj.label      = [];

for itime = 1:size(LLJ.time_n,1)
    for ijet = 1:size(LLJ.lowleveljet,2) - 1
        if LLJ.lowleveljet(itime,ijet) == 1
            % found a jet
            
            % are there any active jets? if not assign to new jet
            if isempty(active_jets)
                next_jet = next_jet + 1;
                active_jets = next_jet;
                
                % new active jet so no need to test thresholds
                aj(next_jet).time_n     = LLJ.time_n(itime);
                aj(next_jet).height     = LLJ.height(itime,ijet);
                aj(next_jet).speed      = LLJ.speed(itime,ijet);
                aj(next_jet).dir        = LLJ.dir(itime,ijet);
                aj(next_jet).height_min = LLJ.height_min(itime,ijet:ijet+1);
                aj(next_jet).speed_min  = LLJ.speed_min(itime,ijet:ijet+1);
                aj(next_jet).dir_min    = LLJ.dir_min(itime,ijet:ijet+1);
                aj(next_jet).label      = next_jet;
            else
                
                is_new_jet = 1;
                
                % test jet against active jets
                for i = 1:length(active_jets)
                    
                    if ( abs(aj(active_jets(i)).time_n(end) - LLJ.time_n(itime)) < thresh_time ...
                            & abs(aj(active_jets(i)).height(end) - LLJ.height(itime,ijet)) < thresh_height  ...
                            & abs((LLJ.speed(itime,ijet) - aj(active_jets(i)).speed(end)) / LLJ.speed(itime,ijet)) < thresh_speed  ...
                            & min(mod(aj(active_jets(i)).dir(end) - LLJ.dir(itime,ijet), 360), ...
                            360 - mod(aj(active_jets(i)).dir(end) - LLJ.dir(itime,ijet), 360)) < thresh_dir )
                        
                        % is a continuation of an active jet
                        aj(active_jets(i)).time_n         = [aj(active_jets(i)).time_n;         LLJ.time_n(itime)];
                        aj(active_jets(i)).height     = [aj(active_jets(i)).height;     LLJ.height(itime,ijet)];
                        aj(active_jets(i)).speed      = [aj(active_jets(i)).speed;      LLJ.speed(itime,ijet)];
                        aj(active_jets(i)).dir        = [aj(active_jets(i)).dir;        LLJ.dir(itime,ijet)];
                        aj(active_jets(i)).height_min = [aj(active_jets(i)).height_min; LLJ.height_min(itime,ijet:ijet+1)];
                        aj(active_jets(i)).speed_min  = [aj(active_jets(i)).speed_min;  LLJ.speed_min(itime,ijet:ijet+1)];
                        aj(active_jets(i)).dir_min    = [aj(active_jets(i)).dir_min;    LLJ.dir_min(itime,ijet:ijet+1)];
                        aj(active_jets(i)).label      = [aj(active_jets(i)).label;      active_jets(i)];
                        
                        is_new_jet = 0;
                        break;
                    end
                end
                
                if is_new_jet
                    next_jet = next_jet + 1;
                    active_jets = [active_jets next_jet];
                    % new active jet so no need to test thresholds
                    aj(next_jet).time_n         = LLJ.time_n(itime);
                    aj(next_jet).height     = LLJ.height(itime,ijet);
                    aj(next_jet).speed      = LLJ.speed(itime,ijet);
                    aj(next_jet).dir        = LLJ.dir(itime,ijet);
                    aj(next_jet).height_min = LLJ.height_min(itime,ijet:ijet+1);
                    aj(next_jet).speed_min  = LLJ.speed_min(itime,ijet:ijet+1);
                    aj(next_jet).dir_min    = LLJ.dir_min(itime,ijet:ijet+1);
                    aj(next_jet).label      = next_jet;
                end
            end
        end
    end
    
    %remove old jets from active_jets
    old_jet_index = [];
    for i = 1:length(active_jets)
        % we can't use the thresh_time here because of the sensitivity analysis
        % so instead let's put (1/24)
        %if ( abs(aj(active_jets(i)).time_n(end) - LLJ.time_n(itime)) >= thresh_time*2.5)
        if ( abs(aj(active_jets(i)).time_n(end) - LLJ.time_n(itime)) >= (1/24)*2.5)
            old_jet_index = [old_jet_index i];
        end
    end
    if ~isempty(old_jet_index)
        active_jets(old_jet_index) = [];
    end
    
end

% test jet duration
cleanindex = [];
for j = 1:length(aj)
    if length(aj(j).time_n) < 3
        % jet not long enough
        cleanindex =  [cleanindex j];
    elseif (aj(j).time_n(end) - aj(j).time_n(1) <= thresh_duration) %kesto v??h. 60min
        % jet not long enough
        cleanindex =  [cleanindex j];
    end
end
% remove short-lived jets
aj(cleanindex) = [];

% relabel jets
for j = 1:length(aj)
    aj(j).label(:) = j;
end


% write out jets in start time order to LLJ struct
for i = 1:length(aj)
    LLJout.time_n     = [LLJout.time_n;         aj(i).time_n];
    LLJout.height     = [LLJout.height;     aj(i).height];
    LLJout.speed      = [LLJout.speed;      aj(i).speed];
    LLJout.dir        = [LLJout.dir;        aj(i).dir];
    LLJout.height_min = [LLJout.height_min; aj(i).height_min];
    LLJout.speed_min  = [LLJout.speed_min;  aj(i).speed_min];
    LLJout.dir_min    = [LLJout.dir_min;    aj(i).dir_min];
    LLJout.label      = [LLJout.label;      aj(i).label];
end

end

