function [ Fubarfield ] = associateTKEwith(Eps,skewn,cloudmask,th,cut_h)
%TKE_incontact_w creates a mask indicating is the turbulence
%in contact with surface, cloud, both, or neither of them.

% rm nans
skewn(isnan(Eps)) = nan;

% Skewness field
%0: no signal
%1: pos skewness connected with surface
%2: neg skewness
%3: in cloud
%4: neg skewness connected with cloud
%5: pos skewness unconnected with surface
%6: no signal

% Initialize
Skewnfield = zeros(size(Eps));
Skewnfield(skewn>0) = 1;
Skewnfield(skewn<0) = 2;
Skewnfield(cloudmask) = 3; % cloud

% Find negative skewness connected with cloud
for ii = 1:size(Skewnfield,1)
    if any(cloudmask(ii,:),2) % cloud in profile?
        for jj = size(Skewnfield,2):-1:cut_h+1
            % cloud top
            if Skewnfield(ii,jj)==3 % if cloud
                jj2 = jj;
                % cloud bottom
                while jj2>cut_h+1 && Skewnfield(ii,jj2)==3
                    jj2 = jj2-1;
                end
                % if neg. skewness
                if jj2>cut_h+1 && (Skewnfield(ii,jj2)==2 || Skewnfield(ii,jj2)==0) 
                    jj3=jj2;
                    % while neg. skewness
                    while jj3>cut_h+1 && (Skewnfield(ii,jj3)==2 || Skewnfield(ii,jj3)==0) 
                        jj3 = jj3-1;
                        % if pos. but next range gate is negative, i.e.
                        % allow one pos. in between neg. gates
                        if Skewnfield(ii,jj3)~=2 && Skewnfield(ii,jj3-1)==2
                            jj3 = jj3-1;
                        end
                    end
                    % From cloud bottom to when skewness changes to pos.
                    Skewnfield(ii,jj2:-1:jj3+1)=4;
                end
                
            end
        end
    end
end

%0: no signal
%1: non turbulent
%2: connected with surface
%3: connectedt with cloud
%4: in cloud
%5: unconnected

Fubarfield = zeros(size(Eps));
Fubarfield(~isnan(Eps)) = 1;
is_hi_eps = real(log10(Eps)) > th;
Fubarfield(is_hi_eps) = 2; % hi turbulence
Fubarfield(is_hi_eps & Skewnfield==4) = 3; % connected w/ cloud
Fubarfield(cloudmask) = 4; % in cloud
Fubarfield(Fubarfield == 3 & repmat(~any(cloudmask,2),1,size(Fubarfield,2))) = 2; % if cloud driven but no clouds in profile
Fubarfield(:,1:cut_h) = nan; % ignore
Epsfield = zeros(size(Eps));
Epsfield(real(log10(Eps))>th)=1;

% Find epsilon not in contact with surface
for ii = 1:size(Epsfield,1)
    izero = find(Fubarfield(ii,cut_h+1:end)~=2,1,'first');
    ione = find(Fubarfield(ii,cut_h+1:end)==2);
    if ~isempty(ione)
        ione(ione<izero) = []; % if '1' higher than '0', remove
        Fubarfield(ii,cut_h+ione) = 5;
    end
end
Fubarfield(isnan(Eps)) = 0;
Fubarfield(isnan(Fubarfield)) = 0; % no signal
end

