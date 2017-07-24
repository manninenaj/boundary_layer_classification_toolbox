function [data,att] = loadModel(site,daten,m_type)
%loadModel
switch site
    case 'hyytiala'
        pathmodel = ['/Users/anjmanni/MEAS_DATA/GDAS/' site '/' ...
            datestr(daten,'yyyymmdd') '_' site '_' m_type '.nc'];
        [data,att] = load_nc_struct_silent(pathmodel);
        
    case 'juelich'
        pathmodel = ['/data/hatpro/jue/cloudnet/' site '/calibrated/' m_type '/' datestr(daten,'yyyy') '/' ...
            datestr(daten,'yyyymmdd') '_' site '_' m_type '.nc'];
        [data,att] = load_nc_struct(pathmodel);
        
        for i=1:size(data.height,2)
           data.range(i,1) = nanmean(data.height(:,i)); 
        end
        
    case 'arm-oliktok'
        pathmodel = ['/data/hatpro/jue/cloudnet/juelich/calibrated/' m_type '/oli/' ...
            'oliecmwfvarX1.c1.' datestr(daten,'yyyymm') '01.000000.cdf'];
        [model,att] = load_nc_struct(pathmodel);
        
        ind = strcmp(datestr(daten,'yyyymmdd'),cellstr(num2str(model.idat)));
        x = find(ind == 1);
        
        if isempty(x)
            pathmodel = ['/data/hatpro/jue/cloudnet/juelich/calibrated/' m_type '/oli/' ...
                'oliecmwfvarX1.c1.' datestr(daten+1,'yyyymm') '01.000000.cdf'];
            [model,att] = load_nc_struct(pathmodel);
            ind = strcmp(datestr(daten,'yyyymmdd'),cellstr(num2str(model.idat)));
            x = find(ind == 1);
        end

        h = zeros(size(model.p,2),size(x,1));
        h(1,:) = 10; 
        for j = 1:size(x,1) 
            for i = 1:size(model.p,2)-1
                h(i+1,j) = h(i,j)+(model.p(x(j),end-i)-model.p(x(j),end-i+1))/(-1.2985*9.80665);
            end
        end
        
        for i = 1:size(h,1)
            data.range(i,1) = mean(h(i,:));
        end
        
        data.uwind = model.u(x,end:-1:1);
        data.vwind = model.v(x,end:-1:1);
        data.time = linspace(0,23,24);
        
    otherwise
        data = []; att = [];
        warning('Site %s not specified yet!',site)
end
end
