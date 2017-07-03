function [data,att] = loadModel(site,daten,m_type)
%loadModel
switch site
    case 'hyytiala'
        pathmodel = ['/Users/anjmanni/MEAS_DATA/GDAS/' site '/' ...
            datestr(daten,'yyyymmdd') '_' site '_' m_type '.nc'];
        [data,att] = load_nc_struct_silent(pathmodel);
    otherwise
        data = []; att = [];
        warning('Site %s not specified yet!',site)
end
end
