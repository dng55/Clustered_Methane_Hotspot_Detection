col1 = FX2;
col2 = FY2(RT:-1:1,:);
% col3 = Lflux;

% combined = [col1;col2;col3]';

% csvwrite('FFP_csv_outputs/BB_fluxMap_x.csv',col1);
cd python_code/data

% endFileName = 'May_Aug2021.csv';

csvwrite([this_site,'_fluxMap_x_',endFileName],col1);

csvwrite([this_site,'_fluxMap_y_',endFileName],col2);

csvwrite([this_site,'_fluxMap_ch4_',endFileName],Lflux_ch4);

csvwrite([this_site,'_fluxMap_h_',endFileName],Lflux_h);

csvwrite([this_site,'_fluxMap_co2_',endFileName],Lflux_co2);

% csvwrite(['Hogg_fluxMap_x_',endFileName],col1);
% 
% csvwrite(['Hogg_fluxMap_y_',endFileName],col2);
% 
% csvwrite(['Hogg_fluxMap_ch4_',endFileName],Lflux_ch4);
% 
% csvwrite(['Hogg_fluxMap_h_',endFileName],Lflux_h);
% 
% csvwrite(['Hogg_fluxMap_co2_',endFileName],Lflux_co2);