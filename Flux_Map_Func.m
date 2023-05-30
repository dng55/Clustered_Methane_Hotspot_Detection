function [endFileName] = Flux_Map_Func(this_site,tstart,tend)
% Turning Fast_camilo_plot.m into a function.
%   Inputs:
%         this_site = Site name ('Young', 'Hogg', 'US-Myb', or 'US-WPT')
%         tstart = starting time ('2020-10-025 00:00')
%         tend = ending time ('2020-11-025 00:00')
%   Outputs:
%         Output savefile name.

% Fast Camilo Plotting

% FP_Opts_BB
cd '/Users/darianng/Documents/MSc_Geography/Methane_Hotspot/FARF_Code/Footprints-master/CalculateFootprint/Run_options'



if strcmp(this_site,'Hogg')
    FP_opts_Hogg;
elseif strcmp(this_site,'Young')
    FP_opts_Young;
elseif strcmp(this_site,'BB1')
    FP_opts_BB;
elseif strcmp(this_site,'BB2')
    FP_opts_BB2;
elseif strcmp(this_site,'US-Myb')
    FP_opts_Myb;
elseif strcmp(this_site,'US-WPT')
    FP_opts_WPT;
end


% Footprint_Run_Continuous
cd '/Users/darianng/Documents/MSc_Geography/Methane_Hotspot/FARF_Code/Footprints-master/CalculateFootprint'
Footprint_Run_Continuous;
% Surface_Flux_Map
cd '/Users/darianng/Documents/MSc_Geography/Methane_Hotspot/FARF_Code/Footprints-master/Analyze_Outputs'
Surface_Flux_Map;
% Saving data as csv
cd '/Users/darianng/Documents/MSc_Geography/Methane_Hotspot/FARF_Code'

% Creating desired file name - Choosing format: 202106041-202106030.csv
endFileName = [erase(extractBefore(tstart,12),'-') '-' erase(extractBefore(tend,12),'-') '.csv'];
disp(['Saved files ending with: ' endFileName])

% Output_as_csv;
col1 = FX2;
col2 = FY2(RT:-1:1,:);

savePath = ['python_code/data/compilation/',this_site];
cd(savePath)

% endFileName = 'May_Aug2021.csv';

csvwrite([this_site,'_fluxMap_x_',endFileName],col1);

csvwrite([this_site,'_fluxMap_y_',endFileName],col2);

csvwrite([this_site,'_fluxMap_ch4_',endFileName],Lflux_ch4);

csvwrite([this_site,'_fluxMap_h_',endFileName],Lflux_h);

csvwrite([this_site,'_fluxMap_co2_',endFileName],Lflux_co2);

end