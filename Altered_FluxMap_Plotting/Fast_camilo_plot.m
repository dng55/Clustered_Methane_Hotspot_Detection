% Fast Camilo Plotting
clear all;

% this_site = 'US-WPT'; % Change this to re-direct all scripts below to the relevant site settings
this_site = 'Young';

% FP_Opts_BB
cd '/Users/darianng/Documents/MSc_Geography/Methane_Hotspot/FARF_Code/Footprints-master/CalculateFootprint/Run_options'

tstart = '2022-07-017 00:00';
tend = '2022-08-017 23:30';


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

Output_as_csv;