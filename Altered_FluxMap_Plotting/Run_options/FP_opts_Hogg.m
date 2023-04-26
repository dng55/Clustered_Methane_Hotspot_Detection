%% Input footprint options

% By Camilo Rey. Jul 2020.

%%%%%% Start loading Footprint options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clearvars -except data % COMMENTED FOR FAST PLOTTING
       
    % What type of data?
    
    opt.DataType=1;% 1: BiometLab Mat files (ready), 2: Ameriflux CSV (pending)
    opt.site='Hogg';
    
    % how?
    opt.models=[1 0 1]; % Select any combination of 3 models [ Hsieh  Kljun  K&M ];
    opt.mscale=0; % The scale of the initial unrotated matrix. Must be zero for now
    opt.mscale2=0;% The scale of the final, rotated matrix. Must be zero for now

    opt.PBLdata=0; % 0 if there is no PBL data available
    opt.calcZo=1; %Calculate roughness length
    opt.calcH=1; % Calculate aerodynamic canopy heigth

    opt.HaveMap=0;% 1= Yes, I have a raster with land covers, 0= I dont have one
    opt.plotYN=0;% Ones for plot, zeros for no plot. % Normally stays at zero if running more than ~30 footprints

    % when?
    opt.BOTH=0; % Both daytime and nighttime at the same time, 1=yes, 0= no,separetely
    opt.DN=1;% daytime= 1. Nighttime=0; (only works if opt.BOTH==0)
    
    % Defining an index to run the footprint in day of year        
%     opt.start='2019-06-152 00:00';
%     opt.end2='2019-12-365 24:00';

    opt.start=tstart;
    opt.end2=tend;

    opt.FPcontour=[50,80];% Choose Footprint Percentages. no more than 60-70 % for nighttime
    
    % Optional: Wind direction adjustment (if needed)
    opt.WDoffset=0;
    
    % Optional: Select a point around the tower to evaluate
    opt.distP=200; % Distance from the tower
    opt.dirP=360; % Azimuth
    
    % Optional: Filter based on wind direction or temperature
    opt.WDlow=0;
    opt.WDhigh=360;
    opt.STlow=15;
    opt.SThigh=20;
    
    % Saving options
    opt.SaveDir=(['../Footprint_Output/' opt.site '/']);
    opt.Sufix=['_T_' num2str(opt.STlow) 'to' num2str(opt.SThigh) '_C2'];% A short description of this run for tracking purposes
    
disp('Options recorded. Run Footprint_Run_xx.m')  

%%%%%%%%% Stop entering footprint options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 