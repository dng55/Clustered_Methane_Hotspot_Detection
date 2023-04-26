%% Main Footprint_Run
% By Camilo Rey. Sep 2020.

% Before running you should have done the following:

% 1-	Go to ‘/CalculateFootprint/Create_Site_Input’ and modify the options for your site in the function ‘load_site_input.m’. Here you should find the desired site or create a new one, and define the data directory.
% 2-	Go to ‘/CalculateFootprint/Run_options’, create a copy of the template file “FP_opts_xx.m’, and create your customized options file. The name of this FP_opts file should reflect the characteristics of this run, for example, “FP_opts_EE_2019.m”.
% 3-	Run Footprint_Run_Continuous.m



%% Initial Adjustments

if opt.DN==0;DT='night';elseif opt.DN==1; DT='day';end  
Fall='';
for i=1:length(opt.FPcontour)
Flab=num2str(opt.FPcontour(i));
Fall=strcat(Fall,'-',Flab);
end

lega=Fall;% String of the customized contours to label files    
%% Load data


addpath('Proc_functions/')
addpath('Create_Site_Input/')

if opt.DataType==1 % 1: BiometLab Mat files (ready), 2: Ameriflux CSV (pending)
    
[data,window,lat_site,lon_site,Tower_HT,Canopy_HT,map_path]=load_site_data(opt.site);

data.WD=data.WD-opt.WDoffset; data.WD(data.WD<0)=data.WD(data.WD<0)+360;

% Load PBL data 
data=Load_PBL_data(opt.site,data,opt.PBLdata);


%% Index desired time

% Daytime vector
% dayvec=data.ze<90;% daytime= 1. Nighttime=0
dayvec = data.daytime;

start_date=datenum(opt.start,'yyyy-mm-dd HH:MM');
end_date=datenum(opt.end2,'yyyy-mm-dd HH:MM');

% wind direction index (adjust if necessary)
wix=data.WD>opt.WDlow | data.WD<opt.WDhigh;
% Soil temp index if desired
temp=data.ST>opt.STlow & data.ST<opt.SThigh;
  

% Index all time within this period
if opt.BOTH==1
    ix=data.Mdate>=start_date & data.Mdate<=end_date & wix & temp;
elseif opt.BOTH==0
    ix=data.Mdate>=start_date & data.Mdate<=end_date & dayvec==opt.DN & wix & temp;
    ix2=data.Mdate>=start_date & data.Mdate<=end_date & dayvec==1 & wix & temp;%daytime index
end

if exist('ix2','var') == 0; ix2=ix;end


   
%% Apply index to input variables:

    Mdate=data.Mdate(ix);
    plottime=datetime(datevec(Mdate));
    daynight=dayvec(ix);
    startH=data.time(ix)-window;
    endH=data.time(ix);
    h=data.PBL(ix);
    ustar=data.ustar(ix);
    ubar=data.ubar(ix);
    mbar=data.mbar(ix);
    vv=data.vv(ix);
    Lo=data.L(ix);
    windir=data.WD(ix);
    DOY=data.DOY(ix);  
    year=data.year(ix);
    TA=data.TA(ix);
    ST=data.ST(ix);

    %Optional
    LE=data.LE(ix);
    FCH4=data.wm(ix);%% 
    FCO2=data.wc(ix);%% 
    FH2O=data.wq(ix);%% 
    Fuv=(data.ustar(ix)).^2;%% Momentum flux
    SH=data.H(ix);%% 
    
    Canopy_ht=Canopy_HT(ix);
    Tower_height=data.z(ix);
    
    Mdate_day=data.Mdate(ix2);
    ustar_day=data.ustar(ix2);
    Lo_day=data.L(ix2);
    ubar_day=data.ubar(ix2);
    Tower_height_day=data.z(ix2);
else
    
% Do Ameriflux csv here

end

%% load map or Define Grid

zm=nanmean(Tower_height)-nanmean(Canopy_ht); % Average tower height
F2H=round(125*zm,-2);

if opt.HaveMap==1
    load(map_path)
    PatchMap=MapB;
    
elseif opt.mscale==0
    
    lgt=F2H*2 ;% Lenght of one side of the square with the tower in the center
    PxSize=lgt/100;% Recommended 100*100 pixels
    FX=-lgt/2:PxSize:lgt/2;
    FY=-lgt/2:PxSize:lgt/2;
    PatchMap=nan;
    PatchCode='NA';
    Cmap=[1,1,1];
    
elseif opt.mscale==1
    
    lgt=F2H*2 ;% Lenght of one side of the square with the tower in the center    
    FX=[flip(-exp(0:0.1:log(lgt/2))),0, exp(0:0.1:log(lgt/2))];
    FY=[flip(-exp(0:0.1:log(lgt/2))),0, exp(0:0.1:log(lgt/2))];
    PxSize=round(nanmean(FX(2:end)-FX(1:end-1)),0);% Recommended no more than 5 m (m)
    PatchMap=nan;
    PatchCode='NA';
    Cmap=[1,1,1];
end

 %% Roughness

z=Tower_height;%   
z2=Tower_height_day;%   
if opt.calcZo==1
    disp('calculating roughness length')
    if length(~isnan(Lo_day))>200 
        
        % Estimate h based on Pennypacker and Baldocchi (2015) paper 
        if opt.calcH==1 
           daysAVG=2;%days to average
           [h_dsk,ZvegSmooth]=CalculateAerodynamicCanopyHeight(Mdate_day,ustar_day,z2,Lo_day,ubar_day,daysAVG,1);
           Canopy_height=nanmean(ZvegSmooth); 
           if opt.DN==1
                d = 0.66*ZvegSmooth; % m - zero plane dispacement height, should vary witht time
            elseif opt.DN==0
                d = ones(length(ubar),1)*nanmean(ZvegSmooth); % m - zero plane dispacement height, should vary witht time
            end
        else
            ZvegSmooth=ones(length(ubar),1)*Canopy_height;
            d = 0.66*ZvegSmooth;
        end
        % Roughness lenght should be constant. 
        % Calculate Roughness Length using Mauer and Bohrer (2016)

        [d3,zo3,n3,stat3]=RoughnessLength(nanmean(Tower_height),Canopy_height,Lo_day,ustar_day,ubar_day);
            zo=ones(length(ubar),1)*zo3;
            
        % Tonzi Case
            if strcmp (opt.site,'Tonzi_under') | strcmp (opt.site,'Vaira')
                disp('Linear increase in canopy height in Tonzi_under or Vaira')
                ZvegSmooth=0.5*ones(length(ubar),1);
                ZvegSmooth(1:24*30)=(CHinit-CHend)./(24*30).*(1:24*30)+CHinit;% y= mx=b
                d = 0.66*ZvegSmooth; % m - zero plane dispacement height, should vary witht time
            end

     else
            ZvegSmooth=Canopy_ht;
            d = 0.66*ZvegSmooth;
            zo=0.1*ZvegSmooth;
            disp('Not enough data to calculate aerodynamic canopy height. Prescribing constant canopy height')
     end
end    

%% Plotting adjustments
if opt.plotYN==1
PlotVec=ones(length(dayvec(ix)),1);% Ones for plot, zeros for no plot
else
PlotVec=zeros(length(dayvec(ix)),1);% Ones for plot, zeros for no plot
end            
    % Set a point of interest
    Tdist=ones(size(windir))*opt.distP; % 
    Twindir=ones(size(windir))*opt.dirP; % 
    

%% Run Footprint Code

    [COUNT,COMP,PATCHSUM_H,PATCHSUM_K,PATCHSUM_KM]=...
    FP_process_SurfaceFlux(ustar,vv,Lo,windir,zo,z,d,ubar,FX,FY,PatchMap,Cmap,PlotVec,opt.site,DOY,year,h,startH,endH,...
    opt.FPcontour,opt.models,opt.mscale,Tdist,Twindir,FCH4,FCO2,FH2O,Fuv,SH,ST);

%% Variable Description:

metadata.variables = {'COUNT','Number of footprints in current subset','';...
             'PATCHSUM','A matrix with the percent contributions of each patch to the footprint (see PatchCode for correspondance)','';...
             'Mdate','Matlab date & time @ end of flux averaging period','';...
             'PatchCode','Code of each land cover within the opt.site (if available)','';...
             'footCount',' A matrix that counts the times that each cell was included in the footprint','';...
             'perF','Percent Footprint calculated from the cumulative FP','';...
             'MapAvg','Equals percF*100','';...
             'MapB','The original Map of patch distribution','';...
             'Cmap','Colormap used in the Map of the opt.site','';...
             'FCH4','CH4 flux','nmonl m-2 s-1';...
             'PxSize','The size of the pixel in the map used','m'};
         
%% Plot Climatology for current subset and save
addpath(genpath('googleearth'))    
% Compute coordinates to plot    

    global R_e Earth_E2 omega_e

    % radius of the Earth
    R_e = 6.37813649e6; % [m]

    % Earth's shape - eccentricity^2
    Earth_E2 = 0.006694385000^2;      % [-]

    % Mean Angular Velocity of the Earth
    omega_e =  7.29211585530e-5;    % [rad/s]
    
    lat_gc = lat_site*pi/180;
    lon_gc = lon_site*pi/180;
    [LATsave,LONsave]=getCartesian(lat_gc,lon_gc,FX,FY);
    [FX2, FY2] = meshgrid(FX,FY);%  
    RT= length(FY);
    
    
modLab={'Hsieh','Kljun','K&M'};
    
for i=1:3
    FileLabel=[opt.site '-' modLab{i} '-' opt.start(1:4) opt.start(9:11) '-' opt.start(13:14) opt.start(16:17) '-to-' opt.end2(1:4) opt.end2(9:11) '-' opt.end2(13:14) opt.end2(16:17) '-' DT lega];
    if opt.models(i)==1
     
    % Plot Base Map
    figure;
    if ~isnan(PatchMap)
        pcolor(FX2,FY2(RT:-1:1,:),PatchMap*100);shading('flat');% This is the same as using flidud, which flips the rows so that the first one is the last one
        colormap(Cmap);
    else
        pcolor(FX2,FY2(RT:-1:1,:),nan(length(FX),length(FY)));shading('flat');% This is the same as using flidud, which flips the rows so that the first one is the last one
    end    
        hold on; xlabel('x [m]');  ylabel('y [m]')
        
    % Cumulative Footprint
    FG(i).percF=CalcPercF_fast(COMP(i).footCUM/COUNT)*100;
    [~,cA]=contour(FX2,FY2(RT:-1:1,:),FG(i).percF,opt.FPcontour,'Fill','off','ShowText','off');%shading('interp') 
    cA.LineColor = 'k'; cA.LineWidth = 1;hold on;   
    title(['Climatology ' FileLabel '%'])
    end
end
    legend(modLab)
for i=1:3
    FileLabel=[opt.site '-' modLab{i} '-' opt.start(1:4) opt.start(9:11) '-' opt.start(13:14) opt.start(16:17) '-to-' opt.end2(1:4) opt.end2(9:11) '-' opt.end2(13:14) opt.end2(16:17) '-' DT lega];
    % Export to kml
    if opt.models(i)==1
%     figure
%     kml_contour(LONsave,LATsave(RT:-1:1,:),FG(i).percF,...
%     [opt.SaveDir 'kml\' FileLabel 'percent.kml'],opt.FPcontour,'r')
    end
end

%% Save all output    

clear data
%     FileLabel = 'Hogg_output';
    FileLabel = [this_site '_output'];
%     save([opt.SaveDir FileLabel opt.Sufix '.mat'])
    save([opt.SaveDir FileLabel '.mat'])




