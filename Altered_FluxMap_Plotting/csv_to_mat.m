% csv_to_mat loads fluxtower csv and creates:
%   1) An Mdate column 
%   2) 0cm soil temp from stefan boltzmann
%   3) Zenith angle
% and saves as .mat file.
clear all;

Site = '2021-Young';
% Site = 'US-WPT';
% Site = 'US-Myb';
% Site = '2021-Hogg';
% Site = '2014 - Burns Bog';
% Site = '2019-Burns Bog 2';
if strcmp(Site,'2019-Burns Bog 2')
    SiteName = 'BB2';
elseif strcmp(Site,'2014-Burns Bog')
    SiteName = 'BB';
elseif strcmp(Site,'2021-Young')
    SiteName = 'Young';
elseif strcmp(Site,'2021-Hogg')
    SiteName = 'Hogg';
else SiteName = Site;
end

% Checks if path name exists in UBC micromet drive. 
if exist(['/Volumes/GoogleDrive/My Drive/Micromet Lab/Projects/' Site '/Flux-tower/Flux_data'], 'dir')
    flux_path = ['/Volumes/GoogleDrive/My Drive/Micromet Lab/Projects/' Site '/Flux-tower/Flux_data'];
    cd(flux_path)
    data = readtable([SiteName '_L3.csv'],'TreatAsEmpty',{'NA'});
else
    disp('- Site does not exist in UBC micromet drive - Load data manually');
    disp('   - Loading manually from: /Users/darianng/Documents/MSc_Geography/Methane_Hotspot/Site_Data');
    cd('/Users/darianng/Documents/MSc_Geography/Methane_Hotspot/Site_Data');
    data = readtable([SiteName '.csv'],'TreatAsEmpty',{'NA'});
end 
% For BB only
% met_path = ['/Volumes/GoogleDrive/My Drive/Micromet Lab/Projects/' Site '/Flux-tower/met_data/met_merged/'];
% cd(met_path)
% met = readtable(['met_corrected_gapfilled' SiteName '.csv'],'TreatAsEmpty',{'NA'});

% %TEMPORARY FIX FOR MISMATCHED BIOMET FILE IN BB2 and "TIME" ISSUE
% data = removevars(data,{'time'});
% number = (32694-32300)+1;
% met = met(number:end,:);
% %TEMPORARY FIX FOR MISMATCHED BIOMET FILE IN BB2

cd '/Users/darianng/Documents/MSc_Geography/Methane_Hotspot/FARF_Code/Footprints-master/CalculateFootprint'
% Converting datetime to matlab time.
N = length(data.DATE);

Mdate = ones(N,1);
TS_0cm = ones(N,1);
ze = ones(N,1);

% Stefan-Boltzmann constant/parameter
sigma = 5.670373e-8;
if any(strcmp('LWOUT_1_1_1',data.Properties.VariableNames))
    E = data.LWOUT_1_1_1;
elseif any(strcmp('LONGWAVE_OUT',data.Properties.VariableNames))
    E = data.LONGWAVEOUT;
elseif any(strcmp('LW_OUT',data.Properties.VariableNames))
    E = data.LW_OUT;
end
e = 0.97;
sigma_adj = sigma * e;

% Zenith parameters
if strcmp(SiteName, 'Hogg')
    location.latitude = 50.371;
    location.longitude = -100.534;
    location.altitude = 0;
elseif strcmp(SiteName,'Young')
    location.latitude = 50.362;
    location.longitude = -100.202;
    location.altitude = 0;
elseif strcmp(SiteName,'BB2')
    location.latitude = 49.11897;
    location.longitude = -122.99511;
    location.altitude = 0;
elseif strcmp(SiteName,'US-Myb')
    location.latitude = 38.0499;
    location.longitude = -121.765;
    location.altitude = 0;
elseif strcmp(SiteName,'US-WPT')
    location.latitude = 41.4646;
    location.longitude = -82.9962;
    location.altitude = 0;
end

% Time offset to get GMT 
if strcmp(SiteName,'BB') || strcmp(SiteName,'BB2')
    gtG = -7/24; % Offset time for PDT
elseif strcmp(SiteName,'Young') || strcmp(SiteName,'Hogg')
    gtG = -5/25; % Offset time for CDT (Manitoba)
elseif strcmp(SiteName,'US-WPT')
    gtG = -5/25; % Offset time for EST (Ohio)
else
    disp('Unspecified timezone: Default to PDT')
    gtG = -7/24;
end

for i = 1:length(Mdate)
%     stringDate = num2str(data.TIMESTAMP(i));
%     years = str2num(stringDate(1:4));
%     month = str2num(stringDate(5:6));
%     day = str2num(stringDate(7:8));
%     hour = str2num(stringDate(9:10)); 
%     minute = str2num(stringDate(11:12));
%     seconds = 00;
%     thisDate = [years,month,day,hour,minute,seconds];
%     Mdate(i) = datenum(thisDate);
%     year(i) = years;

    % Getting datenum, checking if DATE is a string or num.
    if ischar(data.DATE(i))
        Mdate(i) = datenum(data.DATE(i));
    elseif isdatetime(data.DATE(i))
        Mdate(i) = datenum(data.DATE(i));
    else
        Mdate(i) = datenum(num2str(data.DATE(i)),'yyyymmddHHMM');
    end
%     E_num = str2double(E(i));
    if ~isnan(E(i))
        TS_0cm(i) = (E(i)/sigma_adj)^.25 - 273.15;
    else
        TS_0cm(i) = NaN;
    end
    
    % displaying progress bar
    if rem(i,1000) == 0
        disp(['Progress: ' num2str(i) '/' num2str(N)]);
    end
    
    % Calculating Zenith angle
    zx = datestr(Mdate(i)+gtG,'dd-mmm-yyyy HH:MM:SS');
    sun = sun_position(zx, location);
    ze(i)=sun.zenith;
    
end



data.Mdate = Mdate;
% data.year = year;
data.ze = ze;
met.Mdate_met = Mdate;
met.TS_0cm = TS_0cm;


% clear year, clear month, clear day, clear hour,clear minute, clear seconds
% clear thisDate, clear Mdate

if strcmp(Site,'2021-Young') || strcmp(Site,'2021-Hogg')
    % Changing Methane from umol to nmol. ONLY DO THIS WITH OLDER DATA
    columnToModify = data.FCH4_gf_RF;
    modifiedColumn = 1000 * columnToModify;
    data.FCH4_gf_RF = modifiedColumn;
    disp('Converted CH4: umol --> nmol')
    
    % Changing variable names
    data.Properties.VariableNames{'u_'} = 'ustar';
    data.Properties.VariableNames{'wind_speed'} = 'ubar';
    data.Properties.VariableNames{'wind_dir'} = 'WD';
    data.Properties.VariableNames{'hour'} = 'time';
    % data.Properties.VariableNames{'hour_dec'} = 'time'; %For BB2
    % data.Properties.VariableNames{'co2_flux'} = 'wc';
    data.Properties.VariableNames{'NEE_f'} = 'wc';
    % data.Properties.VariableNames{'ch4_flux'} = 'wm';
    data.Properties.VariableNames{'FCH4_gf_RF'} = 'wm';
    data.Properties.VariableNames{'h2o_flux'} = 'wq';
    data.Properties.VariableNames{'co2_molar_density'} = 'cbar';
    data.Properties.VariableNames{'ch4_molar_density'} = 'mbar';
    data.Properties.VariableNames{'h2o_molar_density'} = 'qbar';
    data.Properties.VariableNames{'v_var'} = 'vv';
    data.Properties.VariableNames{'DOY_x'} = 'DOY';
    data.Properties.VariableNames{'TA_1_1_1'} = 'TA';
    % data.Properties.VariableNames{'MO_LENGTH'} = 'L'; %Micromet lab already L
    % data.Properties.VariableNames{'air_temperature'} = 'TA'; %For BB2
end

if strcmp(Site,'US-WPT')
    % ::FOR US-WPT::
    data.Properties.VariableNames{'USTAR'} = 'ustar';
    data.Properties.VariableNames{'WS'} = 'ubar';
    data.Properties.VariableNames{'WD'} = 'WD';
    data.Properties.VariableNames{'FC'} = 'wc';
    data.Properties.VariableNames{'FCH4'} = 'wm';
    data.Properties.VariableNames{'CO2'} = 'cbar';
    data.Properties.VariableNames{'CH4'} = 'mbar';
    data.Properties.VariableNames{'H2O'} = 'qbar';
    data.Properties.VariableNames{'qbar'} = 'wq'; % Dataset doesn't have FH2O, so substituting H2O molar density
    data.Properties.VariableNames{'TA_1_1_1'} = 'TA';
    data.Properties.VariableNames{'MO_LENGTH'} = 'L'; %Micromet lab already L
    
    % Calculating v_var
    data.V_SIGMA(data.V_SIGMA==-9999)=nan;
    data = addvars(data,data.V_SIGMA.^2,'NewVariableNames','vv');
    data.V_SIGMA(isnan(data.V_SIGMA))=-9999;
    data.vv(isnan(data.vv))=-9999;  
    % Calculating DOY, year, month,hour
    doy = day(datetime(data.Mdate,'ConvertFrom','datenum'),'dayofyear');
    year = datetime(data.Mdate,'ConvertFrom','datenum').Year;
    month = datetime(data.Mdate,'ConvertFrom','datenum').Month;
    hour = datetime(data.Mdate,'ConvertFrom','datenum').Hour;
    data = addvars(data,doy,'NewVariableNames','DOY');
    data = addvars(data,year,'NewVariableNames','year');
    data = addvars(data,hour,'NewVariableNames','hour');
    time = datetime(data.Mdate,'ConvertFrom','datenum').Hour+(datetime(data.Mdate,'ConvertFrom','datenum').Minute)/60;
    data = addvars(data,time,'NewVariableNames','time');
    % Calculating daytime (where SWIN>0 = 1)
    daytime = ones(N,1);
    daytime(data.SW_IN==0)=0;
    data = addvars(data,daytime,'NewVariableNames','daytime');
    
    
end

data = standardizeMissing(data,-9999);

save([SiteName '_L3.mat'],'data');
save([SiteName '_met.mat'],'met');