clear all;
filenameList = [];

this_site='US-Myb';

% Young's L8 list
% DONE: 20210527,20210612,20210714,20210808,20210831,20210909,20210925,20211002,20211011
% 20211018,20211103,20220207,...
%     20220216,20220223,20220311,20220421,20220523,20220615,20220701,...
%     20220717,20220726,20220802,20220818,20220903,20220928,20221021,...
%     20221030,20221224,20230102,20230118,20230210,20230307,20230323!!

% MIA: 20211230,20220106,20220207,20220216,20220223

% Hogg's L8 List:
% DONE:20210511,20210527,20210612,20210714,20210730,20210831,20211002,20211018,
% 20220701,20220717,20220802,20220818,20220903
% 20221021,20221224,20220701,20220717,20220802,20220818,20220903,20221021
% MIA: 20230501,20211103,20220106,20220207,20220223,20220311,20221224

% US-WPT L8 List:
% [20130417 20130526 20130604 20130620 20130713 20130823 20130830 20130924
%  20131010]

% US-Myb L8 List:
% Done: 20130409, 20130416, 20130603, 20130619, 20130705,20130721,
% 20130822, 20130907, 20130923,20131212,20140318, 20140419,
% 20140505,20140521, 20140606, 20140622, 20140708, 20140724,
% 20140809,20140825, 20140910, 20140926, 20141012,20150305,
% 20150321,to20150508, 20150524,20150625, 20150711, 20150727, 20150812,
% 20150929, 20151031,20151116, 20151218, 20160204, 20160323, 20160510,
% 20160526,20160611, 20160627, 20160713, 20160729, 20160814,
% 20160830,20160915, 20161001,20170222,20170427, 20170513, 20170614,
% 20170630, 20170716, 20170801,20170817, 20170902, 20171004,
% 20171105,20180209,20180225, 20180329, 20180414, 20180601, 20180617,
% 20180719,20180804, 20180820, 20180905, 20180921, 20181007,
% 20181108,20181226, 20190228, 20190316, 20190417, 20190503,20190604,
% 20190620, 20190706, 20190722, 20190807, 20190823,2019090820200114,
% 20200403,20200419, 20200505, 20200521, 20200606, 20200622,
% 20200708,20200724, 20200809, 20200926, 20201012,

% MIA:
% 20131009,20131025,20131110,20131228,20140113,20141028,20141231,20150406,
% 20150422,20161017,20161102,20161118,20171207,20181124,20190924,20191010,
% 20191026,20201028,20201129,20210116,20210201,

L8dates = [20210217, 20210305, 20210321, 20210406,...
       20210422, 20210508, 20210524, 20210609, 20210625, 20210711,...
       20210828, 20210913, 20210929, 20211015];

% L8dates = [20130417 20130526 20130604 20130620 20130713 20130823 20130830 20130924 20131010];
    

% Creating FFP date ranges
for L8_date = L8dates
    cd /Users/darianng/Documents/MSc_Geography/Methane_Hotspot/FARF_Code
    datestr = num2str(L8_date);
    year = datestr(1:4);
    month = datestr(5:6);
    day = datestr(7:8);
    start_year = year;
    end_year = year;
    disp('Creating FFP start date');
%   Creating starting date
    if str2double(day) <= 15
        if strcmp(month,'01')
            start_month = '12';
            start_year = num2str(str2double(year)-1);
        else
            start_month = sprintf('%02d', str2double(month)-1); % w/ leading 0
        end
        if ismember(start_month,['01','03','05','7','8','10','12']) % 31 day months
            start_day = sprintf('%02d', 31-(15-str2double(day)));
        elseif strcmp(start_month,'02')
            start_day = sprintf('%02d', 28-(15-str2double(day)));
        else
            start_day = sprintf('%02d', 30-(15-str2double(day)));
        end
    else
        start_day = sprintf('%02d', str2double(day)-15);
        start_month = month;
    end
    disp('Creating FFP end date');
%   Creating ending date
    if ismember(month,['01','03','05','7','8','10','12']) % 31 day months
        if str2num(day)> 16
            end_day = sprintf('%02d', 15+str2num(day)-31);
            if strcmp(month,'12')
                end_year = num2str(str2num(year)+1);
                end_month = '01';
            else
                end_month = sprintf('%02d', str2num(month)+1);
            end
        else
            end_day =  sprintf('%02d',str2num(day)+15);
            end_month = month;
            
        end
    elseif strcmp(month,'02')
        end_month = sprintf('%02d',str2num(month)+1);
        if str2num(day)> 13
            end_day = sprintf('%02d', 15+str2num(day)-28);
        else
            end_day =  sprintf('%02d',str2num(day)+15);
        end
    else
        if str2num(day)> 15 % 30 day months
            end_day = sprintf('%02d', 15+str2num(day)-30);
            end_month = sprintf('%02d', str2num(month)+1);
        else
            end_day =  sprintf('%02d',str2num(day)+15);
            end_month = month;
        end
    end
        
    ffp_start = [start_year,start_month,start_day];
    ffp_end = [end_year,end_month,end_day];
    disp(['ffp start: ',ffp_start,', ffp end: ',ffp_end])

%     Starting Flux Map Calculations
    tstart = [start_year,'-',start_month,'-0',start_day,' 00:00'];
    tend = [end_year,'-',end_month,'-0',end_day,' 23:30'];
    disp(['Starting Flux Map Calculations for: ',ffp_start,'-',ffp_end,' || L8: ',datestr]);

    [filename] = Flux_Map_Func(this_site,tstart,tend);
    filenameList = cat(1,filenameList,filename);
    close;close;close;close;close;
end

