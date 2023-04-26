function [FileToLoad]=GetSiteInfo(site,PLT)

%% East End
if strcmp(site,'BB1')
    
if PLT==1
FileToLoad= 'BB1_output';   
end
if PLT==2
FileToLoad= 'EE-K&M-2017250-0000-to-2017290-2400-day-50-80v2';
end
if PLT==3
FileToLoad= 'EE-K&M-2017130-0000-to-2017275-2400-day-50-80v2';
end
if PLT==4
FileToLoad= 'EE-K&M-2018100-0000-to-2018300-2400-day-50-80_T_20to25_C';
end
if PLT==5
FileToLoad= 'EE-K&M-2017100-0000-to-2017300-2400-day-50-80_T_20to25_C2';
end

elseif strcmp(site,'Hogg')
    
if PLT==1
FileToLoad= 'Hogg_output';   
end

elseif strcmp(site,'Young')
    
if PLT==1
FileToLoad= 'Young_output';   
end

elseif strcmp(site,'BB2')
    
if PLT==1
FileToLoad= 'BB2_output';   
end

elseif strcmp(site,'US-Myb')
    
if PLT==1
FileToLoad= 'US-Myb_output';   
end

elseif strcmp(site,'US-WPT')
    
if PLT==1
FileToLoad= 'US-WPT_output';   
end

end

end

