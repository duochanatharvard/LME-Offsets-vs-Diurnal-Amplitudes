addpath(genpath(pwd));

for yr = 1850:2014
    for mon = 1:12
        try
            P.yr = yr;
            P.mon = mon;
            LME_pair_01_Raw_Pairs(P);
        catch
            disp('Year ',num2str(yr),' Mon',num2str(mon),' Failed...')
        end
    end
end