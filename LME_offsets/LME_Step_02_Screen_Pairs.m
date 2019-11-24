addpath(genpath(pwd));

for yr = 1850:2014
    for mon = 1:12
        try
            % bucket compared with ERI - standard analysis  ---------------
            clear('P')
            P.yr                   = yr;
            P.mon                  = mon;
            P.do_connect           = 1;
            P.connect_Kobe         = 1;
            P.use_kent_tracks      = 0;
            P.use_diurnal_points   = 1;
            P.diurnal_QC           = 1;
            P.type                 = 'Bucket_vs_ERI';
            P.all_ERI_in_one_group = 1;
            P.save_app             = 'Bucket_vs_ERI_in_one_group';
            P.relative             = 'mean_SST';
            P.save_app             = [P.save_app,'_diurnal_points_',P.relative];
            LME_pair_02_Screen_Pairs(P);

        catch
            disp(['Year ',num2str(yr),' Mon',num2str(mon),' Failed...'])
        end

    end
end