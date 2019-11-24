function DA_LME_Fig_Step_01_prepare_data(yr_start,yr_end,ct_reg,do_ERI,do_C0_SI_4,app_da,app_lme)

        
    %% *************************************************************************
    % Set directories
    % *************************************************************************
    % dir_save = DIURNAL_OI('Mis');
    % code to call this function

    % *************************************************************************
    % Set parameters
    % *************************************************************************
    do_figure = 1;

    % *********************************************************************
    % Step 1. Plot DA/DP versus LME offsets
    % *********************************************************************
    if do_figure == 1,
        
        P.relative = 'mean_SST';
        da_app     = '';

        dir_load  = DIURNAL_OI('data4figure');

        % Find LME and name parameters for individual regions
        % O.region_name         O.reg_list          O.sea_list
        O = DA_LME_function_get_region_name(ct_reg);

        % load bucket results -------------------------------------------------
        method    = 'bucket';
        file_load = [dir_load,'STATS_',method,'_DA_signals_',num2str(yr_start),...
            '_',num2str(yr_end),'_',O.region_name,da_app,...
            '_relative_to_',P.relative,'.mat'];
        disp(file_load)
        BCK = load(file_load);

        % load ERI results ----------------------------------------------------
        method     = 'ERI';
        if yr_start < 1921,
            file_load = [dir_load,'STATS_',method,'_DA_signals_',num2str(1970),...
                '_',num2str(1989),'_',O.region_name,da_app,...
                '_relative_to_',P.relative,'.mat'];
        else
            file_load = [dir_load,'STATS_',method,'_DA_signals_',num2str(yr_start),...
                '_',num2str(yr_end),'_',O.region_name,da_app,...
                '_relative_to_',P.relative,'.mat'];
        end
        disp(file_load)
        ERI = load(file_load);

        % LME -----------------------------------------------------------------
        if strcmp(P.relative,'mean_SST');
            dir_lme = LME_OI('home');
            if do_C0_SI_4 == 1,
                dir_lme = [dir_lme,'Step_04_LME_output/'];
            else
                dir_lme = [dir_lme,'Step_04_LME_output_v6/'];
            end
            yr_end_2 = yr_end;          
            if do_ERI == 0,
                file_lme = [dir_lme,'LME_Bucket_vs_Bucket_diurnal_points_mean_SST_',num2str(yr_start),'_',num2str(yr_end_2),'_Full_SST_Global.mat'];
            else
                file_lme = [dir_lme,'LME_Bucket_vs_ERI_in_one_group_',app_lme,'_',num2str(yr_start),'_',num2str(yr_end_2),'_Full_SST_Global.mat'];
            end
            disp(file_lme)
            LME = load(file_lme,'out','out_rnd');
        end

        disp('Read data finished!')
        
        % ****************************************************************
        % Put files into a good shape
        % *****************************************************************
        % Diurnal amplitude and phase -------------------------------------
        clear('DATA')
        key = 500;
        clear('a');
        for i = 1:4      a(i,:) = nansum(BCK.num_grp([1:6]+(i-1)*6,1,:),1);      end
        l_qc  = all(a > key,1);
        da_group = double(BCK.groups(l_qc,:));
        fit_out_grp     = BCK.fit_out_grp(:,:,l_qc);
        fit_out_std_grp = BCK.fit_out_std_grp(:,:,l_qc);
        num_bck_DA      = BCK.num_grp(:,:,l_qc);

        % Construct diurnal signals ---------------------------------------
        C0_LCL_hr = [1:24];
        omega = 2*pi/24;
        base_x2 = [ones(numel(C0_LCL_hr),1) sin(C0_LCL_hr'*omega) cos(C0_LCL_hr'*omega)  sin(C0_LCL_hr'*omega*2) cos(C0_LCL_hr'*omega*2)];
        
        clear('diurnal')
        for ct = 1:size(fit_out_grp,3)
            for ct_sea = 1:3
                diurnal(:,ct_sea,ct) = base_x2 * fit_out_grp(4:8,ct_sea,ct);
            end
        end

        % Prepare for data to be plotted: diurnal -------------------------
        clear('aa')
        for i = 1:1000
            aa(:,:,:,i) = normrnd(fit_out_grp,fit_out_std_grp);
        end
        
        clear('amp','amp_rnd','amp_std')
        clear('phs','phs_rnd','phs_std')
        
        [phs,amp]  = cart2pol(fit_out_grp(6,:,:),fit_out_grp(5,:,:));
        phs(phs<0) = (phs(phs<0) + 2*pi);
        phs        = phs / pi * 12;

        [phs_rnd,amp_rnd]  = cart2pol(aa(6,:,:,:),aa(5,:,:,:));
        phs_rnd(phs_rnd<0) = (phs_rnd(phs_rnd<0) + 2*pi);
        phs_rnd    = phs_rnd / pi * 12;

        if ~strcmp(app_da,'raw'),
            l = ismember(da_group,[67 78 781],'rows');
            ll = (phs<9 | phs> 21) & ~repmat(reshape(l,1,1,numel(l)),1,3,1);
            amp(ll) = -amp(ll);
            ll = (phs_rnd<9 | phs_rnd> 21) & ~repmat(reshape(l,1,1,numel(l)),1,3,1,1000);
            amp_rnd(ll) = -amp_rnd(ll);
        end

        amp_std = CDC_std(amp_rnd,4);
        phs_std = CDC_std(phs_rnd,4);
        
        amp = squeeze(amp)';
        phs = squeeze(phs)';
        amp_std = squeeze(amp_std)';
        phs_std = squeeze(phs_std)';
        
        
        % LME offsets -----------------------------------------------------
        lme_group = LME.out.unique_grp(:,1:3);

        offset_fixed  = LME.out.bias_fixed; 
        offset_region = nanmean(LME.out.bias_region(O.reg_list,:),1)';
        offset_season = LME.out.bias_season(O.sea_list,:)';

        offset_fixed_rnd  = LME.out_rnd.bias_fixed_random; 
        offset_region_rnd = squeeze(nanmean(LME.out_rnd.bias_region_rnd(O.reg_list,:,:),1))';
        offset_season_rnd = LME.out_rnd.bias_season_rnd(O.sea_list,:,:);
        
        clear('offset')
        offset(:,1) = offset_fixed + offset_region;
        offset(:,2) = offset(:,1)  + offset_season(:,1);
        offset(:,3) = offset(:,1)  + offset_season(:,2);

        clear('offset_std')
        offset_std(:,1) = CDC_std(offset_fixed_rnd + offset_region_rnd,1);
        offset_std(:,2) = CDC_std(offset_fixed_rnd + offset_region_rnd + squeeze(offset_season_rnd(1,:,:))',1);
        offset_std(:,3) = CDC_std(offset_fixed_rnd + offset_region_rnd + squeeze(offset_season_rnd(2,:,:))',1);
        
        
        % Find common groups ----------------------------------------------
        clear('DATA_pic')
        [group,I_DA,I_LME]  = intersect(da_group,lme_group,'rows');

        DATA_pic.DA_TS_bck   = diurnal(:,:,I_DA);
        DATA_pic.DA_num_bck = num_bck_DA(:,:,I_DA);
        DATA_pic.DA_fit_bck = fit_out_grp(4:8,:,I_DA);
        DATA_pic.DA_fit_std_bck = fit_out_std_grp(4:8,:,I_DA);
        DATA_pic.DA_amp_bck = amp(I_DA,:);
        DATA_pic.DA_amp_std_bck = amp_std(I_DA,:);
        DATA_pic.DA_phase_bck = phs(I_DA,:);
        DATA_pic.DA_phase_std_bck = phs_std(I_DA,:);
        DATA_pic.LME_bck     = offset(I_LME,:);
        DATA_pic.LME_std_bck = offset_std(I_LME,:); 
        DATA_pic.group_bck   = group;

        clear('amp','amp_rnd','amp_std')
        clear('phs','phs_rnd','phs_std')
        
        if do_ERI == 1,
            if all(LME.out.unique_grp(1,:) == [1 1 1 1]);
                DATA_pic.LME_ERI = offset(1,:);
                DATA_pic.LME_std_ERI = offset_std(1,:);
            else
                DATA_pic.LME_ERI = [nan nan nan];
                DATA_pic.LME_std_ERI = [nan nan nan];
            end
        end
        

        % Prepare for ERI and Buoy ---------------------------------------- 
        for ct = 1:2
            
            f_eri = ERI.fit_out_adj(:,:,ct);
            clear('diurnal_eri')
            for ct_sea = 1:3
                switch ct_sea,
                    case 1,
                        diurnal_eri(:,ct_sea) = base_x2 * nanmean(f_eri(4:8,:),2);
                    case 2,
                        diurnal_eri(:,ct_sea) = base_x2 * nanmean(f_eri(4:8,[1 2 12]),2);
                    case 3,
                        diurnal_eri(:,ct_sea) = base_x2 * nanmean(f_eri(4:8,[6 7 8]),2);
                end
            end
            
            if ct == 1,
                DATA_pic.DA_TS_ERI  = diurnal_eri;
            else
                DATA_pic.DA_TS_buoy = diurnal_eri;
            end

            
            clear('temp_save1','temp_save2','temp_save3')
            temp_save1 = [nansum(ERI.num_adj(:,:,ct),2)  ...
                                   nansum(ERI.num_adj(:,[1 2 12],ct),2)  ...
                                   nansum(ERI.num_adj(:,[6 7 8],ct),2)];

            temp_save2 = [nanmean(ERI.fit_out_adj(4:8,:,ct),2) ...
                                   nanmean(ERI.fit_out_adj(4:8,[1 2 12],ct),2) ...
                                   nanmean(ERI.fit_out_adj(4:8,[6 7 8],ct),2)];

            temp_save3 = sqrt([nanmean(ERI.fit_out_std_adj(4:8,[1:12],ct).^2,2)/12 ...
                                            nanmean(ERI.fit_out_std_adj(4:8,[1 2 12],ct).^2,2)/3 ...
                                            nanmean(ERI.fit_out_std_adj(4:8,[6 7 8],ct).^2,2)/3]);
                                        
            if ct == 1,
                DATA_pic.DA_num_ERI     = temp_save1;
                DATA_pic.DA_fit_ERI     = temp_save2;
                DATA_pic.DA_fit_std_ERI = temp_save3;
            else
                DATA_pic.DA_num_buoy     = temp_save1;
                DATA_pic.DA_fit_buoy     = temp_save2;
                DATA_pic.DA_fit_std_buoy = temp_save3;                
            end

            clear('aa')
            for i = 1:1000
                aa(:,:,i) = normrnd(temp_save2,temp_save3);
            end

            clear('amp','amp_rnd','amp_std')
            clear('phs','phs_rnd','phs_std')

            [phs,amp]  = cart2pol(temp_save2(3,:,:),temp_save2(2,:,:));
            phs(phs<0) = (phs(phs<0) + 2*pi);
            phs        = phs / pi * 12;

            [phs_rnd,amp_rnd]  = cart2pol(aa(3,:,:),aa(2,:,:));
            phs_rnd(phs_rnd<0) = (phs_rnd(phs_rnd<0) + 2*pi);
            phs_rnd    = phs_rnd / pi * 12;

            if ~strcmp(app_da,'raw'),
                amp(phs<9 | phs> 21) = -amp(phs<9 | phs> 21);
                amp_rnd(phs_rnd<9 | phs_rnd> 21) = -amp_rnd(phs_rnd<9 | phs_rnd> 21);
            end

            amp_std = CDC_std(amp_rnd,3);
            phs_std = CDC_std(phs_rnd,3);

            if ct == 1,
                DATA_pic.DA_amp_ERI       = amp;
                DATA_pic.DA_amp_std_ERI   = amp_std;
                DATA_pic.DA_phase_ERI     = phs;
                DATA_pic.DA_phase_std_ERI = phs_std;
            else
                DATA_pic.DA_amp_buoy       = amp;
                DATA_pic.DA_amp_std_buoy   = amp_std;
                DATA_pic.DA_phase_buoy     = phs;
                DATA_pic.DA_phase_std_buoy = phs_std;              
            end

            clear('amp','amp_rnd','amp_std')
            clear('phs','phs_rnd','phs_std')
        
        end
        
        % *****************************************************************
        % Save data 1
        % *****************************************************************
        dir_save  = LME_OI('data4figure');
        file_save = [dir_save,'DATA_Plot_DA&LME_',num2str(yr_start),...
            '_',num2str(yr_end),'_',O.region_name,'_da_',app_da,'_lme_',app_lme,'.mat'];
        save(file_save,'DATA_pic','-v7.3');

        % *****************************************************************
        % Data for all buckets
        % *****************************************************************
        clear('DATA_pic')
        DATA_pic.diurnal(:,1) = nanmean(BCK.diurnal_adj(:,1:12,1),2);
        DATA_pic.diurnal_quantile(:,:,1) = nanmean(BCK.diurnal_quantile_adj(:,:,1:12,1),3);
        DATA_pic.num(:,1) = nansum(BCK.num_adj(:,1:12,1),2);

        DATA_pic.diurnal(:,2) = nanmean(BCK.diurnal_adj(:,[1 2 12],1),2);
        DATA_pic.diurnal_quantile(:,:,2) = nanmean(BCK.diurnal_quantile_adj(:,:,[1 2 12],1),3);
        DATA_pic.num(:,2) = nansum(BCK.num_adj(:,[1 2 12],1),2);

        DATA_pic.diurnal(:,3) = nanmean(BCK.diurnal_adj(:,[6 7 8],1),2);
        DATA_pic.diurnal_quantile(:,:,3) = nanmean(BCK.diurnal_quantile_adj(:,:,[6 7 8],1),3);
        DATA_pic.num(:,3) = nansum(BCK.num_adj(:,[6 7 8],1),2);

        DATA_pic.diurnal_buoy(:,1) = nanmean(BCK.diurnal_adj(:,1:12,2),2);
        DATA_pic.diurnal_buoy(:,2) = nanmean(BCK.diurnal_adj(:,1:12,3),2);
        
        % *****************************************************************
        % Save data 2
        % *****************************************************************
        dir_save  = LME_OI('data4figure');
        file_save = [dir_save,'DATA_Plot_All_buckets_',num2str(yr_start),...
            '_',num2str(yr_end),'_',O.region_name,'_da_',app_da,'_lme_',app_lme,'.mat'];
        save(file_save,'DATA_pic','-v7.3');
    end
end