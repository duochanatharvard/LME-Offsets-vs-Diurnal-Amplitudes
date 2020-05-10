function DA_POST_sum_and_fitting(do_sum,yr_start,yr_end,method,do_nation)

    if 0,
        % how to call this function
        % Step 1: summing up all measurements
        DA_POST_sum_and_fitting(1,[],[],'bucket');
        DA_POST_sum_and_fitting(1,[],[],'ERI');
        
        % Step 2: loop over years and fit for individual groups
        for yr_start = 1880:1990
            yr_end   = yr_start + 19;
            method   = 'bucket';
            do_nation  = 0;
            DA_POST_sum_and_fitting(0,yr_start,yr_end,method,do_nation);
            do_nation  = 1;
            DA_POST_sum_and_fitting(0,yr_start,yr_end,method,do_nation);
        end

        for yr_start = 1920:1990
            yr_end   = yr_start + 19;
            method   = 'ERI';
            do_nation  = 0;
            DA_POST_sum_and_fitting(0,yr_start,yr_end,method,do_nation);
            do_nation  = 1;
            DA_POST_sum_and_fitting(0,yr_start,yr_end,method,do_nation);
        end

    end
    
    % *************************************************************************
    % Set directories
    % *************************************************************************
    dir_save  = DIURNAL_OI('data4figure');
    % dir_save = '/Users/duochan/Data/DIURNAL_2019/DATA_for_figures_C0_SI_2/';
    
    % *************************************************************************
    % Set parameters
    % *************************************************************************
    season = 'Annual';                  
    mon_list = [1:12];

    if strcmp(method,'bucket')
        yr_list = 1850:2009;            % for step 1
    else
        yr_list = 1900:2009;            % for step 1
    end

    P.relative = 'mean_SST';            % for step 1 and 2
    P.decimal  = 'normal';              % not changed
    dir_load   = DIURNAL_OI('ship_signal');
    da_app     = '';

    % *************************************************************************
    % Step 1. Prepare for data
    % *************************************************************************
    if do_sum == 1
        % *********************************************************************
        % Read data of diurnal signals from individual ships
        % *********************************************************************
        file_save = [dir_save,'SUM_',method,'_DA_signals_',num2str(yr_list(1)),...
            '_',num2str(yr_list(end)),'_',season,da_app,...
            '_relative_to_',P.relative,'.mat'];
        fid = fopen(file_save);

        if fid <= 0
            var_list = {'C0_LAT','C0_LON','C0_SI_4','D1_EXP','D2_EXP',...
                'D1_EX','D2_EX','C1_DCK','C0_YR','C0_MO','C0_LCL','C98_UID',...
                'Diurnal_signal','Fundemental_SST','Day_indicator','C0_CTY_CRT',...
                'C0_SI_2','C0_SI_3'};

            for var = 1:numel(var_list)
                eval([var_list{var},' = [];']);
            end

            for mon = mon_list
                % *****************************************************************
                % Load in data
                % *****************************************************************
                for yr = yr_list

                    disp(['Starting year ',num2str(yr),'  Month ',num2str(mon)])

                    file_load = [dir_load,'IMMA1_R3.0.0_',num2str(yr),...
                        '-',CDF_num2str(mon,2),'_Ship_Diurnal_Signal',da_app,...
                        '_relative_to_',P.relative,'.mat'];

                    file_icoads3 = [DIURNAL_OI('read_raw'),'ICOADS_QCed/',...
                        'IMMA1_R3.0.0_',num2str(yr),'-',CDF_num2str(mon,2),'_QCed.mat'];

                    try
                        temp = load(file_load,...
                            'C0_LAT','C0_LON','C0_SI_4','D1_EXP','D2_EXP',...
                            'D1_EX','D2_EX','C1_DCK','C0_YR','C0_MO','C0_LCL',...
                            'Diurnal_signal','Fundemental_SST','Day_indicator',...
                            'C0_CTY_CRT','C98_UID');
                        
                        temp_icoads3 = load(file_icoads3,'C0_SI_4','C0_SI_2','C0_SI_3','C98_UID');
                        [~,pst] = ismember(temp.C98_UID,temp_icoads3.C98_UID);
                        temp.C0_SI_2 = temp_icoads3.C0_SI_2(pst);
                        temp.C0_SI_3 = temp_icoads3.C0_SI_3(pst);

                        if strcmp(method,'bucket')
                            l = temp.C0_SI_4 >= 0 &  temp.C0_SI_4 <= 0.05;
                        elseif strcmp(method,'ERI')
                            l = temp.C0_SI_4 >= 0.95 &  temp.C0_SI_4 <= 1;
                        end

                        for var = 1:numel(var_list)
                            if ~ismember(var_list{var},{'C0_ID','C0_CTY_CRT','DCK'})
                                eval([var_list{var},' = [',var_list{var},'  temp.',var_list{var},'(l)];']);
                            else
                                eval([var_list{var},' = [',var_list{var},'; temp.',var_list{var},'(l,:)];']);
                            end
                        end
                    catch
                        disp([file_load,' does not exist'])
                    end
                end
            end

            save(file_save,...
                'C0_LAT','C0_LON','C0_SI_4','C0_SI_2','C0_SI_3','D1_EXP','D2_EXP',...
                'D1_EX','D2_EX','C1_DCK','C0_YR','C0_MO','C0_LCL',...
                'Diurnal_signal','Fundemental_SST','Day_indicator','C0_CTY_CRT','-v7.3')
        else
            fclose(fid);
            disp('file exist!')
        end

    else
        
        % *****************************************************************
        % Step 2. Compute statistics
        % *****************************************************************        
        file_load = [dir_save,'SUM_',method,'_DA_signals_',num2str(yr_list(1)),...
            '_',num2str(yr_list(end)),'_',season,da_app,...
            '_relative_to_',P.relative,'.mat'];
        DATA = load(file_load);

        use_C0_SI_2 = 0;                        % TODO

        l = DATA.C0_YR >= yr_start & DATA.C0_YR <= yr_end & ...
            ismember(DATA.Day_indicator,[1]);
        
        if use_C0_SI_2 == 1
            l = l & ismember(DATA.C0_SI_2,[0 1 -2 3]);
        end
        
        % yr      = DATA.C0_YR(l);
        lat     = DATA.C0_LAT(l);
        lon     = DATA.C0_LON(l);
        month   = DATA.C0_MO(l);
        di_sig  = DATA.Diurnal_signal(l);
        lcl     = DATA.C0_LCL(l);
        dck     = DATA.C1_DCK(l);
        nat     = DATA.C0_CTY_CRT(l,:);
        exp1    = DATA.D1_EXP(l);
        exp2    = DATA.D2_EXP(l);
        clear('l','DATA')
        
        DA_relative_to = 2;            % compute the absolute diurnal amplitude         % TODO
        if DA_relative_to == 2         % standard case: buoy chan 2019
            di_sig = exp1;
        elseif  DA_relative_to == 3    % standard case: buoy MB16
            di_sig = di_sig - exp2;
        elseif  DA_relative_to == 4    % standard case: buoy MB16
            di_sig = exp1;
        end

        mon_adj             = month;
        mon_adj(lat<0)      = mon_adj(lat<0) + 6;
        mon_adj(mon_adj>12) = mon_adj(mon_adj>12) - 12;

        region = LME_lme_effect_regional(lon,lat,5);

        grp           = [nat, dck'];
        % connect decks !!! and rerun the code ...
        P_dck.do_connect   = 1;
        P_dck.connect_Kobe = 1;
        grp(:,1:3) = LME_function_preprocess_deck(double(grp(:,1:3)),P_dck);
        
        if do_nation == 1
            grp = grp(:,1:2);
        end
        
        [grp_uni,~,J] = unique(grp,'rows');
        key           = 1000;
        c             = hist(J,1:1:max(J));
        l_use_grp     = find(c > key);
        groups        = grp_uni(l_use_grp,:);


        for ct_reg = [1 3 7]            % TODO

            O = DA_LME_function_get_region_name(ct_reg);
            l = ismember(region,O.reg_list);

            clear('diurnal','diurnal_quantile','num','fit_out','fit_out_std')
            clear('diurnal_adj','diurnal_quantile_adj','num_adj','fit_out_adj','fit_out_std_adj')
            for ct_mon = 1:12

                % compute diurnal signal for all groups pulled together .......
                for ct_data = 1:3
                    switch ct_data,
                        case 1,
                            data_in = di_sig;
                        case 2,
                            data_in = exp1;
                        case 3,
                            data_in = exp2;
                    end

                    l_use = l & month == ct_mon;
                    [diurnal(:,ct_mon,ct_data),diurnal_quantile(:,:,ct_mon,ct_data),num(:,ct_mon,ct_data),...
                        fit_out(:,ct_mon,ct_data),fit_out_std(:,ct_mon,ct_data)] = ...
                        Compute_and_fit_diurnal_signal(lcl(l_use),data_in(l_use));

                    l_use = l & mon_adj == ct_mon;
                    [diurnal_adj(:,ct_mon,ct_data),diurnal_quantile_adj(:,:,ct_mon,ct_data),num_adj(:,ct_mon,ct_data),...
                        fit_out_adj(:,ct_mon,ct_data),fit_out_std_adj(:,ct_mon,ct_data)] = ...
                        Compute_and_fit_diurnal_signal(lcl(l_use),data_in(l_use));
                end
            end


            % compute diurnal signal for individual groups - annual ...........
            clear('diurnal_grp','diurnal_quantile_grp','num_grp','fit_out_grp','fit_out_std_grp')
            ct_grp = 0;
            for ct = l_use_grp
                ct_grp = ct_grp + 1;

                for ct_sea = 1:3
                    switch ct_sea,
                        case 1,
                            l_use = l & J' == ct;
                        case 2,
                            l_use = l & J' == ct & ismember(mon_adj,[1 2 12]);
                        case 3,
                            l_use = l & J' == ct & ismember(mon_adj,[6 7 8]);
                    end

                    if nnz(l_use) > 10,
                        [diurnal_grp(:,ct_sea,ct_grp),diurnal_quantile_grp(:,:,ct_sea,ct_grp),num_grp(:,ct_sea,ct_grp),...
                            fit_out_grp(:,ct_sea,ct_grp),fit_out_std_grp(:,ct_sea,ct_grp)] = ...
                            Compute_and_fit_diurnal_signal(lcl(l_use),di_sig(l_use));
                    end
                end
            end

            % Save data
            if do_nation == 1
                file_save = [dir_save,'STATS_',method,'_DA_signals_nation_level_',num2str(yr_start),...
                    '_',num2str(yr_end),'_',O.region_name,da_app,...
                    '_relative_to_',P.relative,'.mat'];
            else
                file_save = [dir_save,'STATS_',method,'_DA_signals_',num2str(yr_start),...
                    '_',num2str(yr_end),'_',O.region_name,da_app,...
                    '_relative_to_',P.relative,'.mat'];
            end
            save(file_save,'diurnal','diurnal_quantile','num','fit_out','fit_out_std',...
                'diurnal_adj','diurnal_quantile_adj','num_adj','fit_out_adj','fit_out_std_adj',...
                'diurnal_grp','diurnal_quantile_grp','num_grp','fit_out_grp','fit_out_std_grp','groups','-v7.3')

        end
    end
end