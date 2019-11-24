% Compute diurnal signals from ship SSTs
% 
% This version uses the typical daynight conept for mid latitude, 
% 
% P.relative: diurnal signal is relative to which SST?
% 1. mean_SST  2. fund_SST

function DIURNAL_Step_01_ship_diurnal_signal(yr,P)

    if 0,
        % To run this code, use the following line
        P.relative = 'mean_SST';
        for yr = 1880:1899
            DIURNAL_Step_01_ship_diurnal_signal(yr,P);
        end
    end

    % *********************************************************************
    % Set directories
    % *********************************************************************
    dir_load = DIURNAL_OI('ICOADS3');
    dir_kent = DIURNAL_OI('kent_track');
    dir_save = DIURNAL_OI('ship_signal');
    dir_mis  = DIURNAL_OI('Mis');
    
    % *********************************************************************
    % Set parameters
    % *********************************************************************
    spd_threshold = 100;
    
    % *********************************************************************
    % Load the data and subset for ship measurements with call signs
    % *********************************************************************
    disp('Load Data ...')
    var_list = {'C0_YR','C0_MO','C0_DY','C0_HR','C0_LCL','C0_UTC',...
            'C0_LON','C0_LAT','C0_SST','C0_SI_4','C0_AT','C0_II',...
            'C1_DCK','C98_UID','C0_OI_CLIM',...
            'C0_W','C0_N','QC_FINAL','C0_CTY_CRT','C0_ID'};
    for var = 1:numel(var_list)
        eval([var_list{var},' = [];']);
    end
    clear('var')
    for mon = 1:12
        
        clear('temp','temp_kent','file_load','file_kent')
        file_load = [dir_load,'IMMA1_R3.0.0_',num2str(yr),...
                     '-',CDF_num2str(mon,2),'_QCed.mat'];

        file_kent = [dir_kent,'IMMA1_R3.0.0_',num2str(yr),...
                     '-',CDF_num2str(mon,2),'_Tracks_Kent.mat'];
                 
        temp = load(file_load,...
            'C0_YR','C0_MO','C0_DY','C0_HR','C0_LCL','C0_UTC',...
            'C0_LON','C0_LAT','C0_SST','C0_SI_4','C0_AT','C0_II',...
            'C1_DCK','C98_UID','C0_W','C0_N',...
            'QC_FINAL','C0_CTY_CRT','C0_OI_CLIM');
        
        temp_kent = load(file_kent);

        clear('l_empty','l','l_NA')
        l_empty = all(temp_kent.C0_ID_K == 32,2)' | isnan(temp_kent.C98_UID);
        l_NA    = all(temp_kent.C0_ID_K == repmat(['NA',repmat(' ',1,28)],numel(temp_kent.C98_UID),1),2);
        l       = ~ismember(temp.C0_SI_4,[-3 -2]) & ~l_NA' & ~l_empty & temp.QC_FINAL == 1;
                
        for var = 1:numel(var_list)
            if ~ismember(var_list{var},{'C0_ID','C0_CTY_CRT'}),
                eval([var_list{var},' = [',var_list{var},'  temp.',var_list{var},'(1,l)];']);
            elseif ismember(var_list{var},{'C0_CTY_CRT'}),
                eval([var_list{var},' = [',var_list{var},';  temp.',var_list{var},'(l,:)];']);
            elseif ismember(var_list{var},{'C0_ID'}),
                C0_ID = [C0_ID; temp_kent.C0_ID_K(l,:)];
            end
        end
        
        clear('temp','temp_kent','file_load','file_kent')
        clear('l_empty','l','l_NA')
    end
    clear('mon')
    
    C0_N(C0_N == 9) = nan;
    C0_N = C0_N / 8;
    
    % *********************************************************************
    % Load information of diurnal cycles estimated from buoy
    % *********************************************************************
    clear('Diurnal_Shape','Diurnal_Amplitude','file_amplitude','file_shape')
    file_shape =  [dir_mis,'Diurnal_Shape_buoy_SST.mat'];
    Diurnal_Shape = load(file_shape,'Diurnal_Shape');
    file_amplitude = [dir_mis,'Diurnal_Amplitude_buoy_SST_1990_2014_climatology.mat'];
    Diurnal_Amplitude = load(file_amplitude,'Diurnal_clim_buoy_1990_2014');
    clear('file_amplitude','file_shape')

    % *********************************************************************
    % Load NOCS data for infilling wind and cloud
    % *********************************************************************
    clear('file_wind','file_cloud')
    file_wind  = [dir_mis,'NOCS_5X5_wspd_1973_2002.mat'];
    file_cloud = [dir_mis,'NOCS_5X5_cldc_1973_2002.mat'];
    NOCS_wind  = load(file_wind);
    NOCS_cloud = load(file_cloud);
    clear('file_wind','file_cloud')
    
    clear('x_all','y_all','m_all','id_all')
    m_all = C0_MO;
    x_all = discretize(C0_LON,0:5:360);
    y_all = discretize(C0_LAT,-90:5:90);
    id_all = sub2ind([72,36,12],x_all,y_all,m_all);
    NOCS_wind = NOCS_wind.clim_final(id_all);
    NOCS_cloud = NOCS_cloud.clim_final(id_all)/100;
    clear('x_all','y_all','m_all','id_all')
    
    % *********************************************************************
    % Create variables for results to be saved
    % *********************************************************************
    for var = 1:numel(var_list)-2
        eval([var_list{var},'_save = [ ];']);
    end 
    clear('var')
    C0_CTY_CRT_save = [];
    FD_save = [];
    DA_save = [];
    D1_EXP_save = [];
    D2_EXP_save = [];
    D1_EX_save  = [];
    D2_EX_save  = [];
    Ship_save   = [];
    DAY_save    = [];
    DAY_indicator_save = [];

    clear('DATA_save')
    DATA_save.m            = [];
    DATA_save.lon_m        = [];
    DATA_save.lat_m        = [];
    DATA_save.utc_m        = [];
    DATA_save.fit_ex1      = [];
    DATA_save.fit_ex2      = [];
    DATA_save.fit_raw1     = [];
    DATA_save.fit_raw2     = [];
    DATA_save.fit_ex1_int  = [];
    DATA_save.fit_ex2_int  = [];
    DATA_save.fit_raw1_int = [];
    DATA_save.fit_raw2_int = [];
    DATA_save.deck         = [];
    DATA_save.country      = [];
    DATA_save.sample_size  = [];
    DATA_save.first_uid    = [];
    DATA_save.ship_id      = [];
    DATA_save.day_id       = [];
    DATA_save.qc           = [];
    DATA_save.SST_method   = [];

    % *********************************************************************
    % Loop over individual IDs to compute the diurnal signal
    % *********************************************************************
    clear('ID_uni','J_ID')
    [ID_uni,~,J_ID]=unique(C0_ID,'rows');
    disp(['A total of ',num2str(size(ID_uni,1)),' unique ship IDs in ',num2str(yr)])
     
    disp('Start Computing Diurnal Signals ...')
    for ct_id = 1:1:size(ID_uni,1)                    % TODO: 500 -> 1

        if rem(ct_id,100) == 0,
            disp(['Starting the ',num2str(ct_id),'th track ... ']);
        end

        % *****************************************************************
        % Subset data
        % *****************************************************************
        clear('l','C0_CTY_CRT_id','NOCS_wind_id','NOCS_cloud_id')
        l = J_ID == ct_id;
        for var = 1:numel(var_list)-2
            eval(['clear(''',var_list{var},'_id'')'])
            eval([var_list{var},'_id = ',var_list{var},'(l);']);
        end
        clear('var')
        C0_CTY_CRT_id = C0_CTY_CRT(l,:);
        NOCS_wind_id  = NOCS_wind(l);
        NOCS_cloud_id = NOCS_cloud(l);

        % *****************************************************************
        % Track check
        % *****************************************************************
        % 1. through away duplicate points
        clear('L_dup')
        [~,L_dup] = unique(C0_UTC_id);
        for var = 1:numel(var_list)-2
            eval([var_list{var},'_id = ',var_list{var},'_id(L_dup);']);
        end 
        NOCS_wind_id  = NOCS_wind_id(L_dup);
        NOCS_cloud_id = NOCS_cloud_id(L_dup);
        C0_CTY_CRT_id = C0_CTY_CRT_id(L_dup,:);
        clear('L_dup','var')
        
        % 2. check speed
        clear('dis','spd')
        dis = distance(C0_LAT_id(2:end),C0_LON_id(2:end),...
            C0_LAT_id(1:end-1),C0_LON_id(1:end-1));
        spd = dis*111./(C0_UTC_id(2:end) - C0_UTC_id(1:end-1)); 
        clear('dis')
        
        % *****************************************************************
        % Compute diurnal signal, only if the buoy speed is less than 
        % a threshold in unit of km/hr
        % *****************************************************************        
        if quantile(spd,0.95) <= spd_threshold,
        
            % *************************************************************
            % Loop over individual days to compute diurnal cycles
            % *************************************************************
            clear('day_list','day_uni')
            day_list = datenum([C0_YR_id; C0_MO_id; C0_DY_id]') - ...
                       datenum([C0_YR_id; C0_MO_id*0+1; C0_DY_id*0+1]') + 1;
            day_uni  = unique(day_list);

            for ct_dy = reshape(day_uni,1,numel(day_uni))

                clear('l','NOCS_wind_dy','NOCS_cloud_dy','C0_CTY_CRT_dy')
                l = ismember(day_list,[-1 0 1]+ct_dy);
                for var = 1:numel(var_list)-2
                    eval(['clear(''',var_list{var},'_dy'')'])
                    eval([var_list{var},'_dy = ',var_list{var},'_id(l);']);
                end
                NOCS_wind_dy = NOCS_wind_id(l);
                NOCS_cloud_dy = NOCS_cloud_id(l);
                C0_CTY_CRT_dy = C0_CTY_CRT_id(l,:);
                clear('var','l')
                
                % *********************************************************
                % Truncate for two nights and a day in between
                % *********************************************************
                clear('hr_temp')
                hr_temp = round(C0_UTC_dy - C0_UTC_dy(1)) + C0_LCL_dy(1);
                if(min(hr_temp) > 20)
                    hr_temp = hr_temp - 24;
                end
                if(any(hr_temp > 24 & hr_temp < 56,2)==0)
                    hr_temp = hr_temp + 24;
                end
                if(any(hr_temp > 24 & hr_temp < 56,2)==0)
                    hr_temp = hr_temp - 48;
                end

                clear('l','NOCS_wind_hr','NOCS_cloud_hr','C0_CTY_CRT_hr') 
                
                if strcmp(P.relative,'fund_SST');
                    l = hr_temp >= 24 & hr_temp <= 54;
                else
                    l = hr_temp >= 24 & hr_temp < 48;
                end
                
                for var = 1:numel(var_list)-2
                    eval(['clear(''',var_list{var},'_hr'')'])
                    eval([var_list{var},'_hr = ',var_list{var},'_dy(l);']);
                end
                NOCS_wind_hr = NOCS_wind_dy(l);
                NOCS_cloud_hr = NOCS_cloud_dy(l);
                C0_CTY_CRT_hr = C0_CTY_CRT_dy(l,:);
                hr_temp = hr_temp(l);
                clear('var','l')

                % *********************************************************
                % Check if that day has at least one measurement 
                % in each 6-hourly section
                % *********************************************************
                if strcmp(P.relative,'fund_SST');
                    clear('l_not_miss')
                    l_not_miss = (any(ismember(round(hr_temp),[24:1:29])) | ...
                        any(ismember(round(hr_temp),[24:1:30]+24)) ) & ...
                        any(ismember(round(hr_temp),[24:1:29]+6))  & ...
                        any(ismember(round(hr_temp),[24:1:29]+12)) & ...
                        any(ismember(round(hr_temp),[24:1:29]+18));
                elseif strcmp(P.relative,'mean_SST');
                    clear('l_not_miss')
                    l_not_miss = any(ismember(round(hr_temp),[24:1:29])) & ...
                        any(ismember(round(hr_temp),[24:1:29]+6))  & ...
                        any(ismember(round(hr_temp),[24:1:29]+12)) & ...
                        any(ismember(round(hr_temp),[24:1:29]+18));
                end
                % It should have at least one night of measurements
 
                % *********************************************************
                % If so, then compute fundemental SSTs from two nights
                % and then compute a diurnal signal
                % *********************************************************                     
                if l_not_miss == 1, 
                    
                    if strcmp(P.relative,'fund_SST');
                        
                        clear('logic_night_1','logic_night_2','l_night')
                        logic_night_1 = round(hr_temp) >= 24 & round(hr_temp) <= 30;
                        logic_night_2 = round(hr_temp) >= 48 & round(hr_temp) <= 54;
                        l_night = logic_night_1 | logic_night_2;
                        
                        % Commented because wind and cloud are not used
                        % clear('logic_noon','logic_day')
                        % logic_noon    = hr_temp >=34 & hr_temp <= 39;
                        % logic_day     = hr_temp >=30 & hr_temp <= 42;

                        if nnz(logic_night_1) && nnz(logic_night_2),
                            qc_flag = 1;  % two nights
                            clear('fitted_sst_hr','diurnal_signal_hr','fundemental_sst_hr')
                            b = regress(C0_SST_hr(l_night)'-C0_OI_CLIM_hr(l_night)',...
                                [hr_temp(l_night)' ones(nnz(l_night),1)]);
                            fundemental_sst_hr = b(1) * hr_temp + b(2);
                        else
                            qc_flag = 0;  % one night
                            fundemental_sst_hr = repmat(nanmean(C0_SST_hr(l_night))...
                                -nanmean(C0_OI_CLIM_hr(l_night)),1,numel(l_night));
                        end
                        diurnal_signal_hr  = C0_SST_hr - C0_OI_CLIM_hr - fundemental_sst_hr;
                        clear('logic_night_1','logic_night_2','l_night')
                        
                    elseif strcmp(P.relative,'mean_SST');
                        
                        qc_flag = 1;   % pass the quality control
                        fundemental_sst_hr = ones(size(C0_SST_hr)) .* nanmean(C0_SST_hr - C0_OI_CLIM_hr);
                        diurnal_signal_hr  = C0_SST_hr - C0_OI_CLIM_hr - fundemental_sst_hr;
                        
                    end
                    % *************************************************
                    % Do the analyses: remove buoy diurnal cycles
                    % *************************************************
                    clear('lat_m','lon_m','x','y','m')
                    lon_m = LME_function_mean_period(C0_LON_hr(:),360);
                    lat_m = nanmean(C0_LAT_hr);
                    y     = discretize(lat_m,-90:5:90);
                    x = discretize(lon_m,0:5:360);
                    m     = mode(C0_MO_hr);

                    % -------------------------------------------------
                    % Find the expected diurnal signal
                    % Method 1: from buoy estimated by Chan
                    % -------------------------------------------------
                    clear('shape_exp','amplitude_exp','diurnal_exp_1')
                    shape_exp = Diurnal_Shape.Diurnal_Shape(:,y,m);
                    amplitude_exp = Diurnal_Amplitude.Diurnal_clim_buoy_1990_2014(x,y,m);
                    diurnal_exp_1 = shape_exp(C0_LCL_hr)' * amplitude_exp;
                    clear('shape_exp','amplitude_exp')

                    if all(~isnan(diurnal_exp_1)),
                        % -------------------------------------------------
                        % Find the expected diurnal signal
                        % Method 2: compute from MB16
                        % -------------------------------------------------
                        clear('m_in','lat_in','hr_in','u_in',...
                            'diurnal_exp_2','C_in','do_trd')
                        m_in   = C0_MO_hr;
                        lat_in = C0_LAT_hr;
                        hr_in  = C0_LCL_hr;
                        u_in   = C0_W_hr;
                        u_in(isnan(u_in)) = nanmean(u_in(~isnan(u_in)));
                        u_in(isnan(u_in)) = NOCS_wind_hr(isnan(u_in));
                        u_in(isnan(u_in)) = 5;
                        C_in   = C0_N_hr;
                        C_in(isnan(C_in)) = nanmean(C_in(~isnan(C_in)));
                        C_in(isnan(C_in)) = NOCS_cloud_hr(isnan(C_in));
                        C_in(isnan(C_in)) = 0.5;
                        do_trd = 1 - qc_flag;  % only compute trend if it is one night
                        diurnal_exp_2 = DIURNAL_function_compute_morak_DA...
                            (m_in, lat_in, hr_in, u_in, C_in, do_trd);
                        clear('m_in','lat_in','hr_in','u_in','C_in','do_trd')
                        
                        % -------------------------------------------------
                        % Compute the excessive diurnal signal
                        % and fit a sinusoidal function as Carella et al., 2018
                        % -------------------------------------------------
                        clear('diurnal_ex_1','diurnal_ex_2','base_x')
                        diurnal_ex_1 = diurnal_signal_hr - diurnal_exp_1;
                        diurnal_ex_2 = diurnal_signal_hr - diurnal_exp_2;
                        omega = 2*pi/24;
                        base_x = [ones(numel(C0_LCL_hr),1) sin(C0_LCL_hr'*omega) cos(C0_LCL_hr'*omega)];
                        base_x2 = [ones(numel(C0_LCL_hr),1) sin(C0_LCL_hr'*omega) cos(C0_LCL_hr'*omega),...
                            sin(C0_LCL_hr'*omega*2) cos(C0_LCL_hr'*omega*2)];
                        
                        clear('fit_1','fit_1_int','fit_2','fit_2_int','omega')
                        [fit_1,fit_1_int] = regress(diurnal_ex_1(:),base_x);
                        [fit_2,fit_2_int] = regress(diurnal_ex_2(:),base_x);

                        clear('fit_3','fit_3_int','fit_4','fit_4_int','omega')
                        [fit_3,fit_3_int] = regress(diurnal_signal_hr(:),base_x);
                        [fit_4,fit_4_int] = regress(diurnal_signal_hr(:),base_x2);
                        clear('base_x','base_x2','omega')
                        
                        % -------------------------------------------------
                        % Put data into a vector
                        % -------------------------------------------------
                        clear('uni_cty','L_cty','temp_data')
                        [uni_cty,~,L_cty] = unique(C0_CTY_CRT_hr,'rows');
                        DATA_save.m            = [DATA_save.m             m];
                        DATA_save.lon_m        = [DATA_save.lon_m         lon_m];
                        DATA_save.lat_m        = [DATA_save.lat_m         lat_m];
                        DATA_save.utc_m        = [DATA_save.utc_m         nanmean(C0_UTC_hr)];
                        DATA_save.fit_ex1      = [DATA_save.fit_ex1       fit_1];
                        DATA_save.fit_ex2      = [DATA_save.fit_ex2       fit_2];
                        DATA_save.fit_raw1     = [DATA_save.fit_raw1      fit_3];
                        DATA_save.fit_raw2     = [DATA_save.fit_raw2      fit_4];
                        DATA_save.fit_ex1_int  = cat(3,DATA_save.fit_ex1_int,fit_1_int);
                        DATA_save.fit_ex2_int  = cat(3,DATA_save.fit_ex2_int,fit_2_int);
                        DATA_save.fit_raw1_int = cat(3,DATA_save.fit_raw1_int,fit_3_int);
                        DATA_save.fit_raw2_int = cat(3,DATA_save.fit_raw2_int,fit_4_int);
                        DATA_save.deck         = [DATA_save.deck          mode(C1_DCK_hr)];
                        DATA_save.country      = [DATA_save.country;      double(uni_cty(mode(L_cty),:))];
                        DATA_save.sample_size  = [DATA_save.sample_size   numel(diurnal_ex_1)];
                        DATA_save.first_uid    = [DATA_save.first_uid     C98_UID_hr(1)];
                        DATA_save.ship_id      = [DATA_save.ship_id       ct_id];
                        DATA_save.day_id       = [DATA_save.day_id        ct_dy];
                        DATA_save.qc           = [DATA_save.qc            qc_flag];
                        DATA_save.SST_method   = [DATA_save.SST_method    mode(C0_SI_4_hr)];
                        
                        clear('uni_cty','L_cty','temp_data')
                        clear('lon_m','lat_m','m','y','x')
                        clear('fit_1','fit_1_int','fit_2','fit_2_int')
                        clear('fit_3','fit_3_int','fit_4','fit_4_int')
                        %  Format of DATA:
                        %  1. Month  2. Lon  3. Lat  4. UTC
                        %  Method 1: Buoy diurnal from Chan
                        %  5-7. excessive diurnal fitting (mean - sine - cosine)
                        %  8-9, 10-11, 12-13.    95% coverage for 5-7
                        %  Method 2: Buoy diurnal from MB16
                        %  14-16. excessive diurnal fitting (mean - sine - cosine)
                        %  17-18, 19-20, 21-22.  95% coverage for 14-16
                        %  23. Deck  24-25. Country
                        %  26. # sample in fitting  27. first UID in each day
                        %  28. Ship ID  29. Day ID  30. QC of diurnal signal
                        %  31. Measurement method
                        
                        % *************************************************
                        % Prepare for individual measurements to be saved
                        % *************************************************
                        clear('l')
                        l = round(hr_temp) >= 24 & round(hr_temp) < 48;
                        for var = 1:numel(var_list)-2
                            eval([var_list{var},'_save = [',var_list{var},'_save  ',var_list{var},'_hr];']);
                        end
                        
                        C0_CTY_CRT_save = [C0_CTY_CRT_save; C0_CTY_CRT_hr];
                        FD_save     = [FD_save     fundemental_sst_hr];
                        DA_save     = [DA_save     diurnal_signal_hr];
                        D1_EXP_save = [D1_EXP_save diurnal_exp_1];
                        D2_EXP_save = [D2_EXP_save diurnal_exp_2];
                        D1_EX_save  = [D1_EX_save  diurnal_ex_1];
                        D2_EX_save  = [D2_EX_save  diurnal_ex_2];
                        Ship_save   = [Ship_save ones(1,numel(l))*ct_id];
                        DAY_save    = [DAY_save  ones(1,numel(l))*ct_dy];
                        
                        if strcmp(P.relative,'fund_SST');
                            
                            if qc_flag == 1,
                                qc_temp = double(l)*qc_flag;
                            else
                                qc_temp = double(l)*(1-qc_flag) + 2;
                            end
                            
                        elseif strcmp(P.relative,'mean_SST');
                            
                            qc_temp = double(l)*qc_flag;
                            
                        end
                        DAY_indicator_save = [DAY_indicator_save  qc_temp];
                        clear('l','var','qc_flag','qc_temp')
                    end
                end
                clear('l_not_miss')
            end
        end
        clear('spd')
    end

    % *********************************************************************
    % Clear data
    % *********************************************************************
    for var = 1:numel(var_list)
        eval(['clear(''',var_list{var},''')'])
        eval(['clear(''',var_list{var},'_id'')'])
        eval(['clear(''',var_list{var},'_dy'')'])
        eval(['clear(''',var_list{var},'_hr'')'])
    end

    clear('C0_CTY_CRT_id','NOCS_wind_id','NOCS_cloud_id')
    clear('NOCS_wind_dy','NOCS_cloud_dy','C0_CTY_CRT_dy')
    clear('NOCS_wind_hr','NOCS_cloud_hr','C0_CTY_CRT_hr')
    
    clear('l_empty','l_excld','l_gener','l_night',...
        'l_not_miss','logic_night_1','logic_night_2',...
        'a','b','J_ID','L_dup','ct_id','ct_dy','day_list',...
        'day_uni','diurnal_signal_hr','fundemental_sst_hr',...
        'i','l','lcl','mon','spd','spd_shtreshold','temp',...
        'Diurnal_Amplitude','Diurnal_Shape','NOCS_cloud','NOCS_wind')

    % *********************************************************************
    % Save Data
    % *********************************************************************
    disp('Saving Data ...')
    if ~isempty(C98_UID_save),
        
        for mon = 1:12
            
            clear('l_ind','l_day')
            l_ind = C0_MO_save(1,:) == mon;
            l_day = DATA_save.m     == mon;
            
            file_save = [dir_save,'IMMA1_R3.0.0_',num2str(yr),...
                                  '-',CDF_num2str(mon,2),'_Ship_Diurnal_Signal',...
                                  '_relative_to_',P.relative,'.mat'];            

            if nnz(l_ind) ~=0,
                
                for var = 1:numel(var_list)-2
                    eval(['clear(''',var_list{var},''')'])
                    eval([var_list{var},' = ',var_list{var},'_save(l_ind);']);
                end
                
                clear('C0_CTY_CRT','Fundemental_SST','Diurnal_signal')
                clear('D1_EXP','D2_EXP','D1_EX','D2_EX')
                clear('Ship_indices','Day_indices','Day_indicator')
                C0_CTY_CRT      = C0_CTY_CRT_save(l_ind,:);
                Fundemental_SST = FD_save(l_ind);
                Diurnal_signal  = DA_save(l_ind);
                D1_EXP          = D1_EXP_save(l_ind);
                D2_EXP          = D2_EXP_save(l_ind);
                D1_EX           = D1_EX_save(l_ind);
                D2_EX           = D2_EX_save(l_ind);
                Ship_indices    = Ship_save(l_ind);
                Day_indices     = DAY_save(l_ind);
                Day_indicator   = DAY_indicator_save(l_ind);
                % DAY_indicator: whether measurements belong to the second day (0-6hr).
                
                clear('DATA')
                % DATA = DATA_save(:,l_day);
                DATA.m            = DATA_save.m(l_day);
                DATA.lon_m        = DATA_save.lon_m(l_day);
                DATA.lat_m        = DATA_save.lat_m(l_day);
                DATA.utc_m        = DATA_save.utc_m(l_day);
                DATA.fit_ex1      = DATA_save.fit_ex1(:,l_day);
                DATA.fit_ex2      = DATA_save.fit_ex2(:,l_day);
                DATA.fit_raw1     = DATA_save.fit_raw1(:,l_day);
                DATA.fit_raw2     = DATA_save.fit_raw2(:,l_day);
                DATA.fit_ex1_int  = DATA_save.fit_ex1_int(:,:,l_day);
                DATA.fit_ex2_int  = DATA_save.fit_ex2_int(:,:,l_day);
                DATA.fit_raw1_int = DATA_save.fit_raw1_int(:,:,l_day);
                DATA.fit_raw2_int = DATA_save.fit_raw2_int(:,:,l_day);
                DATA.deck         = DATA_save.deck(l_day);   
                DATA.country      = DATA_save.country(l_day,:);
                DATA.sample_size  = DATA_save.sample_size(l_day);
                DATA.first_uid    = DATA_save.first_uid(l_day);
                DATA.ship_id      = DATA_save.ship_id(l_day);
                DATA.day_id       = DATA_save.day_id(l_day);
                DATA.qc           = DATA_save.qc(l_day);
                DATA.SST_method   = DATA_save.SST_method(l_day);
                
                % format_of_DATA = [...
                %     '1. Month  2. Lon  3. Lat  4. UTC  ',char(10),...
                %     'Method 1: Buoy diurnal from Chan  ',char(10),...
                %     '5-7. excessive diurnal fitting (mean - sine - cosine)  ',char(10),...
                %     '8-9, 10-11, 12-13.    95% coverage for 5-7  ',char(10),...
                %     'Method 2: Buoy diurnal from MB16  ',char(10),...
                %     '14-16. excessive diurnal fitting (mean - sine - cosine)  ',char(10),...
                %     '17-18, 19-20, 21-22.  95% coverage for 14-16  ',char(10),...
                %     '23. Deck  24-25. Country  ',char(10),...
                %     '26. # sample in fitting  27. first UID in each day  ',char(10),...
                %     '28. Ship ID  29. Day ID  30. QC of diurnal signal  ',char(10),...
                %     '31. Measurement method'];

                save(file_save,'C0_YR','C0_MO','C0_DY','C0_HR','C0_LCL','C0_UTC',...
                     'C0_LON','C0_LAT','C0_SST','C0_SI_4','C0_AT','C0_II',...
                     'C1_DCK','C98_UID','C0_OI_CLIM','ID_uni','C0_CTY_CRT',...
                     'Fundemental_SST','Diurnal_signal',...
                     'D1_EXP','D2_EXP','D1_EX','D2_EX',...
                     'Ship_indices','Day_indices','Day_indicator',...
                     'DATA','-v7.3');
            end
        end
    end
end