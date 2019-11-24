function [da_out,bias_out] = DA_LME_function_bucket_model(ct_reg,t_id,shade_id,...
                      alpha_id,size_id,mixing_id,wind_id,thick_id,eri_bias)

    %% *********************************************************************
    % Generating the case name and run models
    % *********************************************************************
    O = DA_LME_function_get_region_name(ct_reg);
    save_app = ['_',O.region_name_m];
    col = lines(7);  col(7,:) = [0 .6 .6];
    
    if numel(t_id) == 1,
        save_app = [save_app,'_int_time_',num2str((t_id-1)/2)];
    else
        save_app = [save_app,'_int_time_',num2str((t_id(1)-1)/2),'-',num2str((t_id(end)-1)/2)];
        lgd = {'1min','3min','5min','7min','9min','11min    '};
        col_model = col(1,:);
        ct_fig = 1;
        num_var = numel(t_id);
    end

    if numel(shade_id) == 1,
        save_app = [save_app,'_shading_',num2str(shade_id/10-0.1,'%6.1f')];
    else
        save_app = [save_app,'_shading_',num2str(shade_id(1)/10-0.1,'%6.1f'),'-',num2str(shade_id(end)/10-0.1,'%6.1f')];
        lgd = {'0%','20%','40%','60%','80%','100%    '};
        col_model = col(2,:);
        ct_fig = 2;
        num_var = numel(shade_id);
    end

    if numel(alpha_id) == 1,
        save_app = [save_app,'_air_signal_',num2str(alpha_id/20-0.05,'%6.2f')];
    else
        save_app = [save_app,'_air_signal_',num2str(alpha_id(1)/20-0.05,'%6.2f'),'-',num2str(alpha_id(end)/20-0.05,'%6.2f')];
        lgd = {'0%','5%','10%','15%','20%    '};
        col_model = col(3,:);
        ct_fig = 3;
        num_var = numel(alpha_id);
    end

    if numel(size_id) == 1,
        save_app = [save_app,'_size_',num2str(size_id)];
    else
        save_app = [save_app,'_size_',num2str(size_id(1)),'-',num2str(size_id(end))];
        lgd = {'large','medium','small    '};
        col_model = col(4,:);
        ct_fig = 4;
        num_var = numel(size_id);
    end
    
    if numel(mixing_id) == 1,
        save_app = [save_app,'_bck_ratio_',num2str(mixing_id/10 -0.1,'%6.1f')];
    else
        save_app = [save_app,'_bck_ratio_',num2str(mixing_id(1)/10 -0.1,'%6.1f'),'-',num2str(mixing_id(end)/10 -0.1,'%6.1f')];
        lgd = {'20%','40%','60%','80%','100%    '};
        col_model = col(5,:);
        ct_fig = 5;
        num_var = numel(mixing_id);
    end

    if numel(wind_id) == 1,
        save_app = [save_app,'_wind_',num2str(wind_id)];
    else
        save_app = [save_app,'_wind_',num2str(wind_id(1)),'-',num2str(wind_id(end))];
        lgd = {'x0','x0.2','x0.4','x0.6','x0.8','x1','x1.2    '};
        col_model = col(6,:);
        ct_fig = 6;
        num_var = numel(wind_id);
    end

    if numel(thick_id) == 1,
        save_app = [save_app,'_thickness_',num2str(thick_id)];
    else
        save_app = [save_app,'_thickness_',num2str(thick_id(1)),'-',num2str(thick_id(end))];
        lgd = {'0.2cm','0.35cm','0.5cm','0.75cm','1cm','1.5cm','2cm    '};
        lgd = fliplr(lgd);
        col_model = col(7,:);
        ct_fig = 7;
        num_var = numel(thick_id);
    end    
     
    save_app = [save_app,'_eri_bias_',num2str(eri_bias,'%6.1f')];

    % ---------------------------------------------------------------------
    % Model Simulations
    % ---------------------------------------------------------------------
    % Change to the directory of bucket models                       % TODO
    P_input.average_forcing = 1;
    [true_SST,true_AT,e_air,u_environment,Qs,direct_ratio,zenith_angle] = ...
                                                  BKT_MD_PREP(P_input);
    clear('P')
    % Set variables with default values                              % TODO
    P.deck_time       = 240;
    P.s_environment   = 10;
    P.solar_shading   = 0.4;
    P.thickness       = 0.01; 
    P.wind_experience = 0.6;
    P.diamter         = 0.25;   
    P.depth           = 0.2;
    alpha             = (1-1) * 0.05;
    init_SST          = true_AT * alpha + true_SST * (1-alpha);
    thick_list        = [0.2 0.35 0.5 0.75 1 1.5 2]*0.01;

    PP.do_sensible    = 1; 
    PP.do_latent      = 1;
    PP.do_long        = 1;
    PP.do_solar       = 1;

    switch ct_reg,
        case 1,
            ct = 1;     % 20S-20N
        case 3,
            ct = 2;     % 20N-40N
        case 7,
            ct = 3;     % 40N-60N
    end
        
    if ct_fig ~=5,
        for ct_var = 1:num_var

            switch ct_fig,
                case 1,
                    P.deck_time     = (t_id(ct_var)-3)*30;
                case 2,
                    P.solar_shading = (shade_id(ct_var)-1)*0.1;
                case 3,
                    alpha           = (alpha_id(ct_var)-1)*0.05;
                    init_SST        = true_AT * alpha + true_SST * (1-alpha);
                case 4,
                    switch size_id(ct_var),
                        case 1,
                             P.diamter = 0.25;   P.depth = 0.2;
                        case 2,
                             P.diamter = 0.163;  P.depth = 0.14;
                        case 3,
                             P.diamter = 0.08;   P.depth = 0.12;
                    end
                case 6,
                    P.wind_experience = (wind_id(ct_var) - 1) * 0.2;
                case 7,
                    P.thickness       = thick_list(thick_id(ct_var));
            end

            SST_raw  = BKT_MD_STP_2_MD_WOODEN_GRD_SIZ_for_Chan2020(...
                init_SST(ct,1,:,:),true_AT(ct,1,:,:),e_air(ct,1,:,:),...
                u_environment(ct,1,:,:),Qs(ct,1,:,:),direct_ratio(ct,1,:,:),...
                zenith_angle(ct,1,:,:),P,PP);

            t_w(:,:,ct_var) = SST_raw(1,1,:,:,end);
            t_i(:,:,ct_var) = SST_raw(1,1,:,:,1);
        end
        t_t   = squeeze(true_SST(ct,:,:,:));
    else
        SST_raw  = BKT_MD_STP_2_MD_WOODEN_GRD_SIZ_for_Chan2020(...
            init_SST(ct,1,:,:),true_AT(ct,1,:,:),e_air(ct,1,:,:),...
            u_environment(ct,1,:,:),Qs(ct,1,:,:),direct_ratio(ct,1,:,:),...
            zenith_angle(ct,1,:,:),P,PP);

        t_w(:,:,1) = SST_raw(1,1,:,:,end);
        t_i(:,:,1) = SST_raw(1,1,:,:,1);
        t_t        = squeeze(true_SST(ct,:,:,:));        
    end

    % ---------------------------------------------------------------------
    % fit for the diurnal amplitude
    % ---------------------------------------------------------------------
    C0_LCL_hr = [1:24];
    omega = 2*pi/24;
    base_x = [ones(numel(C0_LCL_hr),1) sin(C0_LCL_hr'*omega) cos(C0_LCL_hr'*omega)];
    base_x2 = [ones(numel(C0_LCL_hr),1) sin(C0_LCL_hr'*omega) cos(C0_LCL_hr'*omega)  sin(C0_LCL_hr'*omega*2) cos(C0_LCL_hr'*omega*2)];

    clear('fit_1','fit_1_int','fit_2','fit_2_int','omega')
    for ct_x = 1:size(t_w,2)
        for ct_y = 1:size(t_w,3)
            f_w(1:3,ct_x,ct_y) = lscov(base_x,t_w(:,ct_x,ct_y));
            f_i(1:3,ct_x,ct_y) = lscov(base_x,t_i(:,ct_x,ct_y));
            f_w(5:9,ct_x,ct_y) = lscov(base_x2,t_w(:,ct_x,ct_y));
            f_i(5:9,ct_x,ct_y) = lscov(base_x2,t_i(:,ct_x,ct_y));
        end
        f_t(1:3,ct_x) = lscov(base_x,t_t(:,ct_x));
        f_t(5:9,ct_x) = lscov(base_x2,t_t(:,ct_x));
    end 
    
    % Load ERI diurnal cycles ---------------------------------------------
    yr_start   = 1990;
    yr_end     = 2009;
    file_load = ['STATS_ERI_DA_signals_',num2str(yr_start),...
        '_',num2str(yr_end),'_',O.region_name,'_relative_to_mean_SST.mat'];
    ERI = load(file_load,'fit_out_adj');
    
    % reconstruct ERI diurnal cycle from fitting --------------------------
    C0_LCL_hr = [1:24];
    omega = 2*pi/24;
    base_x2 = [ones(numel(C0_LCL_hr),1) sin(C0_LCL_hr'*omega) cos(C0_LCL_hr'*omega)  sin(C0_LCL_hr'*omega*2) cos(C0_LCL_hr'*omega*2)];
    for ct = 1:size(ERI.fit_out_adj,3)
        for ct_sea = 1:size(ERI.fit_out_adj,2)
            ERI.diurnal_adj(:,ct_sea,ct) = base_x2 * ERI.fit_out_adj(4:8,ct_sea,ct);
        end
    end
    clear('C0_LCL_hr','base_x2','ct_speed','ct_sea','omega','yr_end','yr_start')
    
    % compute ERI diurnal signals, assuming that ERI is always [eri_bias]^oC warmer
    t_eri = ERI.diurnal_adj(:,:,1) + repmat(nanmean(t_t(:,:,1),1),24,1) + eri_bias;
    f_eri = ERI.fit_out_adj([1:3 1 4:8],:,1);
    f_eri([1 4 5],:) = f_t([1 4 5],:,1) + eri_bias;
    
    clear('t_b','f_b')
    mixing = (mixing_id/10) -0.1;
    if numel(mixing_id) == 1,
        t_b = t_w .* mixing + repmat(t_eri,1,1,size(t_w,3)) .* (1-mixing);
        f_b = f_w .* mixing + repmat(f_eri,1,1,size(f_w,3)) .* (1-mixing);
    else
        for ct = 1:numel(mixing)
            t_b(:,:,ct) = t_w .* mixing(ct) + t_eri * (1-mixing(ct));
            f_b(:,:,ct) = f_w .* mixing(ct) + f_eri * (1-mixing(ct));
        end
    end
    
    if numel(shade_id) > 1 || numel(thick_id) > 1,
        t_b = t_b(:,:,[end:-1:1]);
        f_b = f_b(:,:,[end:-1:1]);
    end
    
    [p_b,a_b]  = cart2pol(f_b(7,:,:),f_b(6,:,:));
    p_b(p_b<0) = (p_b(p_b<0) + 2*pi) / pi * 12;

    [p_i,a_i]  = cart2pol(f_i(7,:,:),f_i(6,:,:));
    p_i(p_i<0) = (p_i(p_i<0) + 2*pi) / pi * 12;

    [p_eri,a_eri]  = cart2pol(f_eri(7,:,:),f_eri(6,:,:));
    p_eri(p_eri<0) = (p_eri(p_eri<0) + 2*pi) / pi * 12;
    
    [p_t,a_t]  = cart2pol(f_t(7,:,:),f_t(6,:,:));
    p_t(p_t<0) = (p_t(p_t<0) + 2*pi) / pi * 12;
    
    % compute biases of individual types of measurements
    bias     = f_b(1,:,:)   - repmat(f_t(1,:,1),1,1,size(f_b,3));
    bias_t   = f_t(1,:,:)   - f_t(1,:,:);
    bias_eri = f_eri(1,:,1) - f_t(1,:,1);
    if ct_fig ~= 5,
        bias_i   = f_i(1,:,:)   - repmat(f_t(1,:,1),1,1,size(f_b,3));
    else
        bias_i   = repmat(f_i(1,:,:) - f_t(1,:,1),1,1,size(f_b,3));
    end

    
    % *********************************************************************
    % Generating Figures...
    % *********************************************************************
    % close all;
    modif_list = linspace(0.2,1,size(f_b,3)+1);
    modif_list = [modif_list(1:2:end) modif_list(end-2:-2:1)];
    
    if ismember(ct_reg,[1 2]),
        sea_list = [1];
    else
        sea_list = [2 3];
    end

    % *********************************************************************
    % 4. Daily mean bias versus amplitude
    % *********************************************************************
    figure(ct_reg); hold on;
    col = lines(3);
    for ct_sea = sea_list
        
        switch ct_sea,
            case 1,
                l_sea = 1:12;
                col_sea = 'k';
                st = 'o';
            case 2,
                l_sea = [1 2 12];
                col_sea = 'b';
                st = 'o';
            case 3,
                l_sea = [6 7 8];
                col_sea = 'r';
                st = 'd';
        end

        mksiz = 13;
        clear('h')
        h(1) = plot(nanmean(a_t(1,l_sea,1)),nanmean(bias_t(1,l_sea,1)),'d','linewi',2,'markersize',mksiz-7,'color',col_sea);
        h(1) = plot(nanmean(a_t(1,l_sea,1)),nanmean(bias_t(1,l_sea,1)),'d','linewi',2,'markersize',mksiz+3,'color',col_sea);
        h(2) = plot(nanmean(a_eri(1,l_sea,1)),nanmean(bias_eri(1,l_sea,1)),'o','linewi',2,'markersize',mksiz-7,'color',col_sea);
        h(2) = plot(nanmean(a_eri(1,l_sea,1)),nanmean(bias_eri(1,l_sea,1)),'o','linewi',2,'markersize',mksiz+3,'color',col_sea);

        mksiz = 6;

        plot1 = plot(squeeze(nanmean(a_b(1,l_sea,:),2)),squeeze(nanmean(bias(1,l_sea,:),2)),'-','color',col_model,'linewi',2);
        plot1.Color(4) = 1;
        clear('h')

        for ct = 1:size(f_b,3)
            modif = modif_list(ct);
            if ct <= size(f_b,3)/2,
                col_use = col_model * modif;
            else
                col_use = 1 - modif + col_model * modif;
            end
            x     = nanmean(a_b(1,l_sea,ct));
            y1    = nanmean(bias(1,l_sea,ct));
            
            figure(ct_reg);
            h(1) = plot(x,y1,st,'color',col_model,'markerfacecolor',...
                col_use,'markersize',mksiz+3,'linewi',2);
            
        end
        
        da_out(:,ct_sea)   = squeeze(nanmean(a_b(1,l_sea,:)));
        bias_out(:,ct_sea) = squeeze(nanmean(bias(1,l_sea,:)));
        
        figure(ct_reg);
        CDF_panel([0 0.4 -0.5 0.2],'','','Diurnal amplitude (^oC)','Bias (^oC)','fontsize',18);
        daspect([0.4 0.7 1])

        if ct_reg == 1,
            figure(10); hold on;
            pic_list = [3 6 3 3 2 2 7];
            ct = pic_list(ct_fig);
            modif = modif_list(ct);
            if ct <= size(f_b,3)/2,
                col_use = col_model * modif;
            else
                col_use = 1 - modif + col_model * modif;
            end
            if ct_fig == 1,
                pic_t = nanmean(t_t(:,l_sea),2) - 273.15;
                plot(1:24,pic_t,'color',[1 1 1]*.7,'linewi',5)
                pic_b = nanmean(t_b(:,l_sea,ct),2) - 273.15;
                plot(1:24,pic_b,'color',col_model,'linewi',5);
            end
            pic_b = nanmean(t_b(:,l_sea,ct),2) - 273.15;
            plot(1:24,pic_b,'-','color',col_model,'linewi',2);
            if ct_fig == 1,
                CDF_panel([0 25 26.5 27.5],'','','Local hour (hr)','Water temperature (^oC)','fontsize',18);
                daspect([25 1 1])
            end
        end
        
    end

    % The following code is for legends -----------------------------------
    if ct_reg == 1,
        figure(ct_fig + 10); hold on;
        for ct_sea = 1

            l_sea = 1:12;
            col_sea = 'k';
            st = 'o';

            mksiz = 7;

            plot(squeeze(nanmean(a_b(1,l_sea,:),2)),squeeze(nanmean(bias(1,l_sea,:),2)),'-','color',col_model,'linewi',4);
            clear('h')

            for ct = 1:size(f_b,3)
                modif = modif_list(ct);
                if ct <= size(f_b,3)/2,
                    col_use = col_model * modif;
                else
                    col_use = 1 - modif + col_model * modif;
                end
                x     = nanmean(a_b(1,l_sea,ct));
                y1    = nanmean(bias(1,l_sea,ct));
                h(ct) = plot(x,y1,st,'color',col_model,'markerfacecolor',...
                    col_use,'markersize',mksiz+3,'linewi',2);
            end

            legend(h,[lgd],'fontsize',18,'location','eastoutside');
        end
    end
end
