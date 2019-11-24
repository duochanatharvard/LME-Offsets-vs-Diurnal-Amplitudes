function output = DA_LME_function_statistics(yr_start,reg_sea)

    data  = load('All_lme_offsets_and_diurnal_amplitudes.mat');
    grp   = data.grp;
    yr0   = 1879;                                  % The first year of data

    x         = data.da(2:end,yr_start-yr0,reg_sea);
    y         = data.lme(2:end,yr_start-yr0,reg_sea);
    x_std     = data.da_std(2:end,yr_start-yr0,reg_sea);
    y_std     = data.lme_std(2:end,yr_start-yr0,reg_sea);
    x_eri     = data.da(1,yr_start-yr0,reg_sea);
    y_eri     = data.lme(1,yr_start-yr0,reg_sea);
    x_eri_std = data.da_std(1,yr_start-yr0,reg_sea);
    y_eri_std = data.lme_std(1,yr_start-yr0,reg_sea);


    % *********************************************************************
    % Panel a. correlation and r^2
    % *********************************************************************
    [R, N_eff, R2, R2_quan] = CDC_corr(x,y,1);
    output.R      = R;
    output.N      = N_eff;
    output.R2     = R2;
    output.R2_quan = R2_quan;
    
    % *********************************************************************
    % Panel b. york fit
    % *********************************************************************
    N = 100;
    P.mute_output = 1;
    [slp, ipt, slp_member, ipt_member] = CDC_yorkfit_bt(y,x,y_std,x_std,0,1,N,P); 
    slp_quan = quantile(slp_member,[0.025 0.25 0.75 0.975],1);
    ipt_quan = quantile(ipt_member,[0.025 0.25 0.75 0.975],1);
    output.slp = slp;
    output.slp_quan = slp_quan;
    
    fitted_eri = slp_member*x_eri + ipt_member;
    output.is_ERI_in = y_eri > min(fitted_eri) & y_eri < max(fitted_eri); 
    
    % *********************************************************************
    % Panel c. quantile for diurnal cycle and LME
    % *********************************************************************
    output.da_quan_ex_ERI  = quantile(x,[0 0.025 0.25 0.5 0.75 0.975 1]);
    output.da_quan         = quantile([x;x_eri],[0 0.025 0.25 0.5 0.75 0.975 1]);
    output.lme_quan_ex_ERI = quantile(y,[0 0.025 0.25 0.5 0.75 0.975 1]);
    output.lme_quan        = quantile([y;y_eri],[0 0.025 0.25 0.5 0.75 0.975 1]);
    
    % *********************************************************************
    % Porject and find the marginal distrubution of individual points
    %                                  and find end point of bucket and ERI
    % *********************************************************************
    [a, b, a_mem, b_mem] =  CDC_yorkfit_bt([y;y_eri],[x;x_eri],[y_std;0],[x_std;0],0,1,N,P); 
    n_mu     = (a.*(y - b).*x_std.^2 + x.*y_std.^2) ./ (y_std.^2 + x_std.^2.*a.^2);
    end_eri  = (a.*(y_eri - b).*x_eri_std.^2 + x_eri.*y_eri_std.^2) ./ (y_eri_std.^2 + x_eri_std.^2.*a.^2);
    end_bck = max(n_mu);
    output.pctg_mu  = (n_mu - end_eri) ./ (end_bck - end_eri);

    % *********************************************************************
    % Compute uncertainties for percentage estimates
    % Account for errors in both boostrapped linear fitting and
    %                                        projected normal distributions
    % *********************************************************************
    clear('p_mem')
    for ct = 1:100
        clear('n_std_mem','n_mu_mem','p_mu_mem','p_std_mem','end_bck_mem','end_eri_mem')
        n_mu_mem    = (a_mem(ct).*(y - b_mem(ct)).*x_std.^2 + x.*y_std.^2) ./ (y_std.^2 + x_std.^2.*a_mem(ct).^2); 
        end_eri_mem = (a_mem(ct).*(y_eri - b_mem(ct)).*x_eri_std.^2 + x_eri.*y_eri_std.^2) ./ (y_eri_std.^2 + x_eri_std.^2.*a_mem(ct).^2); 
        n_std_mem   = sqrt((x_std.^2 .* y_std.^2) ./ (y_std.^2 + x_std.^2.*a_mem(ct).^2));
        end_bck_mem = max(n_mu_mem);
        p_mu_mem    = (n_mu_mem - end_eri_mem) ./ (end_bck - end_eri_mem);
        p_std_mem   = n_std_mem ./ (end_bck_mem - end_eri_mem);
        p_mem(:,ct) = normrnd(p_mu_mem,p_std_mem,numel(p_mu_mem),1);
    end
    output.pctg_std = CDC_std(p_mem,2);
    
end