clear;

% *************************************************************************
% Figure for offseet-diurnal relationship
% in response to varying distinct bucket model parameters.
%
% Reference bucket
% t_id      = 11;     % 1. initial SST        +2 for every 1 minute
% shade_id  = 4;      % 1. all insolation     11. no insolation
% alpha_id  = 1;      % 1. no air temperature  5. 20% air temperature
% size_id   = 1;      % 1. large bucket        3. small bucket
% mixing_id = 11;     % 1. all eri            11. all bucket
% wind_id   = 3;      % 1. less wind (x0)      7. more wind (x3)
% thick_id  = 5;      % 1. thing bucket (2mm)  7. thick bucket (2cm)
%
%  DA_LME_function_bucket_model(ct_reg,t_id,shade_id,alpha_id,size_id,
%                                      mixing_id,wind_id,thick_id,eri_bias)
%
% This is a faster version that uses regional and seasonal mean fields to
% drive the extended bucket models.  We have compaired the results with
% running model at 5^o monthly resolution first and than computing averages.
% Errors are within 2% for the tropics and witin 10\% outside the tropics.
% This yeilds a faster implementation yet does not qulitatively change the
% main conlusion of the paper, that is:
% A negative slope is expected between diurnal amplitudes and LME offsets
% *************************************************************************
for ct_reg = [1 3 7]      % 1: tropics    3: subtropics    7: extra-tropics
    [da_out{1,ct_reg},bias_out{1,ct_reg}] = ...
        DA_LME_function_bucket_model(ct_reg,3:4:23,     5,     1,   1,      11,   3,   5, 0.1);
    [da_out{2,ct_reg},bias_out{2,ct_reg}] = ...
        DA_LME_function_bucket_model(ct_reg,    11,1:2:11,     1,   1,      11,   3,   5, 0.1);
    [da_out{3,ct_reg},bias_out{3,ct_reg}] = ...
        DA_LME_function_bucket_model(ct_reg,    11,     5, 1:1:5,   1,      11,   3,   5, 0.1);
    [da_out{4,ct_reg},bias_out{4,ct_reg}] = ...
        DA_LME_function_bucket_model(ct_reg,    11,     5,     1, 1:3,      11,   3,   5, 0.1);
    [da_out{5,ct_reg},bias_out{5,ct_reg}] = ...
        DA_LME_function_bucket_model(ct_reg,    11,     5,     1,   1,  3:2:11,   3,   5, 0.1);
    [da_out{6,ct_reg},bias_out{6,ct_reg}] = ...
        DA_LME_function_bucket_model(ct_reg,    11,     5,     1,   1,      11, 1:7,   5, 0.1);
    [da_out{7,ct_reg},bias_out{7,ct_reg}] = ...
        DA_LME_function_bucket_model(ct_reg,    11,     5,     1,   1,      11,   3, 1:7, 0.1);
end
save('Model_output.mat','da_out','bias_out','-v7.3');

% Compute slopes for model outputs
for ct_reg = [1:5]
    for ct_para = 1:7
        switch ct_reg
            case 1
                temp = CDC_trend(bias_out{ct_para,1}(:,1),da_out{ct_para,1}(:,1),1);
            case 2
                temp = CDC_trend(bias_out{ct_para,3}(:,2),da_out{ct_para,3}(:,2),1);
            case 3
                temp = CDC_trend(bias_out{ct_para,3}(:,3),da_out{ct_para,3}(:,3),1);
            case 4
                temp = CDC_trend(bias_out{ct_para,7}(:,2),da_out{ct_para,7}(:,2),1);
            case 5
                temp = CDC_trend(bias_out{ct_para,7}(:,3),da_out{ct_para,7}(:,3),1);
        end
        Tab(ct_para,ct_reg) = temp{1};
    end
end

for ct = 11:17  close(ct);  end  % comment this line if legend is required

%% *************************************************************************
% Figure 3: Tropics
% *************************************************************************
% for ct_yr_start = 1890:20:1990
clear;
revision = 0;  % sensitivity tests:
% 0: standard analysis
% 1: bucket identification priortizing WMO No.47 metadata
% 2: do not infer bucket SSTs from country information
% 3: nation-level analysis
% 4: excess diurnal cycles

for ct_yr_start = 1890:20:1990
    for ct_reg_sea = 1
        DA_LME_function_scatter_plot('bucket',ct_yr_start,ct_reg_sea,revision);
    end
end

%% *************************************************************************
% Figure 4: NH Subtropics and Extra-tropics
% *************************************************************************
clear;
close all;
revision = 0;   % sensitivity tests:
% 0: standard analysis
% 1: bucket identification priortizing WMO No.47 metadata
% 2: do not infer bucket SSTs from country information
% 3: nation-level analysis
% 4: excess diurnal cycles

for ct_yr_start = 1970
    for ct_reg_sea = 2:5
        DA_LME_function_scatter_plot('bucket',ct_yr_start,ct_reg_sea,revision);
    end
end

%% *************************************************************************
% Data for Figure 5: First compute statistics and then generate plot
% *************************************************************************
clear;
close all;
revision = 0;   % sensitivity tests:
% 0: standard analysis
% 1: bucket identification priortizing WMO No.47 metadata
% 2: do not infer bucket SSTs from country information
% 3: nation-level analysis
% 4: excess diurnal cycles

for ct_yr_start = 1880:1:1990

    disp(['Year ',num2str(ct_yr_start)])
    for ct_reg_sea = 1:5

        output = DA_LME_function_statistics(ct_yr_start,ct_reg_sea,revision);
        ct_prd = ct_yr_start - 1879;

        Tab_r(ct_prd,ct_reg_sea)                   = output.R;
        Tab_n(ct_prd,ct_reg_sea)                   = output.N;
        Tab_r2(ct_prd,ct_reg_sea)                  = output.R2;
        Tab_r2_quan(ct_prd,ct_reg_sea,:)           = output.R2_quan;

        Tab_slp(ct_prd,ct_reg_sea)                 = output.slp;
        Tab_ipt(ct_prd,ct_reg_sea)                 = output.ipt;
        Tab_slp_quan(ct_prd,ct_reg_sea,:)          = output.slp_quan;
        Tab_is_eri_in(ct_prd,ct_reg_sea)           = output.is_ERI_in;

        Tab_da_quan(ct_prd,ct_reg_sea,:)           = output.da_quan;
        Tab_da_quan_ex_ERI(ct_prd,ct_reg_sea,:)    = output.da_quan_ex_ERI;
        Tab_lme_quan(ct_prd,ct_reg_sea,:)          = output.lme_quan;
        Tab_lme_quan_ex_ERI(ct_prd,ct_reg_sea,:)   = output.lme_quan_ex_ERI;

    end
end

switch revision
    case 0
        app = '';
    case 1
        app = '_sens_wmo';
    case 2
        app = '_sens_not_infer';
    case 3
        app = '_sens_nation_level';
    case 4
        app = '_sens_excess_diurnal';
end
save(['DA_LME_Statistics',app,'.mat'],...
    'Tab_r','Tab_n','Tab_r2','Tab_r2_quan','Tab_slp','Tab_slp_quan',...
    'Tab_da_quan','Tab_da_quan_ex_ERI','Tab_lme_quan','Tab_lme_quan_ex_ERI',...
    'Tab_is_eri_in','Tab_ipt','-v7.3')

%% *************************************************************************
% Figure 5: First compute statistics and then generate plot
% *************************************************************************
clear;
close all;
revision = 0;   % sensitivity tests:
% 0: standard analysis
% 1: bucket identification priortizing WMO No.47 metadata
% 2: do not infer bucket SSTs from country information
% 3: nation-level analysis
% 4: excess diurnal cycles

switch revision
    case 0
        app = '';
    case 1
        app = '_sens_wmo';
    case 2
        app = '_sens_not_infer';
    case 3
        app = '_sens_nation_level';
    case 4
        app = '_sens_excess_diurnal';
end

load(['DA_LME_Statistics',app,'.mat'])
if revision <= 3
    load(['All_lme_offsets_and_diurnal_amplitudes',app,'.mat'])
else
    load(['All_lme_offsets_and_diurnal_amplitudes.mat']);
end

close all;
c = 0;
for ct_fig = [1 2 3 4 3]
    % 1. R^2   2. Diurnal amplitude range   3. LME offset range   4. Slopes

    c = c + 1;
    
    sea_list = [1];                   % Figure 5
    % sea_list = [2 3 4 5];             % Figure A3
    % 1. Tropical annual mean
    % 2. DJF   3. JJA  for 20-40N
    % 4. DJF   5. JJA  for 40-60N

    switch ct_fig
        case 1
            pic     = Tab_r2;
            pic_std = Tab_r2_quan;
            y_label = 'R^2';
            offset  = -1;
            aa = [0.2 0.5 0.8];
            st = '-';
        case 2
            pic     = Tab_da_quan_ex_ERI(:,:,4);
            pic_std = Tab_da_quan_ex_ERI(:,:,[2 3 5 6]);
            y_label = 'Diurnal Amp. (^oC)';
            offset  = - 0.5;
            aa = [0.1 0.4];
            st = '-';
        case 3
            if c == 5
                pic     = Tab_lme_quan(:,:,end) - Tab_lme_quan(:,:,5);
                st = '--';
            else
                pic     = Tab_lme_quan_ex_ERI(:,:,end) - Tab_lme_quan_ex_ERI(:,:,5);
                st = '-';
            end
            pic_std = Tab_lme_quan_ex_ERI(:,:,[2 3 5 6]) * nan;
            y_label = 'LME offset range (^oC)';
            offset  = -1;
            aa = [0.3 .7];
            
        case 4
            pic     = Tab_slp;
            pic_std = Tab_slp_quan;
            y_label = 'Slope (^oC / ^oC)';
            if revision ~= 2
                offset  = -20;
            else
                offset = -60;
            end
            aa = [-5 0 5];
            st = '-';
    end

    if ct_fig == 4
        y_lim   = [offset*(numel(sea_list)-0.5) offset*-.5];
    else
        y_lim   = [offset*(numel(sea_list)-1) offset*-1];
    end
    y_tick = []; for i = numel(sea_list):-1:1 y_tick  = [y_tick aa+offset*(i-1)];  end
    y_tick_label = repmat(aa,1,numel(sea_list));

    figure(ct_fig);
    if numel(sea_list) == 1
        col_list = [0.65, 0.1, 0.88, 0];
    else
        col_list = [0, 0.65, 0.1,  0.44  0.9];
    end
    ct_lines = 0;
    for   ct_sea   =  sea_list
        ct_lines   = ct_lines + 1;
        for ct = 1:1:111
            x          = 1889+ct;
            y          = [squeeze(pic_std(ct,ct_sea,1:2))'  pic(ct,ct_sea)  squeeze(pic_std(ct,ct_sea,3:4))'];
            q_list     = [-0.05 0.25 0.5 0.75 0.95];
            input_type = 2;
            bar_width  = 1;
            if numel(sea_list) == 1
                col        = col_list(ct_fig);
            else
                col        = col_list(ct_sea);
            end
            P.do_black = 0;
            [~,RGB(:,:,ct_sea)]    = CDF_bar_quantile(x,y + offset*(ct_lines-1),col,q_list,input_type,bar_width,P);
        end
    end

    ct_lines = 0;
    for   ct_sea   = sea_list
        x          = 1890:2000;
        y          = pic(:,ct_sea);
        ct_lines   = ct_lines + 1;
        col        = RGB(2,:,ct_sea);
        if ismember(ct_sea,[1 4])
            h(ct_sea)  = CDF_histplot(x,y + offset*(ct_lines-1),st,col*.85,2);
        else
            h(ct_sea)  = CDF_histplot(x,y + offset*(ct_lines-1),st,1-(1-col)*1,2);
        end

        if ct_fig == 2
            if revision ~= 4
                plot(1890:2000,Buoy(:,ct_sea) + offset*(ct_lines-1) ,'k--','linewi',2);
            else
                plot(1890:2000,Buoy(:,ct_sea)*0 + offset*(ct_lines-1) ,'k--','linewi',2);
            end
        end
    end

    CDF_panel([1890 2000 y_lim],'','','Centered year',y_label,'fontsize',20)
    daspect([111 (y_lim(2)-y_lim(1))*1*(5/numel(sea_list))*0.8 1])
    set(gca,'ytick',y_tick,'yticklabel',y_tick_label);
    set(gca,'xtick',[1900:20:1980])

    set(gcf,'position',[.1 1 15 8],'unit','inches')
    set(gcf,'position',[.1 1 15 8],'unit','inches')
    
end 

%% *************************************************************************
% Figure 6: Explained variance
% *************************************************************************
clear;
close all;
revision = 0;   % sensitivity tests:
% 0: standard analysis
% 1: bucket identification priortizing WMO No.47 metadata
% 2: do not infer bucket SSTs from country information
% 3: nation-level analysis
% 4: excess diurnal cycles

switch revision
    case 0
        app = '';
    case 1
        app = '_sens_wmo';
    case 2
        app = '_sens_not_infer';
    case 3
        app = '_sens_nation_level';
    case 4
        app = '_sens_excess_diurnal';
        app2 = '_sens_clim_diurnal';
end

load(['DA_LME_Statistics',app,'.mat'])
if revision <= 3
    load(['All_lme_offsets_and_diurnal_amplitudes',app,'.mat'])
else
    clim = load(['All_lme_offsets_and_diurnal_amplitudes',app2,'.mat']);
    load(['All_lme_offsets_and_diurnal_amplitudes.mat']);
    da = da - clim.da;
    da_std = sqrt(da_std.^2 + clim.da_std.^2);
end

clear('y_hat')
for ct_reg = 1
    for ct_yr = 1:111
        x = da(2:end,ct_yr,ct_reg);
        y = lme(2:end,ct_yr,ct_reg);
        y_hat(:,ct_yr) = x * Tab_slp(ct_yr,ct_reg) + Tab_ipt(ct_yr,ct_reg);
    end
end
y = lme(2:end,:,ct_reg);

figure(2);clf; hold on; clear('h')
aa = quantile(y,0.75,1) - quantile(y,0.25,1);
cc = quantile(y-y_hat,0.75,1) - quantile(y-y_hat,0.25,1);
col = lines(3);
h(1) = CDF_histplot(1890:2000,aa,'-',col(1,:),3);
h(2) = CDF_histplot(1890:2000,cc(1:end),'-',col(3,:),3);
grid on;
CDF_panel([1890 2000 0 0.5],'','','Year','interquartile range of groupwise offsets (^oC)','fontsize',20)
legend(h,{'Raw','Residual'},'fontsize',20)

set(gcf,'position',[1 12 10 7],'unit','inches');
set(gcf,'position',[1 12 10 7],'unit','inches');


%% *************************************************************************
% Sensitivity Table
% *************************************************************************
clear;
ct = 0;
intvl = 20;

% Standard analysis
ct = ct + 1;
a = load('DA_LME_Statistics.mat');
Tab(:,ct,1) = [a.Tab_slp(111:-intvl:1,1)];  
Tab_ci(:,ct*2+[1:2]-2) = [a.Tab_slp_quan(111:-intvl:1,1,1)  a.Tab_slp_quan(111:-intvl:1,1,4)];

% Prefer WMO No.47 metadata
ct = ct + 1;
a = load('DA_LME_Statistics_sens_wmo.mat');
Tab(:,ct,1) = [a.Tab_slp(111:-intvl:1,1)];  
Tab_ci(:,ct*2+[1:2]-2) = [a.Tab_slp_quan(111:-intvl:1,1,1)  a.Tab_slp_quan(111:-intvl:1,1,4)];

% Not using inferred bucket SSTs
ct = ct + 1;
a = load('DA_LME_Statistics_sens_not_infer.mat');
Tab(:,ct,1) = [a.Tab_slp(111:-intvl:1,1)];  
Tab_ci(:,ct*2+[1:2]-2) = [a.Tab_slp_quan(111:-intvl:1,1,1)  a.Tab_slp_quan(111:-intvl:1,1,4)];

% Nation level analysis
ct = ct + 1;
a = load('DA_LME_Statistics_sens_nation_level.mat');
Tab(:,ct,1) = [a.Tab_slp(111:-intvl:1,1)];  
Tab_ci(:,ct*2+[1:2]-2) = [a.Tab_slp_quan(111:-intvl:1,1,1)  a.Tab_slp_quan(111:-intvl:1,1,4)];

% Anomalies relative to local climatology (version 7)
ct = ct + 1;
a = load('DA_LME_Statistics_sens_excess_diurnal.mat');
Tab(:,ct,1) = [a.Tab_slp(111:-intvl:1,1)];  
Tab_ci(:,ct*2+[1:2]-2) = [a.Tab_slp_quan(111:-intvl:1,1,1)  a.Tab_slp_quan(111:-intvl:1,1,4)];

Tab_yr = [[1990:-intvl:1880]' [2009:-intvl:1899]'];


a = num2str([Tab_yr,Tab],' %4.0f-%4.0f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f \\\\');
b = num2str(Tab_ci,' & [%6.2f, %6.2f] & [%6.2f, %6.2f] & [%6.2f, %6.2f] & [%6.2f, %6.2f] & [%6.2f, %6.2f] \\\\');

c = repmat(' ',size(a,1)*2,max(size(a,2),size(b,2)));
for ct = 1:size(a,1)
    c(ct*2-1,1:size(a,2)) = a(ct,:);
    c(ct*2,1:size(b,2))   = b(ct,:);
end

c