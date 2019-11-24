% How to run this code
clear;

% Merge LME and DA data
app_da  = 'raw';
app_lme = 'diurnal_points_mean_SST';
for ct_yr = 1:1:111
    for ct_reg = [1 3 7]
        DA_LME_Fig_Step_01_prepare_data(1879+ct_yr,1898+ct_yr,ct_reg,1,1,app_da,app_lme);
    end
end

% Merge across years and region-season combinations
app_da  = 'raw';
app_lme = 'diurnal_points_mean_SST';
ct = 0;
clear('P_da','P_da_std','P_lme','P_lme_std')
for ct_case = 1:5
    for ct_prd = 1:1:111
        clear('da','lme','da_std','lme_std','grp')
        switch ct_case,
            case 1,
                [da,lme,da_std,lme_std,grp,buoy] = DA_LME_Fig_Step_02_ICOADS_sum_all_data_20191122(1,1879+ct_prd,1898+ct_prd,1,app_da,app_lme);
            case 2,
                [da,lme,da_std,lme_std,grp,buoy] = DA_LME_Fig_Step_02_ICOADS_sum_all_data_20191122(3,1879+ct_prd,1898+ct_prd,2,app_da,app_lme);
            case 3,
                [da,lme,da_std,lme_std,grp,buoy] = DA_LME_Fig_Step_02_ICOADS_sum_all_data_20191122(3,1879+ct_prd,1898+ct_prd,3,app_da,app_lme);
            case 4,
                [da,lme,da_std,lme_std,grp,buoy] = DA_LME_Fig_Step_02_ICOADS_sum_all_data_20191122(7,1879+ct_prd,1898+ct_prd,2,app_da,app_lme);
            case 5,
                [da,lme,da_std,lme_std,grp,buoy] = DA_LME_Fig_Step_02_ICOADS_sum_all_data_20191122(7,1879+ct_prd,1898+ct_prd,3,app_da,app_lme);
        end
        ct = ct + 1;
        P_da{ct}       = [da       grp];
        P_da_std{ct}   = [da_std   grp];
        P_lme{ct}      = [lme      grp];
        P_lme_std{ct}  = [lme_std  grp];
        Buoy(ct_prd,ct_case) = buoy;
    end
end

NN = 3;
P_da_table       = CDC_merge_group(P_da,NN);
P_da_std_table   = CDC_merge_group(P_da_std,NN);
P_lme_table      = CDC_merge_group(P_lme,NN);
P_lme_std_table  = CDC_merge_group(P_lme_std,NN);

grp       = P_da_table(:,end-NN+1:end);
N_grp     = size(grp,1);
da        = reshape(P_da_table(:,1:555),N_grp,111,5);
da_std    = reshape(P_da_std_table(:,1:555),N_grp,111,5);
lme       = reshape(P_lme_table(:,1:555),N_grp,111,5);
lme_std   = reshape(P_lme_std_table(:,1:555),N_grp,111,5);

load('DA_LME_directories.mat','dir_data')
file_save = [dir_data,'All_lme_offsets_and_diurnal_amplitudes.mat'];
save(file_save,'da','da_std','lme','lme_std','grp','Buoy','-v7.3')