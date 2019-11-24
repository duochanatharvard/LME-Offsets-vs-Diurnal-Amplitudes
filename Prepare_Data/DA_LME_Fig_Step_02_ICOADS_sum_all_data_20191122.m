% This code merges all diurnal amplitudes and LME offsets into a single
% Table

function [da,lme,da_std,lme_std,grp,buoy] = DA_LME_Fig_Step_02_ICOADS_sum_all_data_20191122(ct_reg,yr_start,yr_end,ct_sea,app_da,app_lme)

    % *********************************************************************
    % plot for observations
    % *********************************************************************
    dir_save  = LME_OI('data4figure');
    O = DA_LME_function_get_region_name(ct_reg);
    file_load = [dir_save,'DATA_Plot_DA&LME_',num2str(yr_start),'_',...
        num2str(yr_end),'_',O.region_name,'_da_',app_da,'_lme_',app_lme,'.mat'];
    r1p1 = load(file_load);
    disp(file_load)
    
    % *********************************************************************
    % Subset data of interest
    % *********************************************************************
    clear('da','lme','da_std','lme_std','grp')
    da     = squeeze([r1p1.DATA_pic.DA_amp_bck(:,ct_sea);     r1p1.DATA_pic.DA_amp_ERI(1,ct_sea)]);
    da_std = squeeze([r1p1.DATA_pic.DA_amp_std_bck(:,ct_sea); r1p1.DATA_pic.DA_amp_std_ERI(1,ct_sea)]);   

    lme     = squeeze([r1p1.DATA_pic.LME_bck(:,ct_sea);       r1p1.DATA_pic.LME_ERI(1,ct_sea)]);
    lme_std = squeeze([r1p1.DATA_pic.LME_std_bck(:,ct_sea);   r1p1.DATA_pic.LME_std_ERI(1,ct_sea)]); 
    
    buoy = [r1p1.DATA_pic.DA_amp_buoy(1,ct_sea)];
    
    grp   = [r1p1.DATA_pic.group_bck; 1 1 1];    
end