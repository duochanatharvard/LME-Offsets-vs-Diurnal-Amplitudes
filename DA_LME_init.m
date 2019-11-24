% DA_LME_init(dir_data)
% Input dir_data is the directory for storing data

function DA_LME_init(dir_data)

    % ********************************************
    % Make directories
    % ********************************************
    dir_code = [pwd,'/'];

    cd(dir_code)
    dir_home_ICOADS3 = [dir_data,'ICOADS3/'];
    dir_home_diurnal = [dir_data,'DIURNAL/'];
    dir_home_LME     = [dir_data,'LME_intercomparison/'];
    
    save('DA_LME_directories.mat','dir_home_ICOADS3',...
         'dir_home_diurnal','dir_home_LME','dir_data','dir_code','-v7.3');
     
    mkdir(dir_data) 
    cd(dir_data)
    mkdir Miscellaneous
     
    mkdir(dir_home_ICOADS3)
    cd(dir_home_ICOADS3)
    mkdir ICOADS_QCed
    mkdir ICOADS_Tracks_Kent
    
    mkdir(dir_home_diurnal)
    cd(dir_home_diurnal)
    mkdir Step_01_Ship_Signal
    mkdir DATA_for_figures
    
    mkdir(dir_home_LME)
    cd(dir_home_LME)
    mkdir Step_01_All_Pairs
    mkdir Step_02_Screen_Pairs
    mkdir Step_03_Binned_Pairs
    mkdir Step_04_LME_output
    mkdir DATA_for_figures
    
    cd(dir_code)
    copyfile('Miscellaneous',[dir_data,'Miscellaneous'])
    
end
