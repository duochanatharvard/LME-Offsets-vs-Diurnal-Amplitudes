% Output / Input managements

function output = DIURNAL_OI(input)

    if strcmp(input,'home')
        % home directory for Diurnal analysis
        load('DA_LME_directories.mat','dir_home_diurnal')
        output = dir_home_diurnal; 

    elseif strcmp(input,'read_raw')
        % home directory for ICOADS3
        load('DA_LME_directories.mat','dir_home_ICOADS3')
        output = dir_home_ICOADS3;

    elseif strcmp(input,'ICOADS3')
        output = [DIURNAL_OI('read_raw'),'ICOADS_QCed/'];

    elseif strcmp(input,'kent_track')
        output = [DIURNAL_OI('read_raw'),'ICOADS_Tracks_Kent/'];

    elseif strcmp(input,'ship_signal')
        output = [DIURNAL_OI('home'),'Step_01_Ship_Signal/'];

    elseif strcmp(input,'Mis')
        load('DA_LME_directories.mat','dir_data')
        output = [dir_data,'Miscellaneous/'];

    elseif strcmp(input,'data4figure')
        output = [DIURNAL_OI('home'),'DATA_for_figures/'];
        
    end
end
