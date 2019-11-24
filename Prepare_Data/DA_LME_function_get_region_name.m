
function O = DA_LME_function_get_region_name(ct_reg)

    switch ct_reg,
        case 1,      % Tropics
            region_name   = 'Tropics';
            region_name_m = 'Trop';
            reg_list      = [7 8 9 10 11];
            sea_list      = [1 3];
        case 2,      % Tropical Atlantic
            region_name   = 'Tropical_Atlantic';
            region_name_m = 'Trop_Atl';
            reg_list      = [7];
            sea_list      = [1 3];
        case 3,      % Subtropics
            region_name   = 'Subtropics';
            region_name_m = 'SubT';
            reg_list      = [1 2 3];
            sea_list      = [1 3] + 4;
        case 4,      % Subtropical Atlantic
            region_name   = 'Subtropical_Atlantic';
            region_name_m = 'SubT_Atl';
            reg_list      = [1];
            sea_list      = [1 3] + 4;
        case 5,      % Tropical Pacific
            region_name   = 'Tropical_Pacific';
            region_name_m = 'Trop_Pac';
            reg_list      = [8 9 10];
            sea_list      = [1 3];
        case 6,      % Subtropical Pacific
            region_name   = 'Subtropical_Pacific';
            region_name_m = 'SubT_Pac';
            reg_list      = [2 3];
            sea_list      = [1 3] + 4;    
        case 7,      % Extratropics
            region_name   = 'Extratropics';
            region_name_m = 'ExT';
            reg_list      = [4 5];
            sea_list      = [1 3] + 8;
        case 8,      % Extratropical Atlantic
            region_name   = 'Extratropical_Atlantic';
            region_name_m = 'ExT_Atl';
            reg_list      = [4];
            sea_list      = [1 3] + 8;
        case 9,      % Extratropical Pacific
            region_name   = 'Extratropical_Pacific';
            region_name_m = 'ExT_Pac';
            reg_list      = [5];
            sea_list      = [1 3] + 8;
        case 10,      % Tropical Indian Ocean
            region_name   = 'Indian_Ocean';
            region_name_m = 'Ind_Ocn';
            reg_list      = [11];
            sea_list      = [1 3];     
    end

    O.region_name = region_name;
    O.region_name_m = region_name_m;
    O.reg_list = reg_list;
    O.sea_list = sea_list;
end