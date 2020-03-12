function [WM,ST,NUM] = LME_correct(P)
    
    % *********************************************************************
    % O/I 
    % *********************************************************************
    if ~exist('env','var')
        env = 1;            % 1 means on odyssey
    end

    if ~isfield(P,'connect_Kobe')
        P.connect_Kobe = 0;
    end

    if ~isfield(P,'do_add_JP')
        P.do_add_JP = 0;
    end

    if ~isfield(P,'do_rmdup')
        P.do_rmdup = 0;
    end

    % *********************************************************************
    % Read the LME outputs 
    % *********************************************************************
    dir_lme = LME_OI('LME_output');  
    file_lme = [dir_lme,'LME_',P.save_lme,'.mat'];
    lme = load(file_lme,'out','out_rnd');

    % *********************************************************************
    % Set Parameters for gridding 
    % *********************************************************************
    reso_x = P.reso_x;
    reso_y = P.reso_y;
    yr_list = P.yr_list;
    yr_num = numel(yr_list);

    if P.en == 0
        P.do_individual = 0;
    end

    mon_list = 1:12;

    if P.en == 0 || P.do_individual == 1
        P.do_random = 0;
    else
        P.do_random = 1;
    end

    % *********************************************************************
    % Assigning effects to be corrected 
    % *********************************************************************
    clear('E')
    E = LME_correct_assign_effect(lme,P);

    % *********************************************************************
    % Initialize the correction 
    % *********************************************************************
    N = double(P.do_region) + double(P.do_season) + double(P.do_decade) + 3;

    clear('WM','ST','NUM')
    if P.en == 0 && P.do_individual == 0
        WM  = nan(360/reso_x,180/reso_y,N,12,yr_num);
        ST  = nan(360/reso_x,180/reso_y,N,12,yr_num);
        NUM = nan(360/reso_x,180/reso_y,N,12,yr_num);
    else
        WM  = nan(360/reso_x,180/reso_y,12,yr_num);
        ST  = nan(360/reso_x,180/reso_y,12,yr_num);
        NUM = nan(360/reso_x,180/reso_y,12,yr_num);
    end

    % *********************************************************************
    % Start the correction 
    % *********************************************************************
    for yr = yr_list
        for mon = mon_list

            disp(['En :',num2str(P.en),' Starting Year ',...
                                      num2str(yr),'  Month ',num2str(mon)])

            % *************************************************************
            % Read in files 
            % *************************************************************
            disp('Reading data ...')
            clear('DATA')
            PP = P;  PP.yr  = yr;  PP.mon = mon;   PP.mute_read = 1;
            DATA = LME_pair_function_read_data(PP);
            DATA = LME_function_preprocess_SST_method(DATA);

            % remove non-ship and non bucket, ERI, Hull, or missing SSTs
            if strcmp(P.method,'Ship')
                l_rm = ~ismember(DATA.C0_SI_4,[-1 0 1 3 13 14]);
                var_list = fieldnames(DATA);
                for var = 1:numel(var_list)
                    if ~ismember(var_list{var},{'C0_ID','C0_CTY_CRT','DCK'})
                        eval(['DATA.',var_list{var},'(l_rm) = [];']);
                    else
                        eval(['DATA.',var_list{var},'(l_rm,:) = [];']);
                    end
                end
            end
            kind = [DATA.DCK DATA.C0_SI_4'];

            if strcmp(P.type,'Bucket_vs_ERI')
                if isfield(P,'all_ERI_in_one_group')
                    if P.all_ERI_in_one_group == 1
                        kind(kind(:,end) == 1,:) = 1;
                    end
                end
            end

            if strcmp(P.type,'ERI_vs_Bucket') || strcmp(P.type,'ERIex_vs_Bucket')
                if isfield(P,'all_BCK_in_one_group')
                    if P.all_BCK_in_one_group == 1
                        kind(kind(:,end) == 0,:) = 0;
                    end
                end
            end
              
            % *************************************************************
            % Assigning Effects
            % *************************************************************
            disp(['Assigning effects ...'])
            clear('ID')
            if P.do_region == 1
                ID.rid = LME_lme_effect_regional(DATA.C0_LON,DATA.C0_LAT,5);
            end

            if P.do_season == 1
                ID.sid = LME_lme_effect_seasonal(DATA.C0_LAT,DATA.C0_MO);
            end

            if P.do_decade == 1
                ID.did = LME_lme_effect_decadal(DATA.C0_YR,P);
            end
            clear('mx','my')

            % *************************************************************
            % Applying Correction 
            % *************************************************************
            disp(['Applying Correction ...'])
            clear('CORR')
            CORR = LME_correct_find_corr(DATA,E,P,ID,kind,lme.out.unique_grp);
                % The correction follows: 1. all 2.fix 3.reg 4.dcd 5.sea
                % the output is correction but not bias

            % *************************************************************
            % Gridding data 
            % *************************************************************
            disp(['Gridding data ...'])
            yr_id = yr-yr_list(1)+1;

            if P.en == 0 && P.do_individual == 0

                for ct = 1:size(CORR.sst_correction,1)+1
                    
                    if ct == 1
                        SST_in = DATA.C0_SST- DATA.C0_OI_CLIM;
                    else
                        SST_in = nansum([DATA.C0_SST- DATA.C0_OI_CLIM; CORR.sst_correction(ct-1,:)],1);
                    end

                    [WM(:,:,ct,mon,yr_id),ST(:,:,ct,mon,yr_id),NUM(:,:,ct,mon,yr_id)] = ...
                    LME_function_gridding(DATA.C0_LON,DATA.C0_LAT,[],SST_in,[],reso_x,reso_y,[],2,[],[],[]);
                end
                
            else

                SST_in = nansum([DATA.C0_SST- DATA.C0_OI_CLIM; CORR.sst_correction(1,:)],1);

                [WM(:,:,mon,yr_id),ST(:,:,yr_id),NUM(:,:,mon,yr_id)] = ...
                LME_function_gridding(DATA.C0_LON,DATA.C0_LAT,[],SST_in,[],reso_x,reso_y,[],2,[],[],[]);

            end
        end
    end
end
