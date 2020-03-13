% Statistics of number of SSTs from individual nation-deck groups.
% Result of this script is provided as "Stats_HM_SST_Bucket_deck_level_1.mat"

reso_x = 5;
reso_y = 5;
yr_start = 1850;
yr_end   = 2014;

% *************************************************************************
% Set Parameters
% *************************************************************************
P.varname   = 'SST';         % Variable to be analyzed
P.method    = 'Ship';      % Subset of data according to meansurement method
P.yr_list   = yr_start:yr_end;     % Length of the analysis                          % TODO

P.mon_list            = 1:12;    % Seasonal coverage of the analysis
P.select_region       = 1:17;    % To confine the analysis in certain regions?
P.use_kent_tracks     = 0;       % Only use data from tracked ships in the LME?  % TODO
P.use_diurnal_points  = 0;       % only use data that computed diurnal cycles.   % TODO
P.use_fundemental_SST = 0;       % use fundemental SST in the LME analysis.      % TODO

P.do_connect      = 1;        % Connect decks that have the same discreptions?
P.connect_Kobe    = 1;        % Treat deck 119 and 762 as 118?

% *************************************************************************
% Read the observational file 
% *************************************************************************
dir_lme = LME_OI('LME_output');
file_lme = [dir_lme,'LME_Ship_vs_Ship_all_measurements_1850_2014_Full_SST_Global.mat'];

lme = load(file_lme,'out');
kind_ref = lme.out.unique_grp;
clear('lme')

% *************************************************************************
% Start the statistics
% *************************************************************************
Stats_map = zeros(360/reso_x, 180/reso_y, numel(P.yr_list), size(kind_ref,1));
Stats_glb = zeros(numel(P.yr_list),12,size(kind_ref,1));

for yr = P.yr_list

    disp(['Starting Year ',num2str(yr)])
    yr_id = yr-P.yr_list(1)+1;

    % *********************************************************************
    % Read in files
    % *********************************************************************
    disp('Reading data ...')
    clear('LON','LAT','kind','MON')
    LON  = [];
    LAT  = [];
    kind = [];
    MON  = [];
    for mon = 1:12
        PP = P;   PP.yr = yr;  PP.mon = mon;  PP.mute_read = 1;
        P1 = LME_pair_function_read_data(PP);    clear('PP')
        P1 = LME_function_preprocess_SST_method(P1);
        LON = [LON P1.C0_LON];
        LAT = [LAT P1.C0_LAT];
        MON = [MON P1.C0_MO];
        kind = [kind; [P1.DCK P1.C0_SI_4']];
    end

    % *********************************************************************
    % STATISTICS!!!
    % *********************************************************************
    disp('Do Statistics ...')
    [kind_uni,~,J_kind] = unique(kind,'rows');
    
    for ct = 1:size(kind_uni,1)
        
        [~,dck_id] = ismember(kind_uni(ct,:),kind_ref,'rows');

        if nnz(dck_id)

            % -------------------------------------------------------------
            % Stats for annually map 
            % -------------------------------------------------------------
            temp  = [LON(J_kind == ct)' LAT(J_kind == ct)'];
            count = hist2d(temp, [0,reso_x,360; -90,reso_y,90]);
            Stats_map(:,:,yr_id,dck_id) = count;

            % -------------------------------------------------------------
            % Stats for monthly global 
            % -------------------------------------------------------------
            for mon = 1:12
                Stats_glb(yr_id,mon,dck_id) = nnz(MON(J_kind == ct) == mon);
            end
        end
    end
end

unique_grp = kind_ref;
dir_save   = LME_OI('home');
file_save  = [dir_save,'Stats_All_ships.mat'];
save(file_save,'Stats_map','Stats_glb','unique_grp','Ratio','-v7.3')
