clear;

% *************************************************************************
% Set Directries
% *************************************************************************
addpath(genpath(pwd));

% *************************************************************************
% Parsing Parameters
% *************************************************************************
for ct_yr = 1:111                   % Change if running for different years             % TODO
    
    yr_start = 1879+ct_yr;
    yr_end = 1898+ct_yr;
    
    use_kent_tracks     = 0;
    use_diurnal_points  = 1;
    P.relative          = 'mean_SST';
    
    
    % *************************************************************************
    % Set Parameters
    % *************************************************************************
    % parameters for controling the run ---------------------------------------
    P.varname   = 'SST';         % Variable to be analyzed
    P.method    = 'Bucket';      % Subset of data according to meansurement method
    P.yr_list   = yr_start:yr_end;     % Length of the analysis                         
    P.mon_list  = 1:12;          % Seasonal coverage of the analysis
    P.select_region = 1:17;      % To confine the analysis in certain regions?
    %  1. Sub-NA   2. Sub-WNP   3. Sub-EMP  4. Ex-NA   5. Ex-NP
    %  6. Medt     7. TA        8. WTP      9. CTP     10.ETP
    % 11. Indian  12. Sub-SA   13. SIO     14. SP      15. SO
    % 16. SPole   17. NPole
    P.mute_read = 1;             % Whether to turn off the output for debugging?
    P.restart     = 0;           % Whether to restart the summing of bins?
    
    % parameters for the LME model (binning) ----------------------------------
    P.use_kent_tracks     = use_kent_tracks;         % Only use data from tracked ships in the LME?  
    P.use_diurnal_points  = use_diurnal_points;      % only use data that computed diurnal cycles.  
    P.use_fundemental_SST = 0;                       % use fundemental SST in the LME analysis.     
    P.do_connect      = 1;        % Connect decks that have the same discreptions?
    P.connect_Kobe    = 1;        % Treat deck 119 and 762 as 118?
    P.do_simple       = 0;        % Whether to use the simple model?
    P.do_region       = 1;        % whether turn on the regional effect?
    P.do_season       = 1;        % whether turn on the seasonal effect?
    P.do_decade       = 0;        % whether turn on the decadal effect?
    P.yr_start        = yr_start; % when decadal is on, which is the first year?           
    P.yr_interval     = 10;       % when decadal is on, length of the yearly bins
    P.key             = 6000;
    P.buoy_diurnal    = 1;        % Use diurnal cycle estimated from buoy measurements.
    P.nation_id       = [1 2 4];  % which part of ID is regarded as nation - method
    
    % parameters for the LME model (running) ----------------------------------
    P.do_sampling         = 1000;
    P.target              = 0;
    P.do_hierarchy        = 1;
    P.do_hierarchy_random = 1;
    P.var_list = {'C0_YR','C0_MO','C0_UTC','C0_LON','C0_LAT',...
        'C1_DCK','C0_SST','C0_OI_CLIM','C0_SI_4','C0_CTY_CRT','C98_UID'};
    
    % parameters for input and output -----------------------------------------
    if strcmp(P.method,'Bucket')                  % This part should not be touched
        P.type      = 'Bucket_vs_ERI';                      % Which type of pairs to read from step 1 in paring
        % P.save_app  = 'Bucket_vs_ERI_in_one_group';         % file name for screened pairs for step 2 in paring
        
        P.all_ERI_in_one_group = 1;                  % collapse all ERIs into one group            
        
        if P.use_kent_tracks == 0 && P.use_diurnal_points == 1,
            P.save_app       = 'Bucket_vs_ERI_in_one_group';
            P.save_app       = [P.save_app,'_diurnal_points_',P.relative];
        end
    end
    
    % parameters for diurnal cycle analyses -----------------------------------
    P.diurnal_QC      = 1;             % use different QCs for the diurnal cycles
    
    
    % *************************************************************************
    % Set output filenames
    % *************************************************************************
    % file name for screened pairs --------------------------------------------
    P.save_sum  = [P.save_app,'_',num2str(P.yr_list(1)),'_',num2str(P.yr_list(end))];
    P.save_sum  = [P.save_sum,'_Full_SST'];
    
    
    % file name for bined pairs -----------------------------------------------
    P.save_bin  = [P.save_sum,'_Global'];
    P.save_lme  = P.save_bin;                                                  % file name for LME
    clear('temp','ct','c_temp')
    
    
    % *************************************************************************
    % Bin pairs before fitting the LME model
    % *************************************************************************
    if 1,    [BINNED,W_X] = LME_lme_bin(P);                   end
    
    % *************************************************************************
    % Fitting the LME model
    % *************************************************************************
    if 1,      LME_lme_fit_hierarchy(P);                      end
    if 0,      LME_lme_fit(P);                                end
end