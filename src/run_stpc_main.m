function STPC_p_results = run_stpc_main(basicInfo, STPCInfo, slepian_results, ...
    mssa_results);
% RUN_STPC_MAIN  Do STPC for reconstruction on Slepian and gap-filling MSSA results
%
% Usage:
%   STPC_results = run_stpc_main(basicInfo, STPCInfo, slepian_results, ...
%        mssa_results);
%
% Inputs:
%   basicInfo   - struct, same definition as in run_slepian_main:
%       .Code_Version
%       .Area
%       .Lwindow
%       .c11cmn
%       .ddir1, .ddir2
%       .ifilesRoot
%
%   slepianInfo - struct, same as in run_slepian_main:
%       .Dataproduct
%       .buffer_str
%       .Max_S
%       .Radius
%       .use_institu     % (recommended) list of institution names
%
%   mssaInfo    - struct, corresponding to original MSSA_Information.mat:
%       .S_bou          % lower/upper lambda thresholds for S
%       .N_bou          % lower/upper cumulative EV thresholds for N
%       .Turning_number % number of change points to detect
%       .M_gap_default  % window size M for MSSA gap filling
%       .N_gap_default  % number of RCs N for MSSA gap filling
%       .redo           % optional flag to force recomputation (not used here)
%
%   instResults - struct array, length = intit_num or intit_num+1
%       Each element is the output of run_slepian_main for one institution:
%       .Dataproduct
%       .TH, .XY, .BasinArea
%       .S, .S_sig, .S_shannon, .CC, .V, .Radius
%       .time, .fill, .grid, .date, .ocean
%
% Output:
%   STPC_results - struct, containing:
%       .mssa_Sort_sn, .A5_Centers_Fit
%       .common_use_months
%       .turn_V_four, .turn_MSSA_four
%       .M_rec, .N_gap, .M_gap
%       .MSSAf_evalues_sumup, .RC_RMS_Seq, .V, .Max_S, .S_ol, .BasinArea
%       .lon, .lat, .r_record
%
% References: see run_slepian_main.

arguments
    basicInfo       struct
    STPCInfo        struct
    slepian_results     struct
    mssa_results     struct
end

%% 0. Unpack basic information
Code_Version = basicInfo.Code_Version;
Area         = basicInfo.Area;
Lwindow      = basicInfo.Lwindow;
c11cmn       = basicInfo.c11cmn;
land_or_ocean= basicInfo.land_or_ocean;
Max_S        = basicInfo.Max_S;
ddir1        = basicInfo.ddir1;      % output results directory (Results_*)
ddir2        = basicInfo.ddir2;      % output figure directory (Figures_*)
ifilesRoot   = basicInfo.ifilesRoot; % root directory of original IFILES
redo         = basicInfo.redo;
plotProcess  = basicInfo.plotProcess;
saveAddData  = basicInfo.saveAddData;

Turning_number  = STPCInfo.Turning_number;
S_bou           = STPCInfo.S_bou;
N_bou           = STPCInfo.N_bou;
p_use           = STPCInfo.p_use;

Attach       = mssa_results.Attach;
use_institu  = mssa_results.use_institu;
intit_num    = mssa_results.intit_num;
M_rec        = mssa_results.M_rec;
N_gap        = mssa_results.N_gap;
S            = mssa_results.S;
S_ol         = mssa_results.S_ol;
mssa_Sort    = mssa_results.mssa_Sort;

fillM_deltacoffs    = mssa_results.fillM_deltacoffs;

V               = slepian_results(1).V;
lon             = slepian_results(1).grid.lon;
lat             = slepian_results(1).grid.lat;
XY              = slepian_results(1).XY;
XY_ori          = slepian_results(1).XY_ori;
XY_buffer       = slepian_results(1).XY_buffer;
BasinArea       = slepian_results.BasinArea;
r_record        = slepian_results(1).grid.r_record;
fill_months     = slepian_results(1).date.fill_months;

S_bou           = STPCInfo.S_bou;
N_bou           = STPCInfo.N_bou;
p_use           = STPCInfo.p_use;
S_selection     = STPCInfo.S_selection;
N_selection     = STPCInfo.N_selection;
sAdaptive       = strcmp(S_selection.mode, 'stpc');
nStpc           = strcmp(N_selection.mode, 'stpc');
nWcorr          = strcmp(N_selection.mode, 'wcorr');
nAdaptive       = nStpc || nWcorr;

stpcBaseName = ['MainSTPC_V3_' Attach.Attach_ALL];
datasetsFile = fullfile(ddir1, [stpcBaseName '_datasets.mat']);
turningPointFile = fullfile(ddir1, [stpcBaseName '_TP.mat']);

%% should we redo?

% redo=true; % debug only; keep disabled so config.redo controls cache use.
if ~redo && exist(datasetsFile,'file')

    load(datasetsFile, ...
        'STPC_p_results')
    disp(['STPC: Loading ', datasetsFile ...
        ' for turning points of S and N.']);

    return

end

%% Basic calculation

[lonlon,latlat] = meshgrid(lon, lat);

% if the longitude ranges from 0 to 360. Make it consistent to [-180,180]
if XY(1,1)>180
    XY(:,1)=XY(:,1)-360;
end

[in,on] = check_polygon_in(XY,lonlon,latlat); %#ok<ASGLU>

if ~sum(in,'all')
    error('Be caureful on your chosen XY while checking the inner points.')
end

%% Compute area weights (lat/lon cell area)
c11cmn_area = zeros(length(lat), length(lon));
for i = 1:length(lat)
    for j = 1:length(lon)
        c11cmn_area(i,j) = areaquad(lat(i)-0.5, lon(j)-0.5, ...
                                    lat(i)+0.5, lon(j)+0.5);
    end
end

%% 2. Collect institution results and compute mean product

disp('==== Do STPC calculation ====');

% redo=true; % debug only; keep disabled so config.redo controls cache use.
turningPointTimer = tic;
if ~redo && exist(turningPointFile,'file')
    % directly load this variables

    disp(['STPC: Loading ', turningPointFile ...
        ' for turning points of S and N.']);
    load(turningPointFile, ...
        'turn_V_four','turn_MSSA_four','MSSA_evalues', 'used_sta_SN');

else
    used_sta_SN=4; % what kind of turning point do you want? (default is RSS)

%     plotProcess=true;
    %% 1. Determine turning points of Slepian eigenvalues V
    if sAdaptive
        turn_V_four = det_S_TP(V, Max_S, S_bou, ...
            Turning_number, used_sta_SN, Attach, plotProcess);
    else
        turn_V_four = repmat(S, 1, 4);
    end

    %% 2. Determine turning points of MSSA eigenvalues
    N_selection_diagnostics = struct('mode', N_selection.mode);
    if nStpc
        [turn_MSSA_four, MSSA_evalues_sumup, MSSA_evalues] = det_N_TP( ...
            mssa_results, N_bou, Turning_number, used_sta_SN, Attach, ...
            plotProcess);
    elseif nWcorr
        [turn_MSSA_four, MSSA_evalues_sumup, MSSA_evalues, ...
            N_selection_diagnostics] = det_N_WCorr( ...
            mssa_results, N_bou, Turning_number, Attach, plotProcess, ...
            N_selection);
    else
        fixedN = N_selection.value;
        turn_MSSA_four = cell(1, S_ol);
        for ss = 1:S_ol
            turn_MSSA_four{ss} = repmat(fixedN, 1, 4);
        end
        MSSA_evalues_sumup = [];
        MSSA_evalues = [];
    end

    %% 3. Plot MSSA turning points and combined (V + MSSA) figure
    if plotProcess && sAdaptive && nStpc
        plot_mssa_evalues_and_turning_points(mssa_results, V, Max_S, S_bou, ...
            N_bou, turn_V_four, turn_MSSA_four, MSSA_evalues_sumup, ...
            Turning_number, used_sta_SN, Attach);
    end

    save(turningPointFile, ...
        'turn_V_four','turn_MSSA_four','MSSA_evalues', 'used_sta_SN', ...
        'N_selection_diagnostics')

end
turningPointSeconds = toc(turningPointTimer);

%% 4. Calculate for different S and N at significance level of p

Square_Seq=cell([]);
STPC_p_results=struct();

% V3: the full SVD/RC decomposition is independent of p, S candidate and
% N truncation.  Build it once at the maximum candidate S/N; every p still
% receives its own significance decision from the cached p-values.
[decompositionCache, sharedDecompositionSeconds] = ...
    prepare_stpc_decomposition_cache(mssa_results, turn_V_four, ...
    turn_MSSA_four, used_sta_SN, Turning_number, p_use(1));

numSCandidates = size(turn_V_four, 1);
numNCandidates = size(turn_MSSA_four{1}, 1);
fullGridDiagnostics = plotProcess && sAdaptive && nAdaptive;

for p=1:length(p_use)

    Noise_SigLev=p_use(p);

    if fullGridDiagnostics
        % Plotting requires all full candidate grids.  Preserve the legacy
        % diagnostic behavior, but still reuse the shared decomposition.
        candidateOptions = struct('mode', 'legacy', ...
            'decompositionCache', decompositionCache, ...
            'plotCandidateEigenvalues', true, ...
            'plotSelectedDecomposition', false);
        [STPC_p_SN_results, ~, candidateTiming] = run_stpc_pvalue( ...
            basicInfo, mssa_results, slepian_results, turn_V_four, ...
            turn_MSSA_four, Noise_SigLev, Turning_number, MSSA_evalues, ...
            used_sta_SN, Attach, candidateOptions);

        analysisTimer = tic;
        RMS_SN_Results = run_stpc_temspa_RMS( ...
            Turning_number, turn_V_four, turn_MSSA_four, used_sta_SN, ...
            fill_months, in, c11cmn_area, STPC_p_SN_results);
        candidateTiming.full_rms_analysis_seconds = ...
            candidateTiming.full_rms_analysis_seconds + toc(analysisTimer);
    elseif numSCandidates * numNCandidates > 1
        % A fixed S or fixed N reduces the search to one dimension.  Keep
        % only the WRMS needed for selection; do not create the legacy
        % spatial/temporal K-by-K diagnostic figures.
        screenOptions = struct('mode', 'screen', ...
            'decompositionCache', decompositionCache, ...
            'plotCandidateEigenvalues', plotProcess, ...
            'plotSelectedDecomposition', false);
        [~, RMS_SN_Results, candidateTiming] = run_stpc_pvalue( ...
            basicInfo, mssa_results, slepian_results, turn_V_four, ...
            turn_MSSA_four, Noise_SigLev, Turning_number, MSSA_evalues, ...
            used_sta_SN, Attach, screenOptions);
        STPC_p_SN_results = [];
    else
        % Both S and N are fixed: no candidate screening is necessary.
        STPC_p_SN_results = [];
        RMS_SN_Results = struct('s_ss_nn_rms', nan(5, 1, 1));
        candidateTiming = local_zero_timing();
    end

    %% determine the finally suitable S and N at significance level of p

    % which RMS is used as creterion
    use_rms = RMS_SN_Results.s_ss_nn_rms; % [disp, S, N]
    use_rms_title         = 'WRMS';  % Weighted RMS

    if fullGridDiagnostics
        plot_stpc_spatial_temporal_RMS( ...
            STPC_p_SN_results, Turning_number, Area, XY_ori, XY_buffer, ...
            lon, lat, Noise_SigLev, use_rms, use_rms_title, Attach)
    end

    % which disp_data is interpreted as "signal"
    disp_sn = 3; % 1 for MSSA EWH, 2 for MSSA-detected noise
    % 3 for STPC EWH, 4 for periodic terms of STPC EWH, 5 for residual
    % terms of STPC EWH

    % Return matrix used to determine final S and N
    Square_Seq{p} = squeeze(use_rms(disp_sn, :, :));

    if numSCandidates * numNCandidates == 1
        p_S(p) = 1;
        p_N(p) = 1;
    else
        [~, maxLinearIndex] = max(Square_Seq{p}, [], 'all');
        [p_S(p), p_N(p)] = ind2sub( ...
            [numSCandidates, numNCandidates], maxLinearIndex);
    end

    % Reconstruct complete three-dimensional products only once after the
    % winning pair is known.  This also makes the decomposition figure
    % correspond to the final S/N, rather than the last grid candidate.
    finalOptions = struct('mode', 'selected', ...
        'decompositionCache', decompositionCache, ...
        'selectedSN', [p_S(p), p_N(p)], ...
        'plotCandidateEigenvalues', false, ...
        'plotSelectedDecomposition', plotProcess);
    [selectedResults, ~, finalTiming] = run_stpc_pvalue( ...
        basicInfo, mssa_results, slepian_results, turn_V_four, ...
        turn_MSSA_four, Noise_SigLev, Turning_number, MSSA_evalues, ...
        used_sta_SN, Attach, finalOptions);
    STPC_p_result = selectedResults{p_S(p),p_N(p)};
    
    STPC_p_results(p).name=[STPC_p_result.name '_p' ...
        num2str(Noise_SigLev*100) '%'];
    STPC_p_results(p).c11cmn=STPC_p_result.c11cmn;
    STPC_p_results(p).lon=STPC_p_result.lon;
    STPC_p_results(p).lat=STPC_p_result.lat;
    STPC_p_results(p).process=STPC_p_result.process;
    fill_dates_vec=datevec(STPC_p_result.fill_dates);
    STPC_p_results(p).fill_dates=fill_dates_vec(:,1:3);
    use_dates_vec=datevec(STPC_p_result.use_dates);
    STPC_p_results(p).use_dates=use_dates_vec(:,1:3);
    missing_dates_vec=datevec(STPC_p_result.missing_dates);
    STPC_p_results(p).missing_dates=missing_dates_vec(:,1:3);

    % store area-weighted and gridded EWH (cm)
    STPC_p_results(p).EWH.unit='cm';
    STPC_p_results(p).EWH.fillM_Slepian=STPC_p_result.EWH.fillM_Slepian;
    STPC_p_results(p).EWH.fillM_MSSA=STPC_p_result.EWH.fillM_MSSA;
    STPC_p_results(p).EWH.fillM_STPC=STPC_p_result.EWH.fillM_STPC;
    STPC_p_results(p).EWH.fillM_Grid_Slepian=STPC_p_result.EWH.fillM_Grid_Slepian;
    STPC_p_results(p).EWH.fillM_Grid_MSSA=STPC_p_result.EWH.fillM_Grid_MSSA;
    STPC_p_results(p).EWH.fillM_Grid_STPC=STPC_p_result.EWH.fillM_Grid_STPC;
    STPC_p_results(p).EWH.fillM_Grid_MSSA_bothfail_uptoS_RMS= ...
        STPC_p_result.EWH.fillM_Grid_MSSA_bothfail_uptoS_RMS;

    % store area-weighted mass (Gt)
    STPC_p_results(p).MASS.unit='Gt';
    STPC_p_results(p).MASS.fillM_Slepian=STPC_p_result.MASS.fillM_Slepian;
    STPC_p_results(p).MASS.fillM_MSSA=STPC_p_result.MASS.fillM_MSSA;
    STPC_p_results(p).MASS.fillM_STPC=STPC_p_result.MASS.fillM_STPC;

    if strcmp(land_or_ocean,'ocean')

        % store area-weighted mass-term sea level (cm)
        STPC_p_results(p).MSL.unit='cm';
        STPC_p_results(p).MSL.fillM_Slepian=STPC_p_result.MSL.fillM_Slepian;
        STPC_p_results(p).MSL.fillM_MSSA=STPC_p_result.MSL.fillM_MSSA;
        STPC_p_results(p).MSL.fillM_STPC=STPC_p_result.MSL.fillM_STPC;
        STPC_p_results(p).MSL.fillM_Grid_Slepian=STPC_p_result.MSL.fillM_Grid_Slepian;
        STPC_p_results(p).MSL.fillM_Grid_MSSA=STPC_p_result.MSL.fillM_Grid_MSSA;
        STPC_p_results(p).MSL.fillM_Grid_STPC=STPC_p_result.MSL.fillM_Grid_STPC;
        
        % store inverted barometer correction for OBP (cm)
        STPC_p_results(p).MSL.fillM_Slepian_IB=STPC_p_result.MSL.fillM_Slepian_IB;
        STPC_p_results(p).MSL.fillM_MSSA_IB=STPC_p_result.MSL.fillM_MSSA_IB;
        STPC_p_results(p).MSL.fillM_STPC_IB=STPC_p_result.MSL.fillM_STPC_IB;

    end

    STPC_p_results(p).selected_S_index = p_S(p);
    STPC_p_results(p).selected_N_index = p_N(p);
    STPC_p_results(p).selected_S = ...
        turn_V_four(p_S(p), used_sta_SN);
    selectedSWithOcean = STPC_p_results(p).selected_S + ...
        double(strcmp(land_or_ocean, 'ocean'));
    selectedN = zeros(1, selectedSWithOcean);
    for ss = 1:selectedSWithOcean
        selectedN(ss) = turn_MSSA_four{ss}(p_N(p), used_sta_SN);
    end
    STPC_p_results(p).selected_N = selectedN;
    STPC_p_results(p).selection = struct( ...
        'S', S_selection, 'N', N_selection);
    STPC_p_results(p).screening_RMS = Square_Seq{p};
    STPC_p_results(p).timing = local_pack_timing( ...
        candidateTiming, finalTiming, sharedDecompositionSeconds, ...
        turningPointSeconds, p);

end

save(datasetsFile, 'STPC_p_results', 'p_S', 'p_N', 'used_sta_SN')

end

function timing = local_zero_timing()
timing = struct( ...
    'mssa_slice_and_significance_seconds', 0, ...
    'coefficient_reconstruction_seconds', 0, ...
    'grid_reconstruction_seconds', 0, ...
    'screen_rms_analysis_seconds', 0, ...
    'signal_residual_grid_seconds', 0, ...
    'full_rms_analysis_seconds', 0);
end

function timing = local_pack_timing(screen, final, sharedSeconds, ...
    turningPointSeconds, pIndex)
timing = struct();
timing.turning_point_seconds = turningPointSeconds;
timing.turning_point_charged_to_this_p = ...
    double(pIndex == 1) * turningPointSeconds;
timing.shared_decomposition_cache_seconds = sharedSeconds;
timing.shared_decomposition_charged_to_this_p = double(pIndex == 1) * sharedSeconds;
timing.screen = screen;
timing.final = final;

timing.algorithm_seconds = ...
    timing.turning_point_charged_to_this_p + ...
    timing.shared_decomposition_charged_to_this_p + ...
    screen.mssa_slice_and_significance_seconds + ...
    screen.coefficient_reconstruction_seconds + ...
    screen.grid_reconstruction_seconds + ...
    screen.signal_residual_grid_seconds + ...
    final.mssa_slice_and_significance_seconds + ...
    final.coefficient_reconstruction_seconds + ...
    final.grid_reconstruction_seconds + ...
    final.signal_residual_grid_seconds;

timing.analysis_rms_seconds = ...
    screen.screen_rms_analysis_seconds + ...
    screen.full_rms_analysis_seconds + ...
    final.screen_rms_analysis_seconds + ...
    final.full_rms_analysis_seconds;
end

