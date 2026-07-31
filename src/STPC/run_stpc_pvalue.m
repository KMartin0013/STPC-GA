function [STPC_p_SN_results, RMS_SN_Results, timing] = run_stpc_pvalue(basicInfo, mssa_results, ...
    slepian_results, turn_V_four, turn_MSSA_four, Noise_SigLev, ...
    Turning_number, MSSA_evalues, used_sta_SN, Attach, options)

if nargin < 11 || isempty(options)
    options = struct();
end
mode = local_option(options, 'mode', 'legacy');
decompositionCache = local_option(options, 'decompositionCache', []);
selectedSN = local_option(options, 'selectedSN', []);
screenOnly = strcmp(mode, 'screen');
plotCandidateEigenvalues = local_option(options, ...
    'plotCandidateEigenvalues', strcmp(mode, 'legacy'));
plotSelectedDecomposition = local_option(options, ...
    'plotSelectedDecomposition', strcmp(mode, 'legacy'));

c11cmn       = basicInfo.c11cmn;
Area         = basicInfo.Area;
ddir1        = basicInfo.ddir1;
ddir3        = basicInfo.ddir3;      % output figure directory (Results_*/AddData)
redo         = basicInfo.redo;
plotProcess  = basicInfo.plotProcess;
MaxNumChanges_V     = size(turn_V_four, 1);
MaxNumChanges_MSSA  = size(turn_MSSA_four{1}, 1);

V               = slepian_results(1).V;
lon             = slepian_results(1).grid.lon;
lat             = slepian_results(1).grid.lat;
XY              = slepian_results(1).XY;
XY_ori          = slepian_results(1).XY_ori;
XY_buffer       = slepian_results(1).XY_buffer;
BasinArea       = slepian_results.BasinArea;
r_record        = slepian_results(1).grid.r_record;

intit_num           = mssa_results.intit_num;
mssa_Sort           = mssa_results.mssa_Sort;
fillM_deltacoffs    = mssa_results.fillM_deltacoffs;
M_rec               = mssa_results.M_rec;

fill_dates          = mssa_Sort(1).fill_dates;
fill_nmonths        = numel(fill_dates);

used_sta_V          = used_sta_SN; % what kind of turning point do you want? (default is RSS)sed_sta_V=4; % what kind of turning point do you want? (default is RSS)
used_sta_MSSA       = used_sta_SN; % what kind of turning point do you want? (default is RSS)

fitwhat = [3 365.25 365.25/2]; % linear trend and an annual and a semiannual harmonic term
candidateFile = fullfile(ddir3, ['MainSTPC_V3_' ...
    Attach.Attach_ALL '_p' num2str(Noise_SigLev) '.mat']);

% redo=true;
if strcmp(mode, 'legacy') && ~redo && exist(candidateFile,'file')

    disp(['STPC: Loading ', candidateFile, newline, ...
        'for results at the significance level of ' num2str(Noise_SigLev) '.']);

    load(candidateFile, 'STPC_p_SN_results')
    RMS_SN_Results = struct();
    timing = local_zero_timing();

    return

end

[lonlon,latlat] = meshgrid(lon, lat);

% if the longitude ranges from 0 to 360. Make it consistent to [-180,180]
if XY(1,1)>180
    XY(:,1)=XY(:,1)-360;
end

[in,on] = check_polygon_in(XY,lonlon,latlat); %#ok<ASGLU>

if ~sum(in,'all')
    error('Be caureful on your chosen XY while checking the inner points.')
end

% area-weighted should be used to more precisely calculate the regional
% basin-averaged variables
c11cmn_area = zeros(length(lat),length(lon));
for i = 1:length(lat)
    for j = 1:length(lon)
        c11cmn_area(i,j) = areaquad(lat(i)-0.5,lon(j)-0.5, ...
            lat(i)+0.5,lon(j)+0.5);
    end
end

% the TP(S1~Sn)
STPC_p_SN_results=cell(MaxNumChanges_V,MaxNumChanges_MSSA);
RMS_SN_Results = struct();
RMS_SN_Results.s_ss_nn_rms = nan(5, MaxNumChanges_V, MaxNumChanges_MSSA);
timing = struct( ...
    'mssa_slice_and_significance_seconds', 0, ...
    'coefficient_reconstruction_seconds', 0, ...
    'grid_reconstruction_seconds', 0, ...
    'screen_rms_analysis_seconds', 0, ...
    'signal_residual_grid_seconds', 0, ...
    'full_rms_analysis_seconds', 0);

if isempty(selectedSN)
    useTPss = 1:MaxNumChanges_V;
    useTPnn = 1:MaxNumChanges_MSSA;
else
    useTPss = selectedSN(1);
    useTPnn = selectedSN(2);
end

for TP_ss=useTPss
% for TP_ss=MaxNumChanges_V
    %       for  TP_ss=5;
    %           for  TP_ss=S_shannon

    % the TP(N1~Nn)
    for TP_nn=useTPnn
%     for TP_nn=MaxNumChanges_MSSA

        % for slepian
        S_rec=turn_V_four(TP_ss,used_sta_V);
        % only for shannon
        %         S_rec=TP_ss;

        N_rec=[];
        for k=1:S_rec
            N_rec(k)=turn_MSSA_four{k}(TP_nn,used_sta_MSSA);
        end

        if strcmp(mssa_Sort(1).process{3},'ocean')
            S_ol=S_rec+1;
            N_rec(S_rec+1)=turn_MSSA_four{S_rec+1}(TP_nn,used_sta_MSSA);
        else
            S_ol=S_rec;
        end

        %% decompose each coefficient (M-SSA) following (Gauer et al., 2022)

        stageTimer = tic;
        [MSSAcoffs_rec, MSSAcoffs_evalues, MSSAcoffs_RC, MSSAcoffs_RCTest] = ...
            run_MSSA_decompose(mssa_Sort, fillM_deltacoffs, M_rec, N_rec, ...
            Noise_SigLev, fill_nmonths, intit_num, S_ol, S_rec, ...
            decompositionCache);
        timing.mssa_slice_and_significance_seconds = ...
            timing.mssa_slice_and_significance_seconds + toc(stageTimer);

        %% coefficient reconstruction

        stageTimer = tic;
        mssa_Sort_coff_sn = coffs_reconstruction(MSSAcoffs_RCTest, ...
            MSSAcoffs_RC, fillM_deltacoffs, MSSAcoffs_rec, ...
            MSSAcoffs_evalues, S_ol, fill_nmonths, intit_num);
        timing.coefficient_reconstruction_seconds = ...
            timing.coefficient_reconstruction_seconds + toc(stageTimer);

        if plotProcess && plotCandidateEigenvalues
            plot_mssa_eigenvalues(mssa_Sort_coff_sn,intit_num,S_rec,V, ...
                TP_ss,TP_nn,Turning_number,M_rec,Noise_SigLev,Attach)
        end

        %% grid reconstruction

        stageTimer = tic;
        [mssa_Sort_EWHorMSL_sn,c11cmn_area,in,on] = grids_reconstruction_EWHorMSL( ...
            mssa_Sort_coff_sn,intit_num,fill_nmonths,S_rec,S_ol, ...
            r_record, lat, lon, XY, BasinArea);
        timing.grid_reconstruction_seconds = ...
            timing.grid_reconstruction_seconds + toc(stageTimer);

        if screenOnly
            % V3 lightweight screening retains only the exact quantity
            % used by run_stpc_main to select S and N.
            analysisTimer = tic;
            meanResult = mssa_Sort_EWHorMSL_sn(intit_num+1);
            meanGrid = meanResult.fillM_Grid_EWH_MSSA_both;
            gridRms = squeeze(rms(meanGrid, 1));
            RMS_SN_Results.s_ss_nn_rms(3, TP_ss, TP_nn) = ...
                sum(gridRms(in) .* c11cmn_area(in)) / ...
                sum(c11cmn_area(in));
            timing.screen_rms_analysis_seconds = ...
                timing.screen_rms_analysis_seconds + toc(analysisTimer);
            continue
        end

%         mssa_Sort_RMS_sn = grids_reconstruction_EWHorMSL_RMS( ...
%             mssa_Sort_coff_sn, intit_num, fill_nmonths, S_rec,...
%             r_record, lonlon);

        [mssa_Sort_EWHorMSL_sigres_sn,mssa_Sort_EWHorMSL_RMS_sigres_sn, ...
            productTiming] = ...
            grids_reconstruction_EWHorMSL_RMS_sigres( ...
            mssa_Sort_coff_sn,intit_num,fill_nmonths,S_rec,r_record, ...
            c11cmn_area,in,BasinArea,S_ol,fill_dates,fitwhat,lonlon);
        timing.signal_residual_grid_seconds = ...
            timing.signal_residual_grid_seconds + ...
            productTiming.grid_products_seconds;
        timing.full_rms_analysis_seconds = ...
            timing.full_rms_analysis_seconds + ...
            productTiming.rms_analysis_seconds;

        %% calculate the explained eigenvalues for each coefficient for each components (e.g., annual)
        %         load(fullfile(getenv('IFILES'),['Results_' parrent_path],['Main_' Attach_ALL '_pro.mat']));

        if plotProcess && ~plotCandidateEigenvalues && strcmp(mode, 'legacy')
            plot_mssa_eigenvalues(mssa_Sort_coff_sn,intit_num,S_rec,V, ...
                TP_ss,TP_nn,Turning_number,M_rec,Noise_SigLev,Attach)
        end

        %% select specific SSF, coefficient and some RCs (also show the SSF basis)
        if plotProcess && plotSelectedDecomposition && ...
                (~strcmp(mode, 'legacy') || ...
                (TP_ss==MaxNumChanges_V && TP_nn==MaxNumChanges_MSSA))
            plot_mssa_decompose(mssa_Sort, mssa_Sort_coff_sn, ...
                r_record, lonlon, c11cmn, Area, XY_ori, XY_buffer, ...
                V, MSSAcoffs_evalues, S_rec, S_ol, intit_num, Attach);
        end

        %% only save necessary output for (S, N) selection

        STPC_p_result=struct();

        % we choose the mean of the four institutions
        STPC_p_result.name=mssa_Sort(intit_num+1).name;
        STPC_p_result.c11cmn=c11cmn;
        STPC_p_result.lon=lon;STPC_p_result.lat=lat;

        STPC_p_result.data_year_beg=mssa_Sort(intit_num+1).data_year_beg;
        STPC_p_result.process=mssa_Sort(intit_num+1).process;
        STPC_p_result.data_dates=mssa_Sort(intit_num+1).data_dates;
        STPC_p_result.fill_dates=mssa_Sort(intit_num+1).fill_dates;
        STPC_p_result.fill_months=mssa_Sort(intit_num+1).fill_months;
        STPC_p_result.use_dates=mssa_Sort(intit_num+1).use_dates;
        STPC_p_result.use_months=mssa_Sort(intit_num+1).use_months;
        STPC_p_result.missing_dates=mssa_Sort(intit_num+1).missing_dates;
        STPC_p_result.missing_months=mssa_Sort(intit_num+1).missing_months;

        % eigenvalues for signal and noise
        STPC_p_result.fillM_reconcoffs_eva=mssa_Sort_coff_sn(intit_num+1).fillM_reconcoffs_eva;
        STPC_p_result.fillM_reconcoffs_RCTest=mssa_Sort_coff_sn(intit_num+1).fillM_reconcoffs_RCTest;

        %%%%% EWH %%%%%
        % Time series
        EWH=struct();
        EWH.fillM_Slepian=mssa_Sort_EWHorMSL_sn(intit_num+1).fillM_EWH;
        EWH.fillM_MSSA=mssa_Sort_EWHorMSL_sn(intit_num+1).fillM_EWH_MSSA;
        EWH.fillM_STPC=mssa_Sort_EWHorMSL_sn(intit_num+1).fillM_EWH_MSSA_both;

        EWH.fillM_sig_STPC=mssa_Sort_EWHorMSL_sigres_sn(intit_num+1).fillM_EWHsig_MSSA_both;
        EWH.fillM_res_STPC=mssa_Sort_EWHorMSL_sigres_sn(intit_num+1).fillM_EWHres_MSSA_both;

        % Grid
        EWH.fillM_Grid_Slepian=mssa_Sort_EWHorMSL_sn(intit_num+1).fillM_Grid_EWH;
        EWH.fillM_Grid_MSSA=mssa_Sort_EWHorMSL_sn(intit_num+1).fillM_Grid_EWH_MSSA;
        EWH.fillM_Grid_STPC=mssa_Sort_EWHorMSL_sn(intit_num+1).fillM_Grid_EWH_MSSA_both;

        EWH.fillM_Grid_sig_STPC=mssa_Sort_EWHorMSL_sigres_sn(intit_num+1).fillM_Grid_EWHsig_MSSA_both;
        EWH.fillM_Grid_res_STPC=mssa_Sort_EWHorMSL_sigres_sn(intit_num+1).fillM_Grid_EWHres_MSSA_both;

        % RMS
        %         EWH.fillM_Grid_MSSA_uptoS_RMS=mssa_Sort_EWHorMSL_RMS_sigres_sn.fillM_Grid_EWH_MSSA_uptoS_RMS;
        %         EWH.fillM_Grid_STPC_uptoS_RMS=mssa_Sort_EWHorMSL_RMS_sigres_sn.fillM_Grid_EWH_MSSA_both_uptoS_RMS;

%         EWH.fillM_Grid_MSSA_bothfail_uptoS_RMS=mssa_Sort_EWHorMSL_RMS_sigres_sn.fillM_Grid_EWH_MSSA_bothfail_uptoS_RMS;
%         %This is a bug
        EWH.fillM_Grid_MSSA_bothfail_uptoS_RMS=mssa_Sort_EWHorMSL_RMS_sigres_sn(intit_num+1).fillM_Grid_EWH_MSSA_bothfail_uptoS_RMS;
        % Bug fixed (2026/4/14)

        %%%%% Mass %%%%%
        % Time series
        MASS.fillM_Slepian=mssa_Sort_EWHorMSL_sn(intit_num+1).fillM_MASS;
        MASS.fillM_MSSA=mssa_Sort_EWHorMSL_sn(intit_num+1).fillM_MASS_MSSA;
        MASS.fillM_STPC=mssa_Sort_EWHorMSL_sn(intit_num+1).fillM_MASS_MSSA_both;

        MASS.fillM_sig_STPC=mssa_Sort_EWHorMSL_sigres_sn(intit_num+1).fillM_MASSsig_MSSA_both;
        MASS.fillM_res_STPC=mssa_Sort_EWHorMSL_sigres_sn(intit_num+1).fillM_MASSres_MSSA_both;

        STPC_p_result.EWH=EWH;
        STPC_p_result.MASS=MASS;

        if strcmp(mssa_Sort(1).process{3},'ocean')

            %%%%% MSL %%%%%
            % Time series
            MSL.fillM_Slepian=mssa_Sort_EWHorMSL_sn(intit_num+1).fillM_MSL;
            MSL.fillM_MSSA=mssa_Sort_EWHorMSL_sn(intit_num+1).fillM_MSL_MSSA;
            MSL.fillM_STPC=mssa_Sort_EWHorMSL_sn(intit_num+1).fillM_MSL_MSSA_both;

            MSL.fillM_sig_STPC=mssa_Sort_EWHorMSL_sigres_sn(intit_num+1).fillM_MSLsig_MSSA_both;
            MSL.fillM_res_STPC=mssa_Sort_EWHorMSL_sigres_sn(intit_num+1).fillM_MSLres_MSSA_both;

            % Grid
            MSL.fillM_Grid_Slepian=mssa_Sort_EWHorMSL_sn(intit_num+1).fillM_Grid_MSL;
            MSL.fillM_Grid_MSSA=mssa_Sort_EWHorMSL_sn(intit_num+1).fillM_Grid_MSL_MSSA;
            MSL.fillM_Grid_STPC=mssa_Sort_EWHorMSL_sn(intit_num+1).fillM_Grid_MSL_MSSA_both;

            MSL.fillM_Grid_sig_STPC=mssa_Sort_EWHorMSL_sigres_sn(intit_num+1).fillM_Grid_MSLsig_MSSA_both;
            MSL.fillM_Grid_res_STPC=mssa_Sort_EWHorMSL_sigres_sn(intit_num+1).fillM_Grid_MSLres_MSSA_both;

            % IB for OBP calculation
            MSL.fillM_Slepian_IB=mssa_Sort_coff_sn(intit_num+1).fillM_deltacoffs(:,S_ol) ...
                * 1000/1028/10 ; %IB: mmH2O (kg/m^2) -> cm
            MSL.fillM_MSSA_IB=mssa_Sort_coff_sn(intit_num+1).fillM_reconcoffs(:,S_ol) ...
                * 1000/1028/10 ; %IB: mmH2O (kg/m^2) -> cm
            MSL.fillM_STPC_IB=mssa_Sort_coff_sn(intit_num+1).fillM_reconcoffs_both(:,S_ol) ...
                * 1000/1028/10 ; %IB: mmH2O (kg/m^2) -> cm

            STPC_p_result.MSL=MSL;
        end

        STPC_p_SN_results{TP_ss,TP_nn}=STPC_p_result;

%         STPC_results_struct=STPC_result;
    end

end

if strcmp(mode, 'legacy')
    save(candidateFile, 'STPC_p_SN_results')
end

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

function value = local_option(options, name, defaultValue)
if isfield(options, name)
    value = options.(name);
else
    value = defaultValue;
end
end
