function [mssa_Sort_sigres_sn,mssa_Sort_RMS_sigres_sn,timing] = grids_reconstruction_EWHorMSL_RMS_sigres( ...
    mssa_Sort_sn,intit_num,fill_nmonths,S_rec,r_record, ...
    c11cmn_area,in,BasinArea,S_ol,fill_dates,fitwhat,lonlon)
% -------------------------------------------------------------------------
% compute_MSSA_signal_residual_all
%
% Collect all signal/residual-related computations (previously commented)
% into a single function, including:
%   1) Decomposition of MSSA reconstruction coefficients into signal
%      (fit) and residual parts using slept2resid_m.
%   2) Reconstruction of gridded EWH signal/residual and MASS signal/
%      residual (both MSSA and MSSA_both).
%   3) Basin-averaged (regional) EWH/MASS/MSL signal and residual.
%   4) Gridded MSL signal/residual from EWH signal/residual.
%   5) RMS statistics of:
%        - original EWH, signal, residual
%        - upto S_sum EWH (original / signal / residual)
%        - MSSA and MSSA_both (original & upto S_sum)
%        - MSSA_both signal/residual (original & upto S_sum)
%        - MSSA residuals and both-fail terms (original & upto S_sum)
%
% INPUT:
%   mssa_Sort_sn   : struct array with MSSA-related fields
%   intit_num      : number of MSSA iterations, loop 1:intit_num+1
%   fill_nmonths   : number of time steps (months)
%   S_rec          : number of Slepian functions used for reconstruction
%   r_record       : Slepian basis on grid, size [S_rec, nlat, nlon]
%   c11cmn_area    : grid-cell area matrix (nlat x nlon)
%   in             : logical mask (same size as c11cmn_area) for basin area
%   BasinArea      : basin area in m^2
%   S_ol           : index of IB coefficient in Slepian expansion
%   fill_dates     : time vector for fill_nmonths (e.g., datenum / days)
%   fitwhat        : vector specifying fit components [trend, annual, ...]
%   lonlon,latlat  : 2-D lon/lat grids defining spatial resolution
%
% OUTPUT:
%   mssa_Sort_sn     : updated struct with signal/residual fields:
%                      - fillM_uptoS_EWHsig_MSSA, fillM_uptoS_EWHres_MSSA, ...
%                      - fillM_Grid_EWHsig_MSSA, fillM_Grid_EWHres_MSSA, ...
%                      - fillM_Grid_MASSsig_MSSA, fillM_Grid_MASSres_MSSA, ...
%                      - fillM_EWHsig_MSSA, fillM_EWHres_MSSA, ...
%                      - fillM_MASSsig_MSSA, fillM_MASSres_MSSA, ...
%                      - fillM_MSLsig_MSSA, fillM_MSLres_MSSA, ...
%                      and their *_both variants.
%
%   mssa_Sort_RMS_sn : struct with RMS metrics for:
%                      - original, signal, residual EWH (grid)
%                      - upto S_sum (original, signal, residual)
%                      - MSSA, MSSA_both, MSSA residuals, both-fail
%                      - MSSA_both signal/residual
%                      and their upto-S variants.
%
% NOTE:
%   This function assumes the presence of fields in mssa_Sort_sn:
%     fillM_deltacoffs, fillM_reconcoffs, fillM_reconcoffs_both,
%     fillM_IBdeltacoffs, fill_signaldelta,
%     fill_MSSAcoffs, fill_MSSAcoffs_both,
%     fill_MSSA_bothsignaldelta, fill_MSSA_bothresiddelta
%   and an external function slept2resid_m.
% -------------------------------------------------------------------------

%% 1. Signal / residual decomposition and reconstruction (EWH / MASS / MSL)

mssa_Sort_sigres_sn = struct();
gridTimer = tic;
for ins = 1:intit_num+1

    % Original MSSA reconstruction coefficients
    fillM_reconcoffs_ins      = mssa_Sort_sn(ins).fillM_reconcoffs;
    fillM_reconcoffs_both_ins = mssa_Sort_sn(ins).fillM_reconcoffs_both;

    % Decompose MSSA coefficients into signal (fit) and residual (unfit)
    % Becasue we only use this to fit, so so we have temporarily removed 
    % the warnings to check the number of coefficients (in addmoff.m)
    warning off;
    [fillM_reconsignal_ins,fillM_reconresid_ins,~,~] = slept2resid_m( ...
        fillM_reconcoffs_ins, fill_dates, fitwhat, [], [], [], [], [], []); % mm (kg/m^2 equiv)

    [fillM_reconsignal_both_ins,fillM_reconresid_both_ins,~,~] = slept2resid_m( ...
        fillM_reconcoffs_both_ins, fill_dates, fitwhat, [], [], [], [], [], []); % mm (kg/m^2 equiv)
    warning on

    % Loop over time to reconstruct gridded signal / residual
    for i = 1:fill_nmonths
        sp_ewhsig_mssa_reconst      = 0;
        sp_ewhres_mssa_reconst      = 0;
        sp_ewhsig_mssa_both_reconst = 0;
        sp_ewhres_mssa_both_reconst = 0;

        for j = 1:S_rec
            r = squeeze(r_record(j,:,:));

            % signal part (MSSA)
            sp_ewhsig_mssa_reconst = sp_ewhsig_mssa_reconst + ...
                r * fillM_reconsignal_ins(i,j)' / 1000 * 100; % kg/m^2 -> cm

            % residual part (MSSA)
            sp_ewhres_mssa_reconst = sp_ewhres_mssa_reconst + ...
                r * fillM_reconresid_ins(i,j)' / 1000 * 100; % kg/m^2 -> cm

            % signal part (MSSA_both)
            sp_ewhsig_mssa_both_reconst = sp_ewhsig_mssa_both_reconst + ...
                r * fillM_reconsignal_both_ins(i,j)' / 1000 * 100; % kg/m^2 -> cm

            % residual part (MSSA_both)
            sp_ewhres_mssa_both_reconst = sp_ewhres_mssa_both_reconst + ...
                r * fillM_reconresid_both_ins(i,j)' / 1000 * 100; % kg/m^2 -> cm

            % regional (upto S) EWH signal/residual, area-weighted
            mssa_Sort_sigres_sn(ins).fillM_uptoS_EWHsig_MSSA(j,i)      = sum(sp_ewhsig_mssa_reconst(in)      .* c11cmn_area(in)) / sum(c11cmn_area(in));
            mssa_Sort_sigres_sn(ins).fillM_uptoS_EWHres_MSSA(j,i)      = sum(sp_ewhres_mssa_reconst(in)      .* c11cmn_area(in)) / sum(c11cmn_area(in));
            mssa_Sort_sigres_sn(ins).fillM_uptoS_EWHsig_MSSA_both(j,i) = sum(sp_ewhsig_mssa_both_reconst(in) .* c11cmn_area(in)) / sum(c11cmn_area(in));
            mssa_Sort_sigres_sn(ins).fillM_uptoS_EWHres_MSSA_both(j,i) = sum(sp_ewhres_mssa_both_reconst(in) .* c11cmn_area(in)) / sum(c11cmn_area(in));
        end

        % gridded EWH signal / residual (cm)
        mssa_Sort_sigres_sn(ins).fillM_Grid_EWHsig_MSSA(i,:,:)       = sp_ewhsig_mssa_reconst;
        mssa_Sort_sigres_sn(ins).fillM_Grid_EWHres_MSSA(i,:,:)       = sp_ewhres_mssa_reconst;
        mssa_Sort_sigres_sn(ins).fillM_Grid_EWHsig_MSSA_both(i,:,:)  = sp_ewhsig_mssa_both_reconst;
        mssa_Sort_sigres_sn(ins).fillM_Grid_EWHres_MSSA_both(i,:,:)  = sp_ewhres_mssa_both_reconst;

        % gridded MASS signal / residual (Gt)
        mssa_Sort_sigres_sn(ins).fillM_Grid_MASSsig_MSSA(i,:,:)      = sp_ewhsig_mssa_reconst      /100 .* c11cmn_area * 4*pi*6370000^2 * 10^3 / 10^3 / 10^9;
        mssa_Sort_sigres_sn(ins).fillM_Grid_MASSres_MSSA(i,:,:)      = sp_ewhres_mssa_reconst      /100 .* c11cmn_area * 4*pi*6370000^2 * 10^3 / 10^3 / 10^9;
        mssa_Sort_sigres_sn(ins).fillM_Grid_MASSsig_MSSA_both(i,:,:) = sp_ewhsig_mssa_both_reconst /100 .* c11cmn_area * 4*pi*6370000^2 * 10^3 / 10^3 / 10^9;
        mssa_Sort_sigres_sn(ins).fillM_Grid_MASSres_MSSA_both(i,:,:) = sp_ewhres_mssa_both_reconst /100 .* c11cmn_area * 4*pi*6370000^2 * 10^3 / 10^3 / 10^9;

        % basin-averaged EWH signal / residual (cm)
        mssa_Sort_sigres_sn(ins).fillM_EWHsig_MSSA(i)      = sum(sp_ewhsig_mssa_reconst(in)      .* c11cmn_area(in)) / sum(c11cmn_area(in));
        mssa_Sort_sigres_sn(ins).fillM_EWHres_MSSA(i)      = sum(sp_ewhres_mssa_reconst(in)      .* c11cmn_area(in)) / sum(c11cmn_area(in));
        mssa_Sort_sigres_sn(ins).fillM_EWHsig_MSSA_both(i) = sum(sp_ewhsig_mssa_both_reconst(in) .* c11cmn_area(in)) / sum(c11cmn_area(in));
        mssa_Sort_sigres_sn(ins).fillM_EWHres_MSSA_both(i) = sum(sp_ewhres_mssa_both_reconst(in) .* c11cmn_area(in)) / sum(c11cmn_area(in));

        % basin-averaged MASS signal / residual (Gt)
        mssa_Sort_sigres_sn(ins).fillM_MASSsig_MSSA(i)      = mssa_Sort_sigres_sn(ins).fillM_EWHsig_MSSA(i)      /100 * BasinArea * 10^3 / 10^3 / 10^9;
        mssa_Sort_sigres_sn(ins).fillM_MASSres_MSSA(i)      = mssa_Sort_sigres_sn(ins).fillM_EWHres_MSSA(i)      /100 * BasinArea * 10^3 / 10^3 / 10^9;
        mssa_Sort_sigres_sn(ins).fillM_MASSsig_MSSA_both(i) = mssa_Sort_sigres_sn(ins).fillM_EWHsig_MSSA_both(i) /100 * BasinArea * 10^3 / 10^3 / 10^9;
        mssa_Sort_sigres_sn(ins).fillM_MASSres_MSSA_both(i) = mssa_Sort_sigres_sn(ins).fillM_EWHres_MSSA_both(i) /100 * BasinArea * 10^3 / 10^3 / 10^9;
    end

    % If IB index exists, compute MSL signal / residual
    if S_ol > S_rec
        for i = 1:fill_nmonths

            % basin-averaged MSL signal / residual (cm), from MASS
            mssa_Sort_sigres_sn(ins).fillM_MSLsig_MSSA(i) = ...
                mssa_Sort_sigres_sn(ins).fillM_MASSsig_MSSA(i) * 10^3 * 10^9 / 10^3 / BasinArea * 100 - ...
                fillM_reconsignal_ins(i,S_ol) * 1000/1028/10;

            mssa_Sort_sigres_sn(ins).fillM_MSLres_MSSA(i) = ...
                mssa_Sort_sigres_sn(ins).fillM_MASSres_MSSA(i) * 10^3 * 10^9 / 10^3 / BasinArea * 100 - ...
                fillM_reconresid_ins(i,S_ol) * 1000/1028/10;

            mssa_Sort_sigres_sn(ins).fillM_MSLsig_MSSA_both(i) = ...
                mssa_Sort_sigres_sn(ins).fillM_MASSsig_MSSA_both(i) * 10^3 * 10^9 / 10^3 / BasinArea * 100 - ...
                fillM_reconsignal_both_ins(i,S_ol) * 1000/1028/10;

            mssa_Sort_sigres_sn(ins).fillM_MSLres_MSSA_both(i) = ...
                mssa_Sort_sigres_sn(ins).fillM_MASSres_MSSA_both(i) * 10^3 * 10^9 / 10^3 / BasinArea * 100 - ...
                fillM_reconresid_both_ins(i,S_ol) * 1000/1028/10;

            % gridded MSL signal / residual (cm), only *_both versions in original
            mssa_Sort_sigres_sn(ins).fillM_Grid_MSLsig_MSSA_both(i,:,:) = ...
                mssa_Sort_sigres_sn(ins).fillM_Grid_EWHsig_MSSA_both(i,:,:) * 1000/1028 - ...
                fillM_reconsignal_both_ins(i,S_ol) * 1000/1028/10;

            mssa_Sort_sigres_sn(ins).fillM_Grid_MSLres_MSSA_both(i,:,:) = ...
                mssa_Sort_sigres_sn(ins).fillM_Grid_EWHres_MSSA_both(i,:,:) * 1000/1028 - ...
                fillM_reconresid_both_ins(i,S_ol) * 1000/1028/10;
        end
    end
end

%% 2. RMS: include all signal / residual related RMS (originally commented)

mssa_Sort_RMS_sigres_sn = struct();
timing = struct();
timing.grid_products_seconds = toc(gridTimer);
rmsTimer = tic;

for ins = 1:intit_num+1

    % original coefficients
    fillM_deltacoffs_ins      = mssa_Sort_sn(ins).fillM_deltacoffs;

    % Decompose MSSA coefficients into signal (fit) and residual (unfit)
    % Becasue we only use this to fit, so so we have temporarily removed 
    % the warnings to check the number of coefficients (in addmoff.m)
    warning off;
    [fillM_signaldelta_ins,~,~,~] = slept2resid_m( ...
        fillM_deltacoffs_ins, fill_dates, fitwhat, [], [], [], [], [], []); % mm (kg/m^2 equiv)

    % we think that the signal components of coefficients with
    % 8-coefficient interpolation and MSSA interpolation should be very
    % similar. SO fillM_signaldelta_ins ~ fill_signaldelta_ins.

    % Original MSSA reconstruction coefficients
    fillM_reconcoffs_ins      = mssa_Sort_sn(ins).fillM_reconcoffs;
    fillM_reconcoffs_both_ins = mssa_Sort_sn(ins).fillM_reconcoffs_both;

    % Decompose MSSA coefficients into signal (fit) and residual (unfit)
    [fillM_reconsignal_ins,fillM_reconresid_ins,~,~] = slept2resid_m( ...
        fillM_reconcoffs_ins, fill_dates, fitwhat, [], [], [], [], [], []); % mm (kg/m^2 equiv)

    [fillM_reconsignal_both_ins,fillM_reconresid_both_ins,~,~] = slept2resid_m( ...
        fillM_reconcoffs_both_ins, fill_dates, fitwhat, [], [], [], [], [], []); % mm (kg/m^2 equiv)
    warning on;

    % Cumulative prefix fields must persist across j so that element j
    % represents the RMS of the first j Slepian modes, not mode j alone.
    sp_ewh_uptoS        = zeros(fill_nmonths,size(lonlon,1),size(lonlon,2));
    sp_ewh_uptoS_signal = zeros(fill_nmonths,size(lonlon,1),size(lonlon,2));
    sp_ewh_uptoS_mssa   = zeros(fill_nmonths,size(lonlon,1),size(lonlon,2));
    sp_ewh_uptoS_mssa_both = zeros(fill_nmonths,size(lonlon,1),size(lonlon,2));
    sp_ewh_uptoS_mssa_both_signal = zeros(fill_nmonths,size(lonlon,1),size(lonlon,2));
    sp_ewh_uptoS_mssa_both_resid  = zeros(fill_nmonths,size(lonlon,1),size(lonlon,2));

    for j = 1:S_rec
        % Per-mode fields must be reset for every j.
        % time x lat x lon
        sp_ewh          = zeros(fill_nmonths,size(lonlon,1),size(lonlon,2));
        sp_ewh_signal   = zeros(fill_nmonths,size(lonlon,1),size(lonlon,2));
        sp_ewh_mssa_signal      = zeros(fill_nmonths,size(lonlon,1),size(lonlon,2));
        sp_ewh_mssa_resid       = zeros(fill_nmonths,size(lonlon,1),size(lonlon,2));
        sp_ewh_mssa_both_signal = zeros(fill_nmonths,size(lonlon,1),size(lonlon,2));
        sp_ewh_mssa_both_resid  = zeros(fill_nmonths,size(lonlon,1),size(lonlon,2));

        r = squeeze(r_record(j,:,:));

        for i = 1:fill_nmonths
            % original EWH (cm) % Bug 1 fixed 2026/4/14
            sp_ewh(i,:,:) = r * fillM_deltacoffs_ins(i,j)' / 1000 * 100; % Bug 1 fixed 2026/4/14

            % original signal EWH (cm) % Bug 1 fixed 2026/4/14
            sp_ewh_signal(i,:,:) = r * fillM_signaldelta_ins(i,j)' / 1000 * 100; % Bug 1 fixed: 2026/4/14

            % upto-S original and signal
            sp_ewh_uptoS(i,:,:)        = squeeze(sp_ewh_uptoS(i,:,:))       ...
                + r * fillM_deltacoffs_ins(i,j)' / 1000 * 100;
            sp_ewh_uptoS_signal(i,:,:) = squeeze(sp_ewh_uptoS_signal(i,:,:))...
                + r * fillM_signaldelta_ins(i,j)' / 1000 * 100;

            % MSSA reconstruction
            sp_ewh_uptoS_mssa(i,:,:) = squeeze(sp_ewh_uptoS_mssa(i,:,:))    ...
                + r * fillM_reconcoffs_ins(i,j)' / 1000 * 100;

            % MSSA_both (rejected) reconstruction
            sp_ewh_uptoS_mssa_both(i,:,:) = squeeze(sp_ewh_uptoS_mssa_both(i,:,:)) ...
                + r * fillM_reconcoffs_both_ins(i,j)' / 1000 * 100;

            % MSSA_both signal / residual components
            sp_ewh_mssa_signal(i,:,:) = r * fillM_reconsignal_ins(i,j)' / 1000 * 100;
            sp_ewh_mssa_resid(i,:,:)  = r * fillM_reconresid_ins(i,j)'  / 1000 * 100;

            % MSSA_both signal / residual components
            sp_ewh_mssa_both_signal(i,:,:) = r * fillM_reconsignal_both_ins(i,j)' / 1000 * 100;
            sp_ewh_mssa_both_resid(i,:,:)  = r * fillM_reconresid_both_ins(i,j)'  / 1000 * 100;

            sp_ewh_uptoS_mssa_both_signal(i,:,:) = squeeze(sp_ewh_uptoS_mssa_both_signal(i,:,:)) + ...
                r * fillM_reconsignal_both_ins(i,j)' / 1000 * 100;

            sp_ewh_uptoS_mssa_both_resid(i,:,:)  = squeeze(sp_ewh_uptoS_mssa_both_resid(i,:,:))  + ...
                r * fillM_reconresid_both_ins(i,j)'  / 1000 * 100;
        end

        % compute RMS for each grid cell
        for p = 1:size(lonlon,1)
            for q = 1:size(lonlon,2)

                % original EWH / signal / residual RMS
                mssa_Sort_RMS_sigres_sn(ins).fillM_Grid_EWHsig_RMS(j,p,q)   = rms(squeeze(sp_ewh_signal(:,p,q)));
                mssa_Sort_RMS_sigres_sn(ins).fillM_Grid_EWHres_RMS(j,p,q)   = rms(squeeze(sp_ewh(:,p,q)) - squeeze(sp_ewh_signal(:,p,q)));

                % upto-S original / signal / residual RMS
                mssa_Sort_RMS_sigres_sn(ins).fillM_Grid_EWH_uptoS_RMS(j,p,q)        = rms(squeeze(sp_ewh_uptoS(:,p,q)));
                mssa_Sort_RMS_sigres_sn(ins).fillM_Grid_EWHsig_uptoS_RMS(j,p,q)     = rms(squeeze(sp_ewh_uptoS_signal(:,p,q)));
                mssa_Sort_RMS_sigres_sn(ins).fillM_Grid_EWHres_uptoS_RMS(j,p,q)     = rms(squeeze(sp_ewh_uptoS(:,p,q)) - squeeze(sp_ewh_uptoS_signal(:,p,q)));

                % MSSA and MSSA upto-S RMS
                mssa_Sort_RMS_sigres_sn(ins).fillM_Grid_EWH_MSSA_uptoS_RMS(j,p,q)   = rms(squeeze(sp_ewh_uptoS_mssa(:,p,q)));

                % MSSA_both and upto-S RMS
                mssa_Sort_RMS_sigres_sn(ins).fillM_Grid_EWH_MSSA_both_uptoS_RMS(j,p,q) = rms(squeeze(sp_ewh_uptoS_mssa_both(:,p,q)));

                % MSSA_both signal/residual RMS
                mssa_Sort_RMS_sigres_sn(ins).fillM_Grid_EWHsig_MSSA_RMS(j,p,q)      = rms(squeeze(sp_ewh_mssa_signal(:,p,q)));
                mssa_Sort_RMS_sigres_sn(ins).fillM_Grid_EWHres_MSSA_RMS(j,p,q)      = rms(squeeze(sp_ewh_mssa_resid(:,p,q)));
               
                % MSSA_both signal/residual RMS
                mssa_Sort_RMS_sigres_sn(ins).fillM_Grid_EWHsig_MSSA_both_RMS(j,p,q)      = rms(squeeze(sp_ewh_mssa_both_signal(:,p,q)));
                mssa_Sort_RMS_sigres_sn(ins).fillM_Grid_EWHsig_MSSA_both_uptoS_RMS(j,p,q)= rms(squeeze(sp_ewh_uptoS_mssa_both_signal(:,p,q)));
                mssa_Sort_RMS_sigres_sn(ins).fillM_Grid_EWHres_MSSA_both_RMS(j,p,q)      = rms(squeeze(sp_ewh_mssa_both_resid(:,p,q)));
                mssa_Sort_RMS_sigres_sn(ins).fillM_Grid_EWHres_MSSA_both_uptoS_RMS(j,p,q)= rms(squeeze(sp_ewh_uptoS_mssa_both_resid(:,p,q)));

                % MSSA residuals and both-fail RMS
                mssa_Sort_RMS_sigres_sn(ins).fillM_Grid_EWH_MSSA_bothfail_uptoS_RMS(j,p,q)= ...
                    rms(squeeze(sp_ewh_uptoS_mssa(:,p,q) - sp_ewh_uptoS_mssa_both(:,p,q)));
                mssa_Sort_RMS_sigres_sn(ins).fillM_Grid_EWH_MSSA_resid_uptoS_RMS(j,p,q)   = ...
                    rms(squeeze(sp_ewh_uptoS(:,p,q) - sp_ewh_uptoS_mssa(:,p,q)));
            end
        end
    end
end

timing.rms_analysis_seconds = toc(rmsTimer);
end
