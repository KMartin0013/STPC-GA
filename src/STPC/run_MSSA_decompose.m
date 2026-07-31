function [MSSAcoffs_rec, MSSAcoffs_evalues, MSSAcoffs_RC, MSSAcoffs_RCTest] ...
    = run_MSSA_decompose(mssa_Sort_sn, fillM_deltacoffs, ...
    M_rec, N_rec, Noise_SigLev, fill_nmonths, intit_num, S_ol, S_rec, ...
    decompositionCache)

% mssa_Sort_sn        % struct array with settings for each time series
% fillM_deltacoffs    % [fill_nmonths × intit_num × S_ol] filled Slepian coefficients
% M_rec               % embedding window M
% N_rec               % vector of length S_ol: N_rec(ss) for each Slepian index
% Noise_SigLev        % noise significance level
% fill_nmonths       % number of filled months
% intit_num          % number of input time series
% S_ol              % total number of Slepian indices to process
% S_rec                  % number of Slepian indices treated as standard MSSA

% run_MSSA_reconstruction
% Wrapper for the M-SSA reconstruction loop (following Gauer et al., 2022).
% For each Slepian index ss = 1:S_ol:
%   - Set mssa_Sort_sn(ins).MSSA_TS_order
%   - Call MSSA_final_noCDF_freqsort
%   - Collect reconstructed coefficients, eigenvalues, RCs and RC tests.

    if nargin < 10
        decompositionCache = [];
    end

    % Reuse the p-independent full decomposition when supplied.  The
    % significance level only changes rows 1--2 of RCTest; rows 8--9 store
    % the K-S and Lilliefors p-values and are independent of alpha.
    if ~isempty(decompositionCache)
        MSSAcoffs_evalues = decompositionCache.evalues(:, 1:S_ol);
        MSSAcoffs_rec     = nan(fill_nmonths, intit_num, S_ol);
        MSSAcoffs_RC      = cell(S_ol, 1);
        MSSAcoffs_RCTest  = cell(S_ol, 1);

        for ss = 1:S_ol
            availableN = size(decompositionCache.RC{ss}, 2);
            if N_rec(ss) > availableN
                error('V3:FixedNOutOfRange', ...
                    ['Requested N=%d for Slepian coefficient %d, but only ', ...
                     '%d reconstructed components are available.'], ...
                    N_rec(ss), ss, availableN);
            end
            useN = N_rec(ss);
            useRC = decompositionCache.RC{ss}(:, 1:useN, :);
            useTest = decompositionCache.RCTest{ss}(:, 1:useN);

            % Preserve the original test semantics for every requested p.
            useTest(1, :) = useTest(8, :) < Noise_SigLev;
            useTest(2, :) = useTest(9, :) < Noise_SigLev;

            MSSAcoffs_rec(:, :, ss) = squeeze(sum(useRC, 2));
            MSSAcoffs_RC{ss} = useRC;
            MSSAcoffs_RCTest{ss} = useTest;
        end
        return
    end

    % Preallocate outputs
    MSSAcoffs_rec     = nan(fill_nmonths, intit_num, S_ol);
    MSSAcoffs_evalues = [];              % will be expanded column-wise
    MSSAcoffs_RC      = cell(S_ol, 1);
    MSSAcoffs_RCTest  = cell(S_ol, 1);

    % Main loop over Slepian index
    for ss = 1:S_ol

        % Set MSSA_TS_order for each input time series
        if ss < S_rec + 1
            for ins = 1:intit_num
                mssa_Sort_sn(ins).MSSA_TS_order = ss;
            end
        else
            disp('we need to consider IB for ocean study.')
            for ins = 1:intit_num
                mssa_Sort_sn(ins).MSSA_TS_order = 'IB';
            end
        end

        % Call M-SSA core routine for this Slepian index
        [MSSAcoff_rec, MSSAcoff_RC, MSSAcoff_evalues, RCcoff_Test] = ...
            MSSA_final_noCDF_freqsort( ...
                mssa_Sort_sn, ...
                fillM_deltacoffs(:,:,ss), ...
                M_rec, ...
                N_rec(ss), ...
                [], ...           % placeholder (unchanged from original code)
                Noise_SigLev);

        % Store outputs
        MSSAcoffs_rec(:,:,ss)   = MSSAcoff_rec;
        MSSAcoffs_evalues(:,ss) = MSSAcoff_evalues;
        MSSAcoffs_RC{ss}        = MSSAcoff_RC;
        MSSAcoffs_RCTest{ss}    = RCcoff_Test;
    end
end
