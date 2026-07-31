function [turn_MSSA_four, MSSA_evalues_sumup, MSSA_evalues, ...
    wcorrDiagnostics] = det_N_WCorr(mssa_results, N_bou, ...
    Turning_number, Attach, plotProcess, N_selection)
%DET_N_WCORR Generate N candidates from joint-channel w-correlation.
%
% W-correlation supplies alternative candidate cutoffs. It does not choose
% the final N: the existing STPC spatial-temporal screening remains intact.

intit_num = mssa_results.intit_num;
fillM_deltacoffs = mssa_results.fillM_deltacoffs;
S_ol = mssa_results.S_ol;
M_rec = mssa_results.M_rec;
timeCount = size(fillM_deltacoffs, 1);
availableModes = intit_num * M_rec;

requestedMaxModes = N_selection.max_modes;
if isempty(requestedMaxModes)
    scanModes = min(M_rec, availableModes);
else
    scanModes = min(requestedMaxModes, availableModes);
end
if scanModes < 5
    error('V3:WcorrTooFewModes', ...
        'At least five MSSA modes are required for w-correlation selection.');
end

opts = struct( ...
    'useSvds', true, ...
    'returnFullEvalues', true, ...
    'verbose', false);

MSSA_evalues = zeros(availableModes, S_ol);
MSSA_evalues_sumup = zeros(availableModes, S_ol);
turn_MSSA_four = cell(1, S_ol);
coefficientDiagnostics = cell(1, S_ol);
weights = [];

for ss = 1:S_ol
    [evalues, ~, ~, RC] = MSSA_fast( ...
        1:timeCount, fillM_deltacoffs(:, 1:intit_num, ss), ...
        M_rec, 2, scanModes, opts);
    evalues = evalues(:);
    if numel(evalues) ~= availableModes
        error('V3:UnexpectedWcorrSpectrumLength', ...
            ['MSSA returned %d eigenvalues for coefficient %d; ' ...
             '%d were expected.'], numel(evalues), ss, availableModes);
    end
    MSSA_evalues(:, ss) = evalues;
    MSSA_evalues_sumup(:, ss) = cumsum(evalues);

    rankTolerance = numel(evalues) * eps(max(evalues));
    numericalRank = find(evalues > rankTolerance, 1, 'last');
    if isempty(numericalRank)
        error('V3:ZeroWcorrSpectrum', ...
            'MSSA coefficient %d has no positive numerical rank.', ss);
    end
    modeCount = min([scanModes, numericalRank, size(RC, 2)]);
    [wcorr, coefficientWeights] = compute_mssa_wcorrelation( ...
        RC(:, 1:modeCount, 1:intit_num), M_rec);
    if isempty(weights)
        weights = coefficientWeights;
    end

    [candidates, selectionDiagnostics] = select_wcorr_N_candidates( ...
        wcorr, evalues(1:modeCount), N_bou, Turning_number, ...
        'blockWidth', N_selection.block_width, ...
        'minSpacing', N_selection.min_spacing);
    turn_MSSA_four{ss} = repmat(candidates, 1, 4);

    selectionDiagnostics.ssf = ss;
    selectionDiagnostics.numerical_rank = numericalRank;
    selectionDiagnostics.scan_modes = scanModes;
    selectionDiagnostics.wcorr = wcorr;
    coefficientDiagnostics{ss} = selectionDiagnostics;
end

wcorrDiagnostics = struct();
wcorrDiagnostics.version = 'V3-wcorr-1';
wcorrDiagnostics.method = N_selection.algorithm;
wcorrDiagnostics.channel_aggregation = 'joint_weighted_inner_product';
wcorrDiagnostics.requested_max_modes = requestedMaxModes;
wcorrDiagnostics.scan_modes = scanModes;
wcorrDiagnostics.requested_block_width = N_selection.block_width;
wcorrDiagnostics.min_spacing = N_selection.min_spacing;
wcorrDiagnostics.weights = weights;
wcorrDiagnostics.coefficients = coefficientDiagnostics;

if plotProcess
    plot_mssa_wcorrelation(wcorrDiagnostics, turn_MSSA_four, Attach);
end
end
