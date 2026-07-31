function decision = buffer_selection_diagnostics(group_buffer, buffer_pcc, ...
    buffer_rmse, trend_est, amplitude_est, phase_est, true_trend, ...
    true_amplitude, true_phase, index_pcc, index_rmse)
% Compare buffer candidates and provide an automatic conflict resolution.
%
% PCC and RMSE remain the primary selectors. If they disagree, the
% automatic choice is made between those two candidates by counting which
% one has the smaller absolute error in trend, annual amplitude, and annual
% phase. A normalized-error score breaks a tied vote; RMSE is the final
% deterministic fallback.

group_buffer = group_buffer(:);
buffer_pcc = buffer_pcc(:);
buffer_rmse = buffer_rmse(:);
trend_est = trend_est(:);
amplitude_est = amplitude_est(:);
phase_est = phase_est(:);

nCandidates = numel(group_buffer);
assert(all([numel(buffer_pcc), numel(buffer_rmse), numel(trend_est), ...
    numel(amplitude_est), numel(phase_est)] == nCandidates), ...
    'All buffer diagnostic vectors must have the same length.');
assert(isscalar(index_pcc) && ismember(index_pcc, 1:nCandidates), ...
    'index_pcc is outside the candidate range.');
assert(isscalar(index_rmse) && ismember(index_rmse, 1:nCandidates), ...
    'index_rmse is outside the candidate range.');

trendError = abs(trend_est - true_trend);
amplitudeError = abs(amplitude_est - true_amplitude);
phaseError = abs(mod(phase_est - true_phase + 180, 360) - 180);

candidateNumber = (1:nCandidates)';
isPccBest = candidateNumber == index_pcc;
isRmseBest = candidateNumber == index_rmse;
diagnosticTable = table(candidateNumber, group_buffer, buffer_pcc, ...
    buffer_rmse, trend_est, repmat(true_trend, nCandidates, 1), ...
    trendError, amplitude_est, repmat(true_amplitude, nCandidates, 1), ...
    amplitudeError, phase_est, repmat(true_phase, nCandidates, 1), ...
    phaseError, isPccBest, isRmseBest, ...
    'VariableNames', {'candidate', 'buffer_deg', 'pcc', 'rmse_cm', ...
    'trend_cm_per_year', 'true_trend_cm_per_year', ...
    'trend_abs_error_cm_per_year', 'annual_amplitude_cm', ...
    'true_annual_amplitude_cm', 'annual_amplitude_abs_error_cm', ...
    'annual_phase_deg', 'true_annual_phase_deg', ...
    'annual_phase_abs_error_deg', 'is_pcc_best', 'is_rmse_best'});

decision = struct();
decision.table = diagnosticTable;
decision.has_conflict = index_pcc ~= index_rmse;
decision.pcc_index = index_pcc;
decision.rmse_index = index_rmse;

if ~decision.has_conflict
    decision.auto_index = index_pcc;
    decision.auto_buffer_deg = group_buffer(index_pcc);
    decision.pcc_candidate_wins = 3;
    decision.rmse_candidate_wins = 3;
    decision.reason = sprintf( ...
        'PCC and RMSE agree on buffer %.6g deg.', ...
        decision.auto_buffer_deg);
    return
end

pccErrors = [trendError(index_pcc), amplitudeError(index_pcc), ...
    phaseError(index_pcc)];
rmseErrors = [trendError(index_rmse), amplitudeError(index_rmse), ...
    phaseError(index_rmse)];

pccWins = sum(pccErrors < rmseErrors);
rmseWins = sum(rmseErrors < pccErrors);

if pccWins > rmseWins
    autoIndex = index_pcc;
    tieBreak = '';
elseif rmseWins > pccWins
    autoIndex = index_rmse;
    tieBreak = '';
else
    errorScale = max([pccErrors; rmseErrors], [], 1);
    errorScale(errorScale == 0 | ~isfinite(errorScale)) = 1;
    pccScore = sum(pccErrors ./ errorScale);
    rmseScore = sum(rmseErrors ./ errorScale);
    if pccScore < rmseScore
        autoIndex = index_pcc;
        tieBreak = sprintf(' Normalized-error tie-break: %.6g < %.6g.', ...
            pccScore, rmseScore);
    elseif rmseScore < pccScore
        autoIndex = index_rmse;
        tieBreak = sprintf(' Normalized-error tie-break: %.6g < %.6g.', ...
            rmseScore, pccScore);
    else
        autoIndex = index_rmse;
        tieBreak = ' Exact tie: use the RMSE candidate.';
    end
end

decision.auto_index = autoIndex;
decision.auto_buffer_deg = group_buffer(autoIndex);
decision.pcc_candidate_wins = pccWins;
decision.rmse_candidate_wins = rmseWins;
decision.reason = sprintf([ ...
    'Secondary diagnostics vote (trend/amplitude/phase): PCC candidate ' ...
    'buffer %.6g deg wins %d, RMSE candidate buffer %.6g deg wins %d. ' ...
    'Automatic choice: %.6g deg.%s'], ...
    group_buffer(index_pcc), pccWins, group_buffer(index_rmse), rmseWins, ...
    decision.auto_buffer_deg, tieBreak);
end
