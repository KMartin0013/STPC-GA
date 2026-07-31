function [cache, timing] = prepare_stpc_decomposition_cache( ...
    mssa_results, turn_V_four, turn_MSSA_four, used_sta_SN, ...
    Turning_number, Noise_SigLev)
%PREPARE_STPC_DECOMPOSITION_CACHE Compute p-independent MSSA/SVD once.
%
% The full SVD, RCs, spectra and p-values do not depend on the selected
% significance level.  V3 computes them once at the maximum candidate
% S and N.  run_MSSA_decompose then slices the cached RCs and reapplies the
% requested alpha independently for every p. Candidate counts are derived
% from the actual arrays, so V3 supports K-by-K, 1-by-K, K-by-1 and 1-by-1.

mssa_Sort        = mssa_results.mssa_Sort;
fillM_deltacoffs = mssa_results.fillM_deltacoffs;
M_rec            = mssa_results.M_rec;
intit_num        = mssa_results.intit_num;
fill_nmonths     = size(fillM_deltacoffs, 1);

numSCandidates = size(turn_V_four, 1);
S_rec_max = max(turn_V_four(1:numSCandidates, used_sta_SN));
isOcean = strcmp(mssa_Sort(1).process{3}, 'ocean');
S_ol_max = S_rec_max + double(isOcean);

N_max = zeros(1, S_ol_max);
for ss = 1:S_ol_max
    numNCandidates = size(turn_MSSA_four{ss}, 1);
    N_max(ss) = max(turn_MSSA_four{ss}( ...
        1:numNCandidates, used_sta_SN));
end

availableN = intit_num * M_rec;
if any(N_max > availableN)
    badIndex = find(N_max > availableN, 1);
    error('V3:FixedNOutOfRange', ...
        ['Requested N=%d for Slepian coefficient %d, but at most %d ', ...
         'MSSA modes are available for M=%d and %d institutions.'], ...
        N_max(badIndex), badIndex, availableN, M_rec, intit_num);
end

lowState = warning('query', 'stats:lillietest:OutOfRangePLow');
highState = warning('query', 'stats:lillietest:OutOfRangePHigh');
cleanupWarnings = onCleanup(@() local_restore_warnings(lowState, highState)); %#ok<NASGU>
warning('off', 'stats:lillietest:OutOfRangePLow');
warning('off', 'stats:lillietest:OutOfRangePHigh');

timer = tic;
[~, evalues, RC, RCTest] = run_MSSA_decompose( ...
    mssa_Sort, fillM_deltacoffs, M_rec, N_max, Noise_SigLev, ...
    fill_nmonths, intit_num, S_ol_max, S_rec_max);
timing = toc(timer);

cache = struct();
cache.version = 'V3';
cache.S_rec_max = S_rec_max;
cache.S_ol_max = S_ol_max;
cache.N_max = N_max;
cache.evalues = evalues;
cache.RC = RC;
cache.RCTest = RCTest;
end

function local_restore_warnings(lowState, highState)
warning(lowState.state, 'stats:lillietest:OutOfRangePLow');
warning(highState.state, 'stats:lillietest:OutOfRangePHigh');
end
