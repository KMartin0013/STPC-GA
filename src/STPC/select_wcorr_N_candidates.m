function [candidates, diagnostics] = select_wcorr_N_candidates( ...
    wcorr, eigenvalues, NBound, candidateCount, options)
%SELECT_WCORR_N_CANDIDATES Select contiguous RC cutoffs from w-correlation.
%
% A candidate boundary separates the last local block of retained RCs from
% the first local block of residual RCs. A strong boundary has relatively
% low cross-block w-correlation and relatively high correlation inside the
% following block. The strongest distinct local maxima are retained as N
% candidates; final selection remains the responsibility of STPC.

arguments
    wcorr (:,:) double
    eigenvalues (:,1) double
    NBound (1,1) double {mustBeNonnegative}
    candidateCount (1,1) double {mustBeInteger, mustBePositive}
    options.blockWidth (1,1) double {mustBeInteger, mustBePositive} = 6
    options.minSpacing (1,1) double {mustBeInteger, mustBePositive} = 2
end

if size(wcorr, 1) ~= size(wcorr, 2)
    error('V3:WcorrNotSquare', 'The w-correlation matrix must be square.');
end
if numel(eigenvalues) < size(wcorr, 1)
    error('V3:WcorrEigenvalueLength', ...
        'The eigenvalue vector is shorter than the w-correlation matrix.');
end
if ~(isfinite(NBound) && NBound < 0.5)
    error('V3:InvalidWcorrNBound', ...
        'N_bound must be finite and satisfy 0 <= N_bound < 0.5.');
end
if any(~isfinite(wcorr), 'all') || any(~isfinite(eigenvalues))
    error('V3:WcorrNonfiniteInput', ...
        'W-correlation candidate inputs must be finite.');
end

modeCount = size(wcorr, 1);
blockWidth = min(options.blockWidth, max(2, floor((modeCount - 1) / 3)));
domain = (blockWidth:(modeCount - blockWidth))';
if isempty(domain)
    error('V3:WcorrDomainTooShort', ...
        'Only %d modes are available for w-correlation candidate selection.', ...
        modeCount);
end

absWcorr = abs(wcorr);
absWcorr = min(absWcorr, 1);
cumulative = cumsum(max(eigenvalues(1:modeCount), 0));
if cumulative(end) > 0
    cumulative = cumulative / cumulative(end);
end

score = nan(size(domain));
crossRms = nan(size(domain));
tailRms = nan(size(domain));
globalCrossRms = nan(size(domain));
validDomain = false(size(domain));

for i = 1:numel(domain)
    cutoff = domain(i);
    left = cutoff-blockWidth+1:cutoff;
    right = cutoff+1:cutoff+blockWidth;

    crossBlock = absWcorr(left, right);
    tailBlock = absWcorr(right, right);
    tailBlock(1:blockWidth+1:end) = NaN;

    crossRms(i) = sqrt(mean(crossBlock.^2, 'all'));
    tailRms(i) = sqrt(mean(tailBlock.^2, 'all', 'omitnan'));
    globalBlock = absWcorr(1:cutoff, cutoff+1:modeCount);
    globalCrossRms(i) = sqrt(mean(globalBlock.^2, 'all'));
    score(i) = tailRms(i) - crossRms(i);
    validDomain(i) = cumulative(cutoff) > NBound && ...
        cumulative(cutoff) < 1-NBound;
end

validIndex = find(validDomain & isfinite(score));
if isempty(validIndex)
    error('V3:NoWcorrCandidateDomain', ...
        ['No w-correlation boundary falls inside N_bound. ' ...
         'Reduce config.N_bound or adjust config.Wcorr.maxModes.']);
end

localMaximum = false(size(score));
for i = 1:numel(score)
    previous = -Inf;
    following = -Inf;
    if i > 1, previous = score(i-1); end
    if i < numel(score), following = score(i+1); end
    localMaximum(i) = score(i) >= previous && score(i) >= following;
end

localIndex = find(validDomain & localMaximum & score > 0);
[~, localOrder] = sort(score(localIndex), 'descend');
rankedIndex = localIndex(localOrder);

[~, allOrder] = sort(score(validIndex), 'descend');
rankedIndex = [rankedIndex; validIndex(allOrder)];
rankedIndex = unique(rankedIndex, 'stable');

selectedIndex = zeros(0, 1);
for i = 1:numel(rankedIndex)
    candidateIndex = rankedIndex(i);
    candidateN = domain(candidateIndex);
    if isempty(selectedIndex) || all(abs(candidateN-domain(selectedIndex)) >= ...
            options.minSpacing)
        selectedIndex(end+1, 1) = candidateIndex; %#ok<AGROW>
    end
    if numel(selectedIndex) == candidateCount
        break
    end
end

if numel(selectedIndex) < candidateCount
    for i = 1:numel(rankedIndex)
        candidateIndex = rankedIndex(i);
        if ~ismember(candidateIndex, selectedIndex)
            selectedIndex(end+1, 1) = candidateIndex; %#ok<AGROW>
        end
        if numel(selectedIndex) == candidateCount
            break
        end
    end
end

if isempty(selectedIndex)
    error('V3:NoWcorrCandidates', ...
        'No w-correlation N candidate could be selected.');
end

[bestScore, bestLocalIndex] = max(score(validIndex));
bestIndex = validIndex(bestLocalIndex);
selectedN = sort(domain(selectedIndex));
distinctCount = numel(selectedN);
if distinctCount < candidateCount
    selectedN(end+1:candidateCount, 1) = selectedN(end);
    warning('V3:FewerWcorrCandidates', ...
        ['Only %d distinct w-correlation N candidates were available; ' ...
         'the final candidate is repeated to keep turningNumber=%d.'], ...
        distinctCount, candidateCount);
end
candidates = selectedN(:);

diagnostics = struct();
diagnostics.method = 'fixed_block_contrast_v1';
diagnostics.mode_count = modeCount;
diagnostics.block_width = blockWidth;
diagnostics.min_spacing = options.minSpacing;
diagnostics.domain = domain;
diagnostics.valid_domain = validDomain;
diagnostics.score = score;
diagnostics.cross_rms = crossRms;
diagnostics.tail_rms = tailRms;
diagnostics.global_cross_rms = globalCrossRms;
diagnostics.cumulative_eigenvalues = cumulative;
diagnostics.candidates = candidates;
diagnostics.best_candidate = domain(bestIndex);
diagnostics.best_score = bestScore;
diagnostics.best_cross_rms = crossRms(bestIndex);
diagnostics.best_tail_rms = tailRms(bestIndex);
diagnostics.low_confidence = bestScore <= 0;
end
