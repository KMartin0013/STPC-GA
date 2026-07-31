function [wcorr, weights] = compute_mssa_wcorrelation(RC, windowLength)
%COMPUTE_MSSA_WCORRELATION Joint-channel w-correlation of MSSA RCs.
%
% RC is time x component x channel. The time weights equal the number of
% appearances of each sample in the SSA trajectory matrix. Inner products
% are summed across channels, matching the joint MSSA trajectory space.

arguments
    RC (:,:,:) double
    windowLength (1,1) double {mustBeInteger, mustBePositive}
end

timeCount = size(RC, 1);
componentCount = size(RC, 2);
if windowLength > timeCount
    error('V3:WcorrWindowTooLong', ...
        'MSSA window M=%d exceeds the RC time length %d.', ...
        windowLength, timeCount);
end
if componentCount < 1
    error('V3:WcorrNoComponents', ...
        'At least one reconstructed component is required.');
end
if any(~isfinite(RC), 'all')
    error('V3:WcorrNonfiniteRC', ...
        'MSSA reconstructed components contain nonfinite values.');
end

laggedColumnCount = timeCount - windowLength + 1;
weights = conv(ones(windowLength, 1), ones(laggedColumnCount, 1));

weightedRC = RC .* reshape(sqrt(weights), [], 1, 1);
jointRC = reshape(permute(weightedRC, [1, 3, 2]), [], componentCount);
gram = jointRC' * jointRC;
norms = sqrt(max(diag(gram), 0));
denominator = norms * norms';

wcorr = zeros(componentCount);
valid = denominator > 0;
wcorr(valid) = gram(valid) ./ denominator(valid);
wcorr = max(min(wcorr, 1), -1);
wcorr = (wcorr + wcorr') / 2;
wcorr(1:componentCount+1:end) = double(norms > 0);
end
