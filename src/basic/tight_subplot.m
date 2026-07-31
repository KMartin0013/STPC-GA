function [axesHandles, positions] = tight_subplot( ...
        nRows, nColumns, gap, marginHeight, marginWidth)
% Create a compact grid of axes with explicit gaps and outer margins.
%
% Inputs follow the established tight_subplot interface:
%   gap          - scalar or [vertical horizontal]
%   marginHeight - scalar or [bottom top]
%   marginWidth  - scalar or [left right]

arguments
    nRows (1,1) double {mustBeInteger, mustBePositive}
    nColumns (1,1) double {mustBeInteger, mustBePositive}
    gap (1,:) double = 0.02
    marginHeight (1,:) double = 0.05
    marginWidth (1,:) double = 0.05
end

gap = local_expand_pair(gap, 'gap');
marginHeight = local_expand_pair(marginHeight, 'marginHeight');
marginWidth = local_expand_pair(marginWidth, 'marginWidth');

axesHeight = (1 - sum(marginHeight) - ...
    (nRows - 1) * gap(1)) / nRows;
axesWidth = (1 - sum(marginWidth) - ...
    (nColumns - 1) * gap(2)) / nColumns;
if axesHeight <= 0 || axesWidth <= 0
    error('tight_subplot:InvalidLayout', ...
        'Gaps and margins leave no positive area for the axes.');
end

nAxes = nRows * nColumns;
axesHandles = gobjects(nAxes, 1);
positions = zeros(nAxes, 4);
axesIndex = 0;
yPosition = 1 - marginHeight(2) - axesHeight;

for rowIndex = 1:nRows
    xPosition = marginWidth(1);
    for columnIndex = 1:nColumns
        axesIndex = axesIndex + 1;
        positions(axesIndex, :) = [ ...
            xPosition, yPosition, axesWidth, axesHeight];
        axesHandles(axesIndex) = axes( ...
            'Units', 'normalized', ...
            'Position', positions(axesIndex, :)); %#ok<LAXES>
        xPosition = xPosition + axesWidth + gap(2);
    end
    yPosition = yPosition - axesHeight - gap(1);
end
end

function pair = local_expand_pair(value, name)
if isscalar(value)
    pair = [value value];
elseif numel(value) == 2
    pair = reshape(value, 1, 2);
else
    error('tight_subplot:InvalidSize', ...
        '%s must be a scalar or a two-element vector.', name);
end
if any(~isfinite(pair)) || any(pair < 0)
    error('tight_subplot:InvalidValue', ...
        '%s values must be finite and nonnegative.', name);
end
end
