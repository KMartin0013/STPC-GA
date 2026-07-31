function [S_shannon, S, S_sig] = ...
    choose_S_from_eigenvalues(selection, Max_S, V, Lwindow, XY_buffer)
% CHOOSE_S_FROM_EIGENVALUES Determine V3 S and selective-smoothing bounds.
%
% Inputs:
%   selection  - normalized V3 S-selection struct
%   Max_S      - working S limit when selection.mode is 'stpc'
%   V          - eigenvalues of Slepian basis
%   Lwindow    - maximum spherical harmonic degree
%   XY_buffer  - buffered region polygon
% Outputs:
%   S_shannon  - Shannon number
%   S          - chosen number of Slepian functions
%   S_sig      - either S or [N1, N2] (for smoothing cases)

    % Direct callers from older code may still pass S_choice.  V3 uses the
    % renumbered scheme requested for compatibility: 0 Shannon, 1 STPC,
    % 2 eigenvalue threshold, 3 selective Gaussian.  V2p3's 6 remains an
    % unambiguous STPC alias; removed methods 4 and 5 are rejected.
    if isnumeric(selection) && isscalar(selection) && ...
            isfinite(selection) && selection == fix(selection)
        switch double(selection)
            case 0
                selection = struct('mode', 'shannon');
            case {1, 6}
                selection = struct('mode', 'stpc');
            case 2
                selection = struct('mode', 'eigen_threshold', ...
                    'threshold', 0.1);
            case 3
                selection = struct('mode', 'selective_gaussian', ...
                    'threshold', 0.1);
            otherwise
                error('V3:RemovedSChoice', ...
                    'S_choice=%d is not supported in V3.', selection);
        end
    elseif ~isstruct(selection) || ~isfield(selection, 'mode')
        error('V3:InvalidSSelection', ...
            'S selection must be a normalized struct or V3 S_choice number.');
    end

    S_shannon = round((Lwindow + 1)^2 * spharea(XY_buffer));
    switch selection.mode
        case 'stpc'
            S = Max_S;
            S_sig = S;
        case 'shannon'
            S     = S_shannon;
            S_sig = S;
        case 'eigen_threshold'
            S     = sum(V > selection.threshold);
            S_sig = S;
        case 'selective_gaussian'
            S = sum(V > selection.threshold);
            S_sig = [S_shannon, S];
            if S_shannon > S
                error('V3:SelectiveGaussianBounds', ...
                    ['Selective Gaussian requires the eigenvalue-threshold ', ...
                     'upper limit S=%d to be at least Shannon S=%d.'], ...
                    S, S_shannon);
            end
        case 'fixed'
            S = selection.value;
            S_sig = S;
        otherwise
            error('V3:UnknownSMode', ...
                'Unknown normalized S-selection mode: %s', selection.mode);
    end

    if ~(isscalar(S) && isfinite(S) && S > 0 && S == fix(S))
        error('V3:InvalidResolvedS', 'Resolved S must be a positive integer.');
    end
    if S > numel(V)
        error('V3:SExceedsBasis', ...
            'Resolved S=%d exceeds the available %d Slepian functions.', ...
            S, numel(V));
    end
end
