function [evalues, st_eofs, st_pcs, RC] = MSSA_fast(time, data, M, Method, N, opts)
% MSSA_FAST
%
% Compatible fast version of MSSA.m.
%
% Outputs keep the same meanings and dimensions as MSSA.m for st_eofs,
% st_pcs, and RC:
%   st_eofs : (M * channel_num) x N
%   st_pcs  : (time_num - M + 1) x N
%   RC      : time_num x N x channel_num
%
% evalues:
%   opts.returnFullEvalues = true
%       Return all normalized singular values, like MSSA.m.
%       This is the safest mode for code that detects turning points from
%       the eigenvalue spectrum.
%
%   opts.returnFullEvalues = false
%       Return only the first N normalized singular values. The denominator
%       is still the sum of all singular values, so each selected value is
%       its absolute fraction of the full spectrum.
%
% Method == 2:
%   T = d' * d / N1, symmetric positive semi-definite.
%   sum(all singular values) = trace(T).
%
% Method == 1:
%   T may not be strictly symmetric positive semi-definite numerically.
%   To keep the original SVD-based definition exactly, denominator uses
%   sum(svd(T,'econ')).

if nargin < 6
    opts = struct();
end
if ~isfield(opts, 'verbose'), opts.verbose = false; end
if ~isfield(opts, 'useSvds'), opts.useSvds = true; end
if ~isfield(opts, 'returnFullEvalues'), opts.returnFullEvalues = true; end
if ~isfield(opts, 'svdsTolerance'), opts.svdsTolerance = []; end
if ~isfield(opts, 'svdsMaxIterations'), opts.svdsMaxIterations = []; end

tStart = tic;

[H, L] = size(data);
K = H - M + 1;

if K < 1
    error('MSSA_fast:InvalidM', ...
        'M must be smaller than or equal to the time-series length.');
end

T = GrandMatrix(data, M, Method);

matrix_size = min(size(T));
N = min(N, matrix_size);

% ===== singular values =====
% In full-spectrum mode, use complete SVD. This preserves the old MSSA.m
% evalues output exactly enough for downstream turning-point logic.
if opts.returnFullEvalues

    [Uall, Sall, ~] = svd(T, 'econ');

    all_singular_values = diag(Sall);
    [all_singular_values, full_sort_index] = sort(all_singular_values, 'descend');

    total_singular_sum = sum(all_singular_values);
    if total_singular_sum > 0
        evalues = all_singular_values / total_singular_sum;
    else
        evalues = all_singular_values;
    end

    U = Uall(:, full_sort_index(1:N));

else
    % ===== denominator: sum of all singular values =====
    if Method == 2
        % For T = d' * d / K, T is covariance-like and PSD.
        total_singular_sum = trace(T);

        % Numerical guard.
        if total_singular_sum <= 0 || ~isfinite(total_singular_sum)
            total_singular_sum = sum(svd(T, 'econ'));
        end
    else
        % Exact denominator under the original SVD definition.
        total_singular_sum = sum(svd(T, 'econ'));
    end

    % ===== numerator and EOFs: first N singular components =====
    if opts.useSvds && N < matrix_size
        try
            svds_opts = struct();
            if ~isempty(opts.svdsTolerance)
                svds_opts.tol = opts.svdsTolerance;
            end
            if ~isempty(opts.svdsMaxIterations)
                svds_opts.maxit = opts.svdsMaxIterations;
            end

            if isempty(fieldnames(svds_opts))
                [U, S, ~] = svds(T, N);
            else
                [U, S, ~] = svds(T, N, 'largest', svds_opts);
            end
        catch
            [Uall, Sall, ~] = svd(T, 'econ');
            U = Uall(:,1:N);
            S = Sall(1:N,1:N);
        end
    else
        [Uall, Sall, ~] = svd(T, 'econ');
        U = Uall(:,1:N);
        S = Sall(1:N,1:N);
    end

    singular_values = diag(S);
    [singular_values, sort_index] = sort(singular_values, 'descend');
    U = U(:, sort_index);

    if total_singular_sum > 0
        evalues = singular_values / total_singular_sum;
    else
        evalues = singular_values;
    end
end

st_eofs = U(:,1:N);

% ===== PCs =====
X = zeros(K, M * L);

for l = 1:L
    cols = (l-1)*M + (1:M);

    for lag = 1:M
        X(:, cols(lag)) = data(lag:lag+K-1, l);
    end
end

st_pcs = X * st_eofs;

% ===== RC: diagonal averaging =====
RC = zeros(H, N, L);

idx = (1:K).' + (0:M-1);
counts = accumarray(idx(:), 1, [H, 1]);

for l = 1:L
    rows = (l-1)*M + (1:M);

    for k = 1:N
        eof_lk = st_eofs(rows, k);
        vals = st_pcs(:,k) * eof_lk.';
        RC(:,k,l) = accumarray(idx(:), vals(:), [H, 1]) ./ counts;
    end
end

if opts.verbose
    fprintf(1, 'M-SSA fast computational time = %f seconds\n', toc(tStart));
end

end
