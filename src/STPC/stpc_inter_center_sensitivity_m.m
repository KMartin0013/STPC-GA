function sensitivity = stpc_inter_center_sensitivity_m(center_stpc, fitwhat, params)
%STPC_INTER_CENTER_SENSITIVITY_M Inter-center spread for STPC results.
%
% Inter-center spread quantifies the sensitivity of the STPC result to the
% choice of GRACE/GRACE-FO Level-2 solution. It is not interpreted as a
% formal observational standard error or confidence interval.
%
% Inputs:
%   center_stpc - output of run_stpc_center_specific_selected.
%   fitwhat     - harmonic fit definition, e.g. [3 365.25 365.25/2].
%   params      - optional struct. Set params.compute_grid = true to also
%                 compute grid-cell trend spread maps.

if nargin < 2 || isempty(fitwhat)
    fitwhat = [3 365.25 365.25/2];
end
if nargin < 3 || isempty(params)
    params = struct();
end
if ~isfield(params, 'compute_grid')
    params.compute_grid = false;
end
if ~isfield(params, 'fit_params')
    params.fit_params = struct();
end
params.fit_params.n_param_draws = 0;

Ycenter = center_stpc.EWH.fillM_STPC_center;
dates = center_stpc.fill_dates(:);
[~, n_centers, n_p] = size(Ycenter);
num_harmonics = numel(fitwhat) - 1;

sensitivity = struct();
sensitivity.method = 'inter-center sensitivity';
sensitivity.center_names = center_stpc.center_names;
sensitivity.n_centers = n_centers;
sensitivity.p_use = center_stpc.p_use;
sensitivity.fill_dates = dates;
sensitivity.EWH.unit = center_stpc.EWH.unit;

sensitivity.EWH.monthly_mean = local_nanmean(Ycenter, 2);
sensitivity.EWH.monthly_std = local_nanstd(Ycenter, 2);

sensitivity.EWH.trend_values = nan(n_centers, n_p);
sensitivity.EWH.trend_std = nan(1, n_p);
sensitivity.EWH.amplitude_values = nan(num_harmonics, n_centers, n_p);
sensitivity.EWH.amplitude_std = nan(num_harmonics, n_p);
sensitivity.EWH.phase_values = nan(num_harmonics, n_centers, n_p);
sensitivity.EWH.phase_circular_std = nan(num_harmonics, n_p);
sensitivity.EWH.peak_time_values = nan(num_harmonics, n_centers, n_p);
sensitivity.EWH.peak_time_std = nan(num_harmonics, n_p);

for p_i = 1:n_p
    for center_i = 1:n_centers
        y = squeeze(Ycenter(:, center_i, p_i));
        unc = local_harmonic_hac_uncertainty_m(y, dates, fitwhat, params.fit_params);
        sensitivity.EWH.trend_values(center_i, p_i) = unc.trend;
        sensitivity.EWH.amplitude_values(:, center_i, p_i) = unc.amplitude;
        sensitivity.EWH.phase_values(:, center_i, p_i) = unc.phase;
        sensitivity.EWH.peak_time_values(:, center_i, p_i) = unc.peak_time_days;
    end
    sensitivity.EWH.trend_std(1, p_i) = ...
        local_nanstd(sensitivity.EWH.trend_values(:, p_i), 1);
    for harmonic_i = 1:num_harmonics
        sensitivity.EWH.amplitude_std(harmonic_i, p_i) = ...
            local_nanstd(sensitivity.EWH.amplitude_values(harmonic_i, :, p_i), 2);
        sensitivity.EWH.phase_circular_std(harmonic_i, p_i) = ...
            local_circular_std(sensitivity.EWH.phase_values(harmonic_i, :, p_i));
        sensitivity.EWH.peak_time_std(harmonic_i, p_i) = ...
            local_nanstd(sensitivity.EWH.peak_time_values(harmonic_i, :, p_i), 2);
    end
end

if params.compute_grid && isfield(center_stpc.EWH, 'fillM_Grid_STPC_center')
    grid = center_stpc.EWH.fillM_Grid_STPC_center;
    sensitivity.EWH.monthly_grid_mean = local_nanmean(grid, 4);
    sensitivity.EWH.monthly_grid_std = local_nanstd(grid, 4);
    [~, nlat, nlon, ~, ~] = size(grid);
    sensitivity.EWH.trend_grid_values = nan(nlat, nlon, n_centers, n_p);
    sensitivity.EWH.trend_grid_std = nan(nlat, nlon, n_p);

    for p_i = 1:n_p
        for center_i = 1:n_centers
            for lat_i = 1:nlat
                for lon_i = 1:nlon
                    y = squeeze(grid(:, lat_i, lon_i, center_i, p_i));
                    if any(~isfinite(y))
                        continue
                    end
                    unc = local_harmonic_hac_uncertainty_m(y, dates, fitwhat, params.fit_params);
                    sensitivity.EWH.trend_grid_values(lat_i, lon_i, center_i, p_i) = unc.trend;
                end
            end
        end
        sensitivity.EWH.trend_grid_std(:, :, p_i) = ...
            local_nanstd(sensitivity.EWH.trend_grid_values(:, :, :, p_i), 3);
    end
end
end

function m = local_nanmean(x, dim)
valid = isfinite(x);
x(~valid) = 0;
count = sum(valid, dim);
m = sum(x, dim) ./ count;
m(count == 0) = nan;
end

function s = local_nanstd(x, dim)
m = local_nanmean(x, dim);
rep_size = ones(1, ndims(x));
rep_size(dim) = size(x, dim);
delta = x - repmat(m, rep_size);
delta(~isfinite(delta)) = 0;
count = sum(isfinite(x), dim);
variance = sum(delta.^2, dim) ./ max(count - 1, 1);
s = sqrt(variance);
s(count <= 1) = nan;
end

function circular_std = local_circular_std(phase_radians)
phase_radians = phase_radians(isfinite(phase_radians));
if isempty(phase_radians)
    circular_std = nan;
    return
end
mean_resultant_length = abs(mean(exp(1i * phase_radians)));
mean_resultant_length = min(max(mean_resultant_length, realmin), 1);
circular_std = sqrt(-2 * log(mean_resultant_length));
end
