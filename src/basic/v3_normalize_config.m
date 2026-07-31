function config = v3_normalize_config(config)
%V3_NORMALIZE_CONFIG Normalize the public R/S/M/N configuration interface.
%
% Public V3 interface:
%   config.R = "STPC" | numeric buffer in degrees
%   config.S = "STPC" | "Shannon" | ["Fixed", threshold] |
%              ["Sel_Gau_default", radius_km] | positive integer
%   config.M = "STPC" | positive integer MSSA window
%   config.N = "STPC" | "Wcorr" | positive integer retained RC count
%
% Legacy V2p3 fields are translated when the new field is absent.

selection = struct();

%% R: buffer selection
if isfield(config, 'R')
    rValue = config.R;
elseif isfield(config, 'buffer_deg')
    rValue = config.buffer_deg;
else
    error('V3:MissingR', 'Provide config.R as "STPC" or a numeric buffer in degrees.');
end

if local_is_stpc(rValue) || local_is_auto(rValue)
    selection.R = struct('mode', 'stpc', 'value', []);
    config.R = "STPC";
    config.buffer_deg = 'Auto';
elseif local_positive_or_signed_scalar(rValue)
    selection.R = struct('mode', 'fixed', 'value', double(rValue));
    config.R = double(rValue);
    config.buffer_deg = double(rValue);
else
    error('V3:InvalidR', 'config.R must be "STPC" or a finite numeric scalar.');
end

%% S: Slepian truncation / smoothing
if isfield(config, 'S')
    [selection.S, config] = local_parse_new_s(config.S, config);
elseif isfield(config, 'S_choice')
    [selection.S, config] = local_parse_legacy_s(config.S_choice, config);
else
    error('V3:MissingS', 'Provide config.S or a supported legacy config.S_choice.');
end
config.S = local_public_s(selection.S);
if ~isfield(config, 'Radius') || isempty(config.Radius)
    config.Radius = 500;
end

%% M: reconstruction-window selection
if ~isfield(config, 'M')
    error('V3:MissingM', 'Provide config.M as "STPC" or a positive integer.');
end
if local_is_stpc(config.M) || local_is_auto(config.M)
    selection.M = struct('mode', 'stpc', 'value', []);
    config.M = 'Auto';
elseif local_positive_integer_scalar(config.M)
    selection.M = struct('mode', 'fixed', 'value', double(config.M));
    config.M = double(config.M);
else
    error('V3:InvalidM', 'config.M must be "STPC" or a positive integer.');
end

%% N: retained MSSA RC count
if ~isfield(config, 'N')
    nValue = "STPC"; % V2p3-compatible default
else
    nValue = config.N;
end
if local_is_stpc(nValue) || local_is_auto(nValue)
    selection.N = struct('mode', 'stpc', 'value', []);
    config.N = "STPC";
elseif local_is_wcorr(nValue)
    wcorr = local_wcorr_options(config);
    selection.N = struct( ...
        'mode', 'wcorr', ...
        'value', [], ...
        'algorithm', 'fixed_block_contrast_v1', ...
        'block_width', wcorr.blockWidth, ...
        'min_spacing', wcorr.minSpacing, ...
        'max_modes', wcorr.maxModes);
    config.N = "Wcorr";
    config.Wcorr = wcorr;
elseif local_positive_integer_scalar(nValue)
    selection.N = struct('mode', 'fixed', 'value', double(nValue));
    config.N = double(nValue);
else
    error('V3:InvalidN', ...
        'config.N must be "STPC", "Wcorr", or a positive integer.');
end

%% Candidate count
if ~isfield(config, 'turningNumber')
    config.turningNumber = 5;
end
if ~local_positive_integer_scalar(config.turningNumber)
    error('V3:InvalidTurningNumber', ...
        'config.turningNumber must be a positive integer.');
end
config.turningNumber = double(config.turningNumber);
selection.turningNumber = config.turningNumber;

%% Configuration directories
if ~isfield(config, 'caseName') || strlength(string(config.caseName)) == 0
    config.caseName = 'paper_default';
end
caseName = char(string(config.caseName));
if isempty(regexp(caseName, '^[A-Za-z0-9_-]+$', 'once'))
    error('V3:InvalidCaseName', ...
        'config.caseName may contain only letters, numbers, "_" and "-".');
end
config.caseName = caseName;

if ~isfield(config, 'resultRoot') || isempty(config.resultRoot)
    if isfield(config, 'resultDir') && ~isempty(config.resultDir)
        config.resultRoot = config.resultDir;
    else
        error('V3:MissingResultRoot', 'Provide config.resultRoot.');
    end
end
if ~isfield(config, 'figureRoot') || isempty(config.figureRoot)
    if isfield(config, 'figureDir') && ~isempty(config.figureDir)
        config.figureRoot = config.figureDir;
    else
        error('V3:MissingFigureRoot', 'Provide config.figureRoot.');
    end
end

config.sharedResultDir = char(string(config.resultRoot));
config.resultDir = fullfile(config.sharedResultDir, ...
    'Configurations', config.caseName);
config.figureDir = fullfile(char(string(config.figureRoot)), ...
    'Configurations', config.caseName);
config.AddresultDir = fullfile(config.resultDir, 'AddData');
config.v3Selection = selection;
end

function [spec, config] = local_parse_new_s(value, config)
if local_positive_integer_scalar(value)
    spec = struct('mode', 'fixed', 'value', double(value), ...
        'threshold', [], 'radius', []);
    return
end

tokens = string(value);
if isempty(tokens)
    error('V3:InvalidS', 'config.S cannot be empty.');
end
mode = lower(strtrim(tokens(1)));
switch mode
    case "stpc"
        spec = struct('mode', 'stpc', 'value', [], ...
            'threshold', [], 'radius', []);
    case "shannon"
        spec = struct('mode', 'shannon', 'value', [], ...
            'threshold', [], 'radius', []);
    case "fixed"
        if numel(tokens) ~= 2
            error('V3:InvalidSFixed', ...
                'Use config.S = ["Fixed","0.1"] for an eigenvalue threshold.');
        end
        threshold = str2double(tokens(2));
        local_validate_threshold(threshold);
        spec = struct('mode', 'eigen_threshold', 'value', [], ...
            'threshold', threshold, 'radius', []);
    case "sel_gau_default"
        if numel(tokens) ~= 2
            error('V3:InvalidSSelectiveGaussian', ...
                'Use config.S = ["Sel_Gau_default","500"].');
        end
        radius = str2double(tokens(2));
        if ~(isscalar(radius) && isfinite(radius) && radius > 0)
            error('V3:InvalidGaussianRadius', ...
                'Selective-Gaussian radius must be a positive number in km.');
        end
        spec = struct('mode', 'selective_gaussian', 'value', [], ...
            'threshold', 0.1, 'radius', radius);
        config.Radius = radius;
    otherwise
        error('V3:InvalidS', 'Unsupported config.S mode: %s', tokens(1));
end
end

function [spec, config] = local_parse_legacy_s(choice, config)
if ~local_positive_or_zero_integer_scalar(choice)
    error('V3:InvalidLegacySChoice', 'Legacy config.S_choice must be an integer.');
end
switch double(choice)
    case 0
        spec = struct('mode', 'shannon', 'value', [], ...
            'threshold', [], 'radius', []);
    case 1
        spec = struct('mode', 'stpc', 'value', [], ...
            'threshold', [], 'radius', []);
    case 2
        spec = struct('mode', 'eigen_threshold', 'value', [], ...
            'threshold', 0.1, 'radius', []);
    case 3
        if ~isfield(config, 'Radius') || isempty(config.Radius)
            config.Radius = 500;
        end
        if ~(isnumeric(config.Radius) && isscalar(config.Radius) && ...
                isfinite(config.Radius) && config.Radius > 0)
            error('V3:InvalidGaussianRadius', ...
                'Legacy config.Radius must be a positive number in km.');
        end
        spec = struct('mode', 'selective_gaussian', 'value', [], ...
            'threshold', 0.1, 'radius', double(config.Radius));
    case 6
        % Compatibility alias for V2p3 startup files.
        spec = struct('mode', 'stpc', 'value', [], ...
            'threshold', [], 'radius', []);
    otherwise
        error('V3:RemovedLegacySChoice', ...
            ['Legacy S_choice=%d was removed in V3. Use config.S as ', ...
             '"STPC", "Shannon", ["Fixed","threshold"], ', ...
             '["Sel_Gau_default","radius"], or a positive integer.'], choice);
end
end

function value = local_public_s(spec)
switch spec.mode
    case 'stpc'
        value = "STPC";
    case 'shannon'
        value = "Shannon";
    case 'eigen_threshold'
        value = ["Fixed", string(spec.threshold)];
    case 'selective_gaussian'
        value = ["Sel_Gau_default", string(spec.radius)];
    case 'fixed'
        value = spec.value;
    otherwise
        error('V3:InternalSMode', 'Unknown normalized S mode: %s', spec.mode);
end
end

function tf = local_is_stpc(value)
tf = (ischar(value) || (isstring(value) && isscalar(value))) && ...
    strcmpi(strtrim(char(string(value))), 'STPC');
end

function tf = local_is_auto(value)
tf = (ischar(value) || (isstring(value) && isscalar(value))) && ...
    strcmpi(strtrim(char(string(value))), 'Auto');
end

function tf = local_is_wcorr(value)
tf = (ischar(value) || (isstring(value) && isscalar(value))) && ...
    strcmpi(strtrim(char(string(value))), 'Wcorr');
end

function wcorr = local_wcorr_options(config)
wcorr = struct('blockWidth', 6, 'minSpacing', 2, 'maxModes', []);
if ~isfield(config, 'Wcorr') || isempty(config.Wcorr)
    return
end
if ~(isstruct(config.Wcorr) && isscalar(config.Wcorr))
    error('V3:InvalidWcorrOptions', ...
        'config.Wcorr must be a scalar struct when config.N="Wcorr".');
end

names = fieldnames(config.Wcorr);
supported = {'blockWidth', 'minSpacing', 'maxModes'};
unknown = setdiff(names, supported);
if ~isempty(unknown)
    error('V3:UnknownWcorrOption', ...
        'Unsupported config.Wcorr option: %s.', unknown{1});
end
for i = 1:numel(names)
    wcorr.(names{i}) = config.Wcorr.(names{i});
end

if ~(local_positive_integer_scalar(wcorr.blockWidth) && wcorr.blockWidth >= 2)
    error('V3:InvalidWcorrBlockWidth', ...
        'config.Wcorr.blockWidth must be an integer greater than or equal to 2.');
end
if ~local_positive_integer_scalar(wcorr.minSpacing)
    error('V3:InvalidWcorrMinSpacing', ...
        'config.Wcorr.minSpacing must be a positive integer.');
end
if ~(isempty(wcorr.maxModes) || local_positive_integer_scalar(wcorr.maxModes))
    error('V3:InvalidWcorrMaxModes', ...
        'config.Wcorr.maxModes must be empty or a positive integer.');
end

wcorr.blockWidth = double(wcorr.blockWidth);
wcorr.minSpacing = double(wcorr.minSpacing);
if ~isempty(wcorr.maxModes)
    wcorr.maxModes = double(wcorr.maxModes);
end
end

function tf = local_positive_integer_scalar(value)
tf = isnumeric(value) && isscalar(value) && isfinite(value) && ...
    value > 0 && value == fix(value);
end

function tf = local_positive_or_zero_integer_scalar(value)
tf = isnumeric(value) && isscalar(value) && isfinite(value) && ...
    value >= 0 && value == fix(value);
end

function tf = local_positive_or_signed_scalar(value)
tf = isnumeric(value) && isscalar(value) && isfinite(value);
end

function local_validate_threshold(value)
if ~(isscalar(value) && isfinite(value) && value > 0 && value < 1)
    error('V3:InvalidEigenvalueThreshold', ...
        'The eigenvalue threshold must be strictly between 0 and 1.');
end
end
