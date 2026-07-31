function config = stpc_configure_center_list(config)
% STPC_CONFIGURE_CENTER_LIST Normalize center-list fields in a run config.
%
% The manuscript-level center_list corresponds to config.givInstitu in code.
% config.use_institu is the center list plus the derived ensemble tag.

if ~isfield(config, 'givInstitu') || isempty(config.givInstitu)
    error('stpc_configure_center_list:MissingCenters', ...
        'config.givInstitu must list the input GRACE/GRACE-FO products.');
end

givInstitu = reshape(string(config.givInstitu), 1, []);
config.givInstitu = givInstitu;
config.center_list = givInstitu;

registryFile = local_registry_file(config);
[canonicalKey, sortedCenters] = local_canonical_key(givInstitu, config);
config.center_canonical_key = canonicalKey;
config.center_list_sorted = sortedCenters;
config.ensembleRegistryFile = registryFile;

hasExplicitEnsemble = isfield(config, 'ensInstitu') && ~isempty(config.ensInstitu);
if hasExplicitEnsemble
    ensInstitu = string(config.ensInstitu);
elseif isfield(config, 'use_institu') && ...
        numel(config.use_institu) == numel(givInstitu) + 1
    ensInstitu = string(config.use_institu(end));
else
    [ensInstitu, foundInRegistry] = local_lookup_registry(registryFile, canonicalKey);
    if ~foundInRegistry
        ensInstitu = local_default_ensemble_tag(sortedCenters);
    end
end

if numel(ensInstitu) ~= 1
    error('stpc_configure_center_list:InvalidEnsembleTag', ...
        'config.ensInstitu must be a scalar string or character vector.');
end

config.ensInstitu = ensInstitu;
config.use_institu = [givInstitu ensInstitu];

local_update_registry(config, hasExplicitEnsemble);
end

function ensInstitu = local_default_ensemble_tag(givInstitu)
centerLetters = strings(1, numel(givInstitu));
timeTags = strings(1, numel(givInstitu));

for i = 1:numel(givInstitu)
    product = char(givInstitu(i));
    token = regexp(product, '^(.+?)(\d{4})$', 'tokens', 'once');
    if isempty(token)
        centerName = product;
        timeTags(i) = "";
    else
        centerName = token{1};
        timeTags(i) = string(token{2});
    end

    if isempty(centerName)
        error('stpc_configure_center_list:InvalidCenterName', ...
            'Invalid product name: %s', product);
    end
    centerLetters(i) = string(centerName(1));
end

ensInstitu = join(centerLetters, "");
if all(timeTags ~= "") && numel(unique(timeTags)) == 1
    ensInstitu = ensInstitu + timeTags(1);
end
end

function registryFile = local_registry_file(config)
if isfield(config, 'ensembleRegistryFile') && ~isempty(config.ensembleRegistryFile)
    registryFile = char(config.ensembleRegistryFile);
elseif isfield(config, 'ifilesRoot') && ~isempty(config.ifilesRoot)
    registryFile = fullfile(char(config.ifilesRoot), ...
        'STPC_Ensemble_Registry.xlsx');
else
    registryFile = '';
end
end

function [canonicalKey, sortedCenters] = local_canonical_key(givInstitu, config)
sortedCenters = sort(reshape(string(givInstitu), 1, []));
releaseTag = "";
if isfield(config, 'Institu_ver') && ~isempty(config.Institu_ver)
    releaseTag = string(config.Institu_ver);
end
canonicalKey = char(releaseTag + ":" + join(sortedCenters, "|"));
end

function [ensInstitu, found] = local_lookup_registry(registryFile, canonicalKey)
ensInstitu = "";
found = false;
if isempty(registryFile) || exist(registryFile, 'file') ~= 2
    return
end

try
    registry = readtable(registryFile, 'TextType', 'string');
catch
    warning('stpc_configure_center_list:RegistryReadFailed', ...
        'Could not read ensemble registry: %s', registryFile);
    return
end

if ~all(ismember(["canonical_key", "ensInstitu"], string(registry.Properties.VariableNames)))
    return
end

matchIndex = find(registry.canonical_key == string(canonicalKey), 1, 'first');
if ~isempty(matchIndex)
    ensInstitu = string(registry.ensInstitu(matchIndex));
    found = true;
end
end

function local_update_registry(config, hasExplicitEnsemble)
registryFile = config.ensembleRegistryFile;
if isempty(registryFile)
    return
end

registry = table();
if exist(registryFile, 'file') == 2
    try
        registry = readtable(registryFile, 'TextType', 'string');
    catch
        warning('stpc_configure_center_list:RegistryReadFailed', ...
            'Could not read ensemble registry: %s', registryFile);
        return
    end
end

requiredNames = ["canonical_key", "center_list_input", "center_list_sorted", ...
    "release_tag", "ensInstitu", "first_created", "note", "result_dir", ...
    "comment"];
if isempty(registry)
    registry = array2table(strings(0, numel(requiredNames)), ...
        'VariableNames', cellstr(requiredNames));
elseif ~all(ismember(requiredNames, string(registry.Properties.VariableNames)))
    return
end

for name = requiredNames
    registry.(name) = string(registry.(name));
end

if any(string(registry.canonical_key) == string(config.center_canonical_key))
    return
end

releaseTag = "";
if isfield(config, 'Institu_ver') && ~isempty(config.Institu_ver)
    releaseTag = string(config.Institu_ver);
end
noteText = "";
if isfield(config, 'note') && ~isempty(config.note)
    noteText = string(config.note);
end
resultDir = "";
if isfield(config, 'resultDir') && ~isempty(config.resultDir)
    resultDir = string(config.resultDir);
end
commentText = "";
if hasExplicitEnsemble
    commentText = "explicit ensemble tag";
end

newRow = table( ...
    string(config.center_canonical_key), ...
    join(string(config.givInstitu), ","), ...
    join(string(config.center_list_sorted), ","), ...
    releaseTag, string(config.ensInstitu), ...
    string(datetime("now", "Format", "yyyy-MM-dd HH:mm:ss")), ...
    noteText, resultDir, commentText, ...
    'VariableNames', cellstr(requiredNames));

registry = [registry; newRow];
try
    registryDir = fileparts(registryFile);
    if ~isempty(registryDir) && exist(registryDir, 'dir') ~= 7
        mkdir(registryDir);
    end
    writetable(registry, registryFile);
catch
    warning('stpc_configure_center_list:RegistryWriteFailed', ...
        'Could not update ensemble registry: %s', registryFile);
end
end
