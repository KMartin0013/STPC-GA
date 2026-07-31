function dependencyRoot = stpc_dependency_root(requestedRoot)
%STPC_DEPENDENCY_ROOT Locate the folder containing external dependencies.
%
% Search order:
%   1. Explicit requestedRoot
%   2. STPC_DEPENDENCIES environment variable
%   3. <projectRoot>/external
%   4. <sourceRoot>/required_softwares (legacy layout)

if nargin < 1
    requestedRoot = '';
end

sourceRoot = stpc_source_root();
projectRoot = fileparts(sourceRoot);
environmentRoot = getenv('STPC_DEPENDENCIES');

candidates = {
    char(string(requestedRoot))
    environmentRoot
    fullfile(projectRoot, 'external')
    fullfile(sourceRoot, 'required_softwares')
    };

dependencyRoot = '';
for index = 1:numel(candidates)
    candidate = candidates{index};
    if isempty(candidate)
        continue
    end
    if local_contains_dependency(candidate)
        dependencyRoot = candidate;
        break
    end
end

if isempty(dependencyRoot)
    nonemptyCandidates = candidates(~cellfun(@isempty, candidates));
    dependencyRoot = nonemptyCandidates{end};
end

dependencyRoot = char(string(dependencyRoot));
end

function tf = local_contains_dependency(candidate)
dependencyNames = {
    'slepian_alpha-master'
    'slepian_bravo-master'
    'slepian_delta-master'
    'slepian_zero-master'
    'GRACE-filter-master'
    };

tf = exist(candidate, 'dir') == 7 && any(cellfun( ...
    @(name) exist(fullfile(candidate, name), 'dir') == 7, ...
    dependencyNames));
end
