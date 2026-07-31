function setupInfo = setup_paths(ifilesRoot, dependencyRoot)
% SETUP_PATHS  Configure data directories and external software paths.
%
% Usage:
%   setup_paths(ifilesRoot)
%   setup_paths(ifilesRoot, dependencyRoot)
%
% Input:
%   ifilesRoot - char, path to the project root directory
%   dependencyRoot - optional folder containing downloaded dependencies

    arguments
        ifilesRoot (1,:) char
        dependencyRoot (1,:) char = ''
    end

    sourceRoot = stpc_source_root();
    setenv('IFILES', ifilesRoot);
    setenv('STPC_SRC', sourceRoot);

    directoryInfo = initialize_stpc_directories(ifilesRoot);
    dependencyRoot = stpc_dependency_root(dependencyRoot);
    setenv('STPC_DEPENDENCIES', dependencyRoot);

    requiredDependencies = {
        'slepian_alpha-master'
        'slepian_bravo-master'
        'slepian_delta-master'
        'slepian_zero-master'
        };
    missingDependencies = cell(0, 1);
    for index = 1:numel(requiredDependencies)
        dependencyPath = fullfile(dependencyRoot, ...
            requiredDependencies{index});
        if exist(dependencyPath, 'dir') == 7
            addpath(genpath(dependencyPath));
        else
            missingDependencies{end + 1, 1} = ...
                requiredDependencies{index}; %#ok<AGROW>
        end
    end

    if ~isempty(missingDependencies)
        error('STPCGA:MissingDependencies', ...
            ['Missing required Slepian dependencies under "%s": %s. ', ...
             'Pass their parent directory to startup_STPC_GA or set ', ...
             'the STPC_DEPENDENCIES environment variable.'], ...
            dependencyRoot, strjoin(missingDependencies, ', '));
    end

    ddkPath = fullfile(dependencyRoot, 'GRACE-filter-master');
    if exist(ddkPath, 'dir') == 7
        addpath(genpath(fullfile(ddkPath, 'src', 'matlab')));
    else
        warning('STPCGA:MissingDDKDependency', ...
            ['GRACE-filter-master was not found under "%s". ', ...
             'DDK comparison products will be unavailable.'], ...
            dependencyRoot);
    end

    % Keep STPC-GA compatibility functions (for example ls2cell.m) ahead
    % of same-named functions in the external Slepian toolboxes.
    addpath(genpath(sourceRoot), '-begin');

    setupInfo = struct( ...
        'projectRoot', ifilesRoot, ...
        'sourceRoot', sourceRoot, ...
        'dependencyRoot', dependencyRoot, ...
        'createdDirectories', {directoryInfo.createdDirectories}, ...
        'runtimeDirectories', {directoryInfo.allDirectories});
end
