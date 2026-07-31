function setupInfo = startup_STPC_GA(projectRoot, dependencyRoot)
%STARTUP_STPC_GA Configure STPC-GA (Version 1) for the current MATLAB session.
%
% setupInfo = startup_STPC_GA(projectRoot)
% setupInfo = startup_STPC_GA(projectRoot, dependencyRoot)
%
% projectRoot is the STPC-GA folder that contains src, GRACE, COASTS,
% Data, and the other project data directories. dependencyRoot is optional
% and should contain the downloaded folders:
%   slepian_alpha-master, slepian_bravo-master, slepian_delta-master,
%   slepian_zero-master, and optionally GRACE-filter-master.
%
% If dependencyRoot is omitted, STPC-GA checks the STPC_DEPENDENCIES
% environment variable, <projectRoot>/external, and
% <projectRoot>/src/required_softwares, in that order.

if nargin < 1 || isempty(projectRoot)
    projectRoot = fileparts(mfilename('fullpath'));
end
if nargin < 2
    dependencyRoot = '';
end

projectRoot = char(string(projectRoot));
dependencyRoot = char(string(dependencyRoot));
sourceRoot = fullfile(projectRoot, 'src');

if exist(sourceRoot, 'dir') ~= 7
    error('STPCGA:SourceNotFound', ...
        'The STPC-GA source directory does not exist: %s', sourceRoot);
end

addpath(genpath(sourceRoot));
setenv('STPC_SRC', sourceRoot);
setenv('IFILES', projectRoot);

setupInfo = setup_paths(projectRoot, dependencyRoot);

if nargout == 0
    fprintf('STPC-GA (Version 1) is ready.\n');
    fprintf('  Project root:    %s\n', setupInfo.projectRoot);
    fprintf('  Source root:     %s\n', setupInfo.sourceRoot);
    fprintf('  Dependency root: %s\n', setupInfo.dependencyRoot);
    if ~isempty(setupInfo.createdDirectories)
        fprintf('  Created %d runtime director%s.\n', ...
            numel(setupInfo.createdDirectories), ...
            local_plural_suffix(numel(setupInfo.createdDirectories)));
    end
end
end

function suffix = local_plural_suffix(count)
if count == 1
    suffix = 'y';
else
    suffix = 'ies';
end
end
