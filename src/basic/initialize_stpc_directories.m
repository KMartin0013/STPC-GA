function directoryInfo = initialize_stpc_directories(ifilesRoot)
%INITIALIZE_STPC_DIRECTORIES Create writable STPC-GA runtime directories.
%
% Static input folders are not populated here. The function only creates
% directories in which STPC-GA and the Slepian toolboxes may save caches,
% generated models, and transformed coefficients.

arguments
    ifilesRoot (1,:) char
end

if exist(ifilesRoot, 'dir') ~= 7
    [ok, message] = mkdir(ifilesRoot);
    if ~ok
        error('STPCGA:ProjectRootCreateFailed', ...
            'Could not create the project root "%s": %s', ...
            ifilesRoot, message);
    end
end

relativeDirectories = {
    'COASTS'
    'GLMALPHA'
    'KERNELC'
    'KERNELCP'
    'LEGENDRE'
    'DLMB'
    'DLMLMP'
    'HASHES'
    'LOCALIZE'
    'WIGNER'
    'GRACE'
    fullfile('GRACE', 'Originals')
    fullfile('GRACE', 'Degree1')
    fullfile('GRACE', 'SlepianExpansions')
    'MOD'
    fullfile('MOD', 'SlepianExpansions')
    };

created = cell(0, 1);
available = cell(size(relativeDirectories));
for index = 1:numel(relativeDirectories)
    fullPath = fullfile(ifilesRoot, relativeDirectories{index});
    available{index} = fullPath;
    if exist(fullPath, 'dir') ~= 7
        [ok, message] = mkdir(fullPath);
        if ~ok
            error('STPCGA:DirectoryCreateFailed', ...
                'Could not create runtime directory "%s": %s', ...
                fullPath, message);
        end
        created{end + 1, 1} = fullPath; %#ok<AGROW>
    end
end

directoryInfo = struct( ...
    'allDirectories', {available}, ...
    'createdDirectories', {created});
end
