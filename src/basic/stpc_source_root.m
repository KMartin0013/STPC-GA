function sourceRoot = stpc_source_root()
% STPC_SOURCE_ROOT  Return the active STPC source directory.
%
% The environment variable allows tests and development copies to select a
% source tree explicitly. The public STPC-GA package normally uses /src.

sourceRoot = getenv('STPC_SRC');
if isempty(sourceRoot) || exist(sourceRoot, 'dir') ~= 7
    sourceRoot = fileparts(fileparts(mfilename('fullpath')));
    setenv('STPC_SRC', sourceRoot);
end
end
