function suffix = stpc_note_suffix(note)
% STPC_NOTE_SUFFIX Return a filename-safe run-note suffix.

if nargin < 1 || isempty(note)
    suffix = '';
    return
end

suffix = char(string(note));
suffix = strtrim(suffix);
if isempty(suffix)
    return
end

suffix = regexprep(suffix, '[^\w-]', '_');
suffix = regexprep(suffix, '_+', '_');
if ~startsWith(suffix, '_')
    suffix = ['_' suffix];
end
end
