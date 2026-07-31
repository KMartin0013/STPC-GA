function stpc_benchmark_append(logFile, entry)
% STPC_BENCHMARK_APPEND  Append one structured benchmark record.
%
% The CSV/MAT files are append-only across repeated runs, so first-run
% recomputation and later cache reads can be compared directly.

fields = {'run_id','timestamp','code_version','area','stage','item', ...
    'cache_status','elapsed_seconds','input_products','date_start', ...
    'date_end','output_file','note'};

for k = 1:numel(fields)
    if ~isfield(entry, fields{k}) || isempty(entry.(fields{k}))
        entry.(fields{k}) = '';
    end
end

if isnumeric(entry.elapsed_seconds)
    elapsedSeconds = entry.elapsed_seconds;
else
    elapsedSeconds = str2double(entry.elapsed_seconds);
end

T = table(string(entry.run_id), string(entry.timestamp), ...
    string(entry.code_version), string(entry.area), string(entry.stage), ...
    string(entry.item), string(entry.cache_status), elapsedSeconds, ...
    string(entry.input_products), string(entry.date_start), ...
    string(entry.date_end), string(entry.output_file), string(entry.note), ...
    'VariableNames', fields);

logDir = fileparts(logFile);
if ~exist(logDir, 'dir')
    mkdir(logDir);
end

if exist(logFile, 'file') == 2
    writetable(T, logFile, 'WriteMode', 'append');
else
    writetable(T, logFile);
end

matFile = strrep(logFile, '.csv', '.mat');
if exist(matFile, 'file') == 2
    S = load(matFile, 'benchmarkLog');
    benchmarkLog = [S.benchmarkLog; T]; %#ok<AGROW>
else
    benchmarkLog = T;
end
save(matFile, 'benchmarkLog');
end
