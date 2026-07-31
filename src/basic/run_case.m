function STPC_results = run_case(config)
% RUN_CASE  Main controller for the Slepian -> MSSA -> STPC workflow.
%
% Benchmark records are appended to Results_*/Benchmark_Log.csv and .mat.
% Stages are separated into Slepian, MSSA_smooth_extra, MSSA_main, and STPC
% so repeated runs can distinguish recomputation from cache reads.

    arguments
        config struct
    end

    config = v3_normalize_config(config);
    config = stpc_configure_center_list(config);
    if ~isfield(config, 'note')
        config.note = '';
    end
    config.note = stpc_note_suffix(config.note);

    %---------------------------%
    % 1. Paths and output dirs  %
    %---------------------------%
    setup_paths(config.ifilesRoot);

    if ~isfield(config, 'AddresultDir') || isempty(config.AddresultDir)
        config.AddresultDir = fullfile(config.resultDir, 'AddData');
    end

    if ~exist(config.resultDir, 'dir')
        mkdir(config.resultDir);
    end
    if ~exist(config.figureDir, 'dir')
        mkdir(config.figureDir);
    end
    if ~exist(config.AddresultDir, 'dir')
        mkdir(config.AddresultDir);
    end
    v3_validate_configuration(config);

    codeVersion = [config.codeVersionBase config.TH_ori];
    runId = datestr(now, 'yyyymmdd_HHMMSSFFF');
    logFile = fullfile(config.resultDir, 'Benchmark_Log.csv');
    inputProducts = strjoin(cellstr(config.givInstitu), ',');

    %---------------------------%
    % 2. Save Basic information %
    %---------------------------%
    infoSuffix = config.note;

    basicInfo = prepare_basic_info(codeVersion, config);
    save(fullfile(config.resultDir, ['Basic_Information' infoSuffix '.mat']), ...
        '-struct', 'basicInfo');

    %---------------------------------------------%
    % 3. Prepare and run Slepian info procedure   %
    %---------------------------------------------%
    num_givIns = numel(config.givInstitu);

    slepianInfo = prepare_slepian_info(config);

    if strcmp(slepianInfo.buffer_deg,'Auto')
        fprintf(['==== Running Slepian procedure (det_opt_buffer) to determine buffer zone. ====\n']);
        bufferTimer = tic;
        buffer_deg = det_opt_buffer(basicInfo, slepianInfo);
        bufferElapsed = toc(bufferTimer);
        slepianInfo.buffer_deg = buffer_deg;
        fprintf('==== End (buffer zone: %s) ====\n', num2str(slepianInfo.buffer_deg));
        stpc_log_stage(logFile, runId, basicInfo, 'Slepian', 'det_opt_buffer', ...
            'recomputed_or_checked', bufferElapsed, inputProducts, '', '', '', ...
            'Automatic buffer-zone selection.');
    end

    if slepianInfo.buffer_deg >= 0
        buffer_str = num2str(slepianInfo.buffer_deg);
        buffer_str(buffer_str=='.') = 'p';
    else
        buffer_str = ['neg' num2str(abs(slepianInfo.buffer_deg))];
        buffer_str(buffer_str=='.') = 'p';
    end

    slepianInfo.buffer_str = buffer_str;
    save(fullfile(config.resultDir, ['Slepian_Information' infoSuffix '.mat']), ...
        '-struct', 'slepianInfo');

    slepian_results = [];
    for Ins = 1:num_givIns
        productName = char(config.givInstitu(Ins));
        dataProduct = {productName, config.Institu_ver, ...
            config.Lwindow, config.landOrOcean, config.GIA};

        figAttach = sprintf('%s_%s_%s_%s', productName, ...
            config.Institu_ver(3:4), num2str(config.Lwindow), buffer_str);
        figAttach = [figAttach basicInfo.note];
        slepianCache = fullfile(config.resultDir, ['MainSlep_' figAttach '.mat']);
        cacheStatus = stpc_cache_status(config.redo, slepianCache);

        fprintf(['==== Running Slepian procedure (run_slepian_main) for ''%s'' ====\n'], productName);
        stageTimer = tic;
        slepian_result = run_slepian_main(basicInfo, slepianInfo, dataProduct);
        elapsed = toc(stageTimer);

        [dateStart, dateEnd] = stpc_result_date_range(slepian_result);
        stpc_log_stage(logFile, runId, basicInfo, 'Slepian', productName, ...
            cacheStatus, elapsed, productName, dateStart, dateEnd, slepianCache, '');

        disp(['SLEPIAN: The time series ''Code'' and ''Area-weighted'' should be ', ...
            'generally consistent. Otherwise, check c11cmn & buffer zones.']);
        fprintf(['==== End for ''%s'' ====\n'], productName)

        if Ins == 1
            slepian_results = repmat(slepian_result, 1, num_givIns);
        else
            slepian_results(Ins) = slepian_result;
        end
    end
    stpc_close_figures('Slepian', []);

    %---------------------------%
    % 5. MSSA settings & run    %
    %---------------------------%
    mssaInfo = prepare_mssa_info(config);
    save(fullfile(config.resultDir, ['MSSA_Information' infoSuffix '.mat']), ...
        '-struct', 'mssaInfo');

    Attach = build_mssa_paths_and_tags(basicInfo, slepianInfo.buffer_str, mssaInfo);

    smoothCache = fullfile(config.resultDir, ['MainMSSA_SMO2_' Attach.Attach_ALL '.mat']);
    cacheStatus = stpc_cache_status(config.redo, smoothCache);
    stageTimer = tic;
    mssa_smooth_results = run_mssa_smooth(basicInfo, slepianInfo, mssaInfo, ...
        slepian_results);
    elapsed = toc(stageTimer);
    [dateStart, dateEnd] = stpc_result_date_range(slepian_results(1));
    stpc_log_stage(logFile, runId, basicInfo, 'MSSA_smooth_extra', ...
        char(strjoin(cellstr(config.Smooth), ',')), cacheStatus, elapsed, ...
        inputProducts, dateStart, dateEnd, smoothCache, ...
        'Extra smoothing comparison; excluded from main MSSA runtime.');
    stpc_close_figures('MSSA_smooth_extra', []);

    mainMssaCache = fullfile(config.resultDir, ['MainMSSA_' Attach.Attach_ALL '.mat']);
    cacheStatus = stpc_cache_status(config.redo, mainMssaCache);
    stageTimer = tic;
    mssa_results = run_mssa_main(basicInfo, slepianInfo, mssaInfo, ...
        slepian_results);
    elapsed = toc(stageTimer);
    stpc_log_stage(logFile, runId, basicInfo, 'MSSA_main', char(Attach.Attach_ALL), ...
        cacheStatus, elapsed, inputProducts, dateStart, dateEnd, mainMssaCache, ...
        'Main MSSA gap filling and reconstruction.');
    stpc_close_figures('MSSA_main', []);

    %---------------------------%
    % 6. STPC filtering         %
    %---------------------------%
    STPCInfo = prepare_stpc_info(config);
    save(fullfile(config.resultDir, ['STPC_Information' infoSuffix '.mat']), ...
        '-struct', 'STPCInfo');

    stpcCache = fullfile(config.resultDir, ...
        ['MainSTPC_V3_' Attach.Attach_ALL '_datasets.mat']);
    cacheStatus = stpc_cache_status(config.redo, stpcCache);
    stpcFiguresBefore = findall(0, 'Type', 'figure');
    stageTimer = tic;
    STPC_p_results = run_stpc_main(basicInfo, STPCInfo, slepian_results, ...
        mssa_results);
    elapsed = toc(stageTimer);
    stpc_log_stage(logFile, runId, basicInfo, 'STPC', char(Attach.Attach_ALL), ...
        cacheStatus, elapsed, inputProducts, dateStart, dateEnd, stpcCache, ...
        'Main STPC filtering runtime.');
    stpc_close_figures('STPC', stpc_last_new_figure(stpcFiguresBefore));

    STPC_results = stpc_concat_results( ...
        STPC_p_results, mssa_smooth_results);
end

function cacheStatus = stpc_cache_status(redo, cacheFile)
    if ~redo && exist(cacheFile, 'file') == 2
        cacheStatus = 'cache_hit';
    else
        cacheStatus = 'recomputed_or_created';
    end
end

function [dateStart, dateEnd] = stpc_result_date_range(result)
    dateStart = '';
    dateEnd = '';
    if isfield(result, 'date')
        if isfield(result.date, 'fill_dates') && ~isempty(result.date.fill_dates)
            dateStart = datestr(min(result.date.fill_dates), 'yyyy-mm-dd');
            dateEnd = datestr(max(result.date.fill_dates), 'yyyy-mm-dd');
        elseif isfield(result.date, 'data_dates') && ~isempty(result.date.data_dates)
            dateStart = datestr(min(result.date.data_dates), 'yyyy-mm-dd');
            dateEnd = datestr(max(result.date.data_dates), 'yyyy-mm-dd');
        end
    end
end

function stpc_log_stage(logFile, runId, basicInfo, stage, item, cacheStatus, ...
        elapsed, inputProducts, dateStart, dateEnd, outputFile, note)
    entry = struct();
    entry.run_id = runId;
    entry.timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS.FFF');
    entry.code_version = basicInfo.Code_Version;
    entry.area = basicInfo.Area;
    entry.stage = stage;
    entry.item = item;
    entry.cache_status = cacheStatus;
    entry.elapsed_seconds = elapsed;
    entry.input_products = inputProducts;
    entry.date_start = dateStart;
    entry.date_end = dateEnd;
    entry.output_file = outputFile;
    entry.note = stpc_stage_note(basicInfo, note);
    stpc_benchmark_append(logFile, entry);
end

function noteText = stpc_stage_note(basicInfo, stageNote)
    noteText = stageNote;
    if isfield(basicInfo, 'note') && ~isempty(basicInfo.note)
        runNote = ['run_note=' basicInfo.note];
        if isempty(noteText)
            noteText = runNote;
        else
            noteText = [runNote '; ' noteText];
        end
    end
end

function keepFigure = stpc_last_new_figure(beforeFigures)
    keepFigure = [];
    figs = findall(0, 'Type', 'figure');
    if isempty(figs)
        return
    end

    newFigs = figs(~ismember(figs, beforeFigures));
    if isempty(newFigs)
        return
    end

    figNums = arrayfun(@(h) get(h, 'Number'), newFigs);
    [~, order] = sort(figNums, 'descend');
    keepFigure = newFigs(order(1));
end

function stpc_close_figures(stageName, keepFigures)
    if nargin < 2
        keepFigures = [];
    end

    figs = findall(0, 'Type', 'figure');
    if isempty(figs)
        fprintf('FIGURE_CLEANUP %s: no open figures.\n', stageName);
        return
    end

    keepFigures = keepFigures(ishandle(keepFigures));
    closeFigures = figs(~ismember(figs, keepFigures));

    if ~isempty(closeFigures)
        close(closeFigures);
    end

    fprintf('FIGURE_CLEANUP %s: closed %d figure(s), kept %d figure(s).\n', ...
        stageName, numel(closeFigures), numel(keepFigures));
    drawnow;
end
