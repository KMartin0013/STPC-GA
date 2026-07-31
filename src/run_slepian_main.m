function slepian_results = run_slepian_main(basicInfo, slepianInfo, ...
    Dataproduct)
% RUN_SLEPIAN_MAIN  Compute Slepian coefficients, mass changes, and gridded EWH for a given basin
%
% Usage:
%   results = run_slepian_main(basicInfo, slepianInfo, paths, options)
%
% Inputs:
%   basicInfo  - struct, corresponds to fields in original Basic_Information.mat, e.g.:
%       .Code_Version
%       .TH_ori
%       .Area
%       .Lwindow
%       .land_or_ocean
%       .c11cmn
%       .ll
%       .ddir1         % output data directory
%       .ddir2         % output figure directory
%       .ifilesRoot    % root directory of original IFILES
%       .plotRegion    (true/false) whether to plot study region
%       .saveMainData  (true/false) whether to save MainData_*.mat

%   slepianInfo - struct, corresponds to fields in original Slepian_Information.mat, e.g.:
%       .Dataproduct
%       .buffer_deg
%       .buffer_str
%       .S_choice
%       .Max_S
%       .Radius
%       .phi, .theta, .omega
%       .artificial_months
%
%
% Output:
%   results    - struct, contains most variables that were saved in MainData_*.mat:
%       .Dataproduct
%       .use_sleptdelta, .use_signalcoffs, .use_residcoffs
%       .use_MASS, .use_MASSfit, .use_MASS_CC, ...
%       .fill_signaldelta, .fill_residdelta, .fill_deltacoffs
%       .fill_MASS, .fill_MASSsig, .fill_MASSres
%       .fill_Grid_EWH, .fill_Grid_EWHsig, .fill_Grid_EWHres
%       .CC, .V, .S, .S_sig, .S_shannon, .Radius, ...
%       and other variables needed for subsequent MSSA steps.
%
% Main Reference:
% Ref1: Harig, C., & Simons, F. J. (2012). Mapping GreenlandAs mass loss in space and time. Proceedings of the National Academy of Sciences, 109(49), 19934-19937.
% Ref2: Ma, Z., Fok, H. S., Tenzer, R., & Chen, J. (2024). A novel Slepian approach for determining mass-term sea level from GRACE over the South China Sea. International Journal of Applied Earth Observation and Geoinformation, 132, 104065.
% Ref3: Gauer, L. M., Chanard, K., & Fleitout, L. (2023). Data‐driven gap filling and spatio‐temporal filtering of the GRACE and GRACE‐FO records. Journal of Geophysical Research: Solid Earth, 128(5), e2022JB025561.
%
% Note:
%   In the original scripts, many intermediate variables were saved to mat files.
%   In this functional version, we only keep variables in 'results' that are
%   needed by later steps. Other intermediates can be saved to disk when needed.

    arguments
        basicInfo   struct
        slepianInfo struct
        Dataproduct cell
    end

    % ---------------------------------------------------------------
    % 0. Unpack basic fields
    % ---------------------------------------------------------------
    Code_Version = basicInfo.Code_Version; %#ok<NASGU>
    TH_ori       = basicInfo.TH_ori;
    Area         = basicInfo.Area;
    Lwindow      = basicInfo.Lwindow;
    land_or_ocean= basicInfo.land_or_ocean;
    c11cmn       = basicInfo.c11cmn;
    Max_S        = basicInfo.Max_S;
    ddir1        = basicInfo.ddir1;
    ddir2        = basicInfo.ddir2;
    ifilesRoot   = basicInfo.ifilesRoot;
    
    redo         = basicInfo.redo;
    plotProcess  = basicInfo.plotProcess;
    saveAddData  = basicInfo.saveAddData;

    buffer_deg        = slepianInfo.buffer_deg
    buffer_str        = slepianInfo.buffer_str;
    S_selection       = slepianInfo.S_selection;
    Radius            = slepianInfo.Radius;
    phi               = slepianInfo.phi;
    theta             = slepianInfo.theta;
    omega             = slepianInfo.omega;
    artificial_months = slepianInfo.artificial_months;

    % ---------------------------------------------------------------
    % 1. Build name tags, IB filename, region polygons, and choose S
    % ---------------------------------------------------------------
    [FIG_Attach, fnpl_IB] = build_name_tags_and_ib_path( ...
        Dataproduct, Lwindow, buffer_str, ifilesRoot);
    if isfield(basicInfo, 'note')
        FIG_Attach = [FIG_Attach stpc_note_suffix(basicInfo.note)];
    end

    [XY_buffer, BasinArea_buffer, XY_ori, BasinArea_ori, Earth_radius] = ...
        compute_regions(TH_ori, buffer_deg);

    maybe_plot_region(XY_buffer, XY_ori, FIG_Attach, buffer_deg, ...
        ddir2, plotProcess);

    if redo || ~exist(fullfile(ddir1, ['MainSlep_' FIG_Attach '.mat']), 'file')

        disp('I''m working, not stuck. If stuck, I will tell you.')

        % ---------------------------------------------------------------
        % 2. Project GRACE to Slepian basis & choose S
        % ---------------------------------------------------------------
        [data_slepcoffs, data_dates, CC, V, TH_buffer] = ...
            project_grace_to_slepian(Dataproduct, TH_ori, buffer_deg, ...
            Lwindow, phi, theta, omega, ifilesRoot);

        [S_shannon, S, S_sig] = ...
            choose_S_from_eigenvalues(S_selection, Max_S, V, ...
            Lwindow, XY_buffer);

        % ---------------------------------------------------------------
        % 3. Selection region polygons
        % ---------------------------------------------------------------

        [TH, XY, BasinArea] = select_integration_region( ...
            Dataproduct, TH_ori, TH_buffer, XY_ori, XY_buffer, ...
            BasinArea_buffer, BasinArea_ori);

        % ---------------------------------------------------------------
        % 4. Handle missing months and build time metadata
        % ---------------------------------------------------------------
        DateInfo = handle_missing_months(data_dates, Dataproduct);

        % ---------------------------------------------------------------
        % 5. Temporal fitting of Slepian coefficients
        % ---------------------------------------------------------------
        fitwhat = [3 365.25 365.25/2];  % trend + annual + semi-annual

        [timeRes, fillRes] = fit_slepian_timeseries( ...
            data_slepcoffs, CC, TH, S, S_sig, Radius, ...
            DateInfo, fitwhat);

        % ---------------------------------------------------------------
        % 6. Build gridded EWH and basin-averaged mass time series
        % ---------------------------------------------------------------
        gridRes = compute_grid_EWH_and_mass( ...
            CC, S, S_sig, Radius, Lwindow, c11cmn, ...
            XY, DateInfo, fillRes);

        % ---------------------------------------------------------------
        % 7. Handle artifical months for sensitivity analysis
        % ---------------------------------------------------------------
        if any(artificial_months)
            ArtInfo = handle_artificial_months(data_slepcoffs, DateInfo, ...
                artificial_months, CC, TH, S_sig, Radius, fitwhat);
        else
            ArtInfo = struct();
        end

        % Validation
        figure
        plot(1:DateInfo.fill_nmonths,fillRes.fill_MASS);
        hold on
        plot(1:DateInfo.fill_nmonths,gridRes.fill_awMASS);

        legend('Code','Area-weighted');

        % ---------------------------------------------------------------
        % 8. If ocean region, apply IB correction (MSL fields)
        % ---------------------------------------------------------------
        if strcmp(Dataproduct{4}, 'ocean')
            oceanRes = apply_IB_correction( ...
                fnpl_IB, fitwhat, DateInfo, gridRes, timeRes, fillRes, S);
        else
            oceanRes = struct();
        end

        % ---------------------------------------------------------------
        % 9. Collect outputs into a single struct
        % ---------------------------------------------------------------
        slepian_results = struct();
        slepian_results.Dataproduct = Dataproduct;
        slepian_results.TH          = TH;
        slepian_results.XY          = XY;
        slepian_results.XY_ori      = XY_ori;
        slepian_results.XY_buffer   = XY_buffer;
        slepian_results.BasinArea   = BasinArea;
        slepian_results.S           = S;
        slepian_results.S_sig       = S_sig;
        slepian_results.S_shannon   = S_shannon;
        slepian_results.S_selection = S_selection;
        slepian_results.CC          = CC;
        slepian_results.V           = V;
        slepian_results.Radius      = Radius;
        slepian_results.FIG_Attach  = FIG_Attach;

        slepian_results.time  = timeRes;
        slepian_results.fill  = fillRes;
        slepian_results.grid  = gridRes;
        slepian_results.ocean = oceanRes;
        slepian_results.date  = DateInfo;
        slepian_results.arti  = ArtInfo;
        % ---------------------------------------------------------------
        % 10. Save MainData_*.mat
        % ---------------------------------------------------------------

        matFile = fullfile(ddir1, ['MainSlep_' FIG_Attach '.mat']);
        save(matFile, '-struct', 'slepian_results');

    else

        disp(['SLEPIAN: Loading ',fullfile(ddir1, ['MainSlep_' FIG_Attach '.mat'])]);

        slepian_results=load(fullfile(ddir1, ['MainSlep_' FIG_Attach '.mat']));

    end

end

%% ------------------------------------------------------------------------
%% Subfunctions
%% ------------------------------------------------------------------------

function [FIG_Attach, fnpl_IB] = build_name_tags_and_ib_path( ...
    Dataproduct, Lwindow, buffer_str, ifilesRoot)
% BUILD_NAME_TAGS_AND_IB_PATH  Build figure/file name prefix and IB file path.

%     FIG_Attach = sprintf('%s_%s_%s_%s_%s_%s', ...
%         Dataproduct{1}, Dataproduct{2}(3:4), Area, ...
%         num2str(Lwindow), buffer_str, num2str(Radius));
    FIG_Attach = sprintf('%s_%s_%s_%s', ...
        Dataproduct{1}, Dataproduct{2}(3:4), ...
        num2str(Lwindow), buffer_str);

    fnpl_IB = fullfile(ifilesRoot, 'GRACE', ...
        sprintf('%s_%s_IB_%s_%s_%s.mat', Dataproduct{1}, Dataproduct{2}, ...
                num2str(Lwindow), 'SD', 'GSHHC'));
end

function [TH, XY, BasinArea] = select_integration_region( ...
    Dataproduct, TH_ori, TH_buffer, XY_ori, XY_buffer, ...
        BasinArea_buffer, BasinArea_ori)
% SELECT_INTEGRATION_REGION  Decide which polygon is used for integration.
%
% For land basins (except Greenland), use the original region (no buffer).
% For ocean basins, include buffer in the integration region.
% For regions with severe leakage-out effect (e.g., Greenland), also
% include buffer zones to maximize the restored signals

    if strcmp(Dataproduct{4},'land') && ~strcmp(TH_ori,'greenland')
        TH        = TH_ori;
        XY        = XY_ori;
        BasinArea = BasinArea_ori;
    else
        TH        = TH_buffer;    
        XY        = XY_buffer;
        BasinArea = BasinArea_buffer;
    end
end

function [data_slepcoffs, data_dates, CC, V, TH_buffer] = ...
    project_grace_to_slepian(Dataproduct, TH_ori, buffer_deg, ...
                             Lwindow, phi, theta, omega, ifilesRoot)
% PROJECT_GRACE_TO_SLEPIAN  Project GRACE fields onto Slepian basis.
%
% Outputs:
%   data_slepcoffs  - Slepian coefficients (kg/m^2)
%   data_dates      - datenum time vector
%   CC, V           - Slepian basis objects / eigenvalues
%   TH_buffer       - region handle used when building Slepian basis

    % Load Earth constants (kept for completeness)
    load(fullfile(ifilesRoot, 'EARTHMODELS', 'CONSTANTS', 'Earth.mat')); %#ok<LOAD>

    [data_slepcoffs, calerrors, data_dates, TH_buffer, G, CC, V] = ...
        grace2slept_m(Dataproduct, TH_ori, buffer_deg, Lwindow, ...
                      phi, theta, omega, [], 'SD', 1); %#ok<ASGLU>
end

function DateInfo = handle_missing_months(data_dates, Dataproduct)
% HANDLE_MISSING_MONTHS  Build time axes, detect gaps, and define dates.
%
% This function:
%   * Converts calendar dates to "month index" from the start of dataset
%   * Removes duplicate months
%   * Detects missing months (gaps)
%   * Builds Slept_dates (including gaps) and fill_dates

    InsYear_str = char(Dataproduct{1});
    year_start  = str2double(InsYear_str(end-3:end-2));
    year_end    = str2double(InsYear_str(end-1:end));
    data_year_beg = 2000 + year_start;
    data_year_end = 2000 + year_end;

    % Month indices for each data_dates entry
    data_months = (str2num(datestr(data_dates,'yyyy')) - data_year_beg) * 12 + ...
                  str2num(datestr(data_dates,'mm'));
    data_nmonths = numel(data_months);

    % Full reference months (no gap)
    fill_nmonths = data_months(end)-data_months(1)+1;
    fill_months  = (1:fill_nmonths)';

    % Remove duplicate months; keep first occurrence
    [use_months, data2use_monIndex] = unique(data_months, 'stable');
    use_nmonths = numel(use_months);
    use_dates   = data_dates(data2use_monIndex);

    % Gaps between GRACE/GRACE-FO
    missing_months = setdiff(fill_months, use_months)';
    missing_dates  = [];

    if any(missing_months)
        missing_dates = zeros(1, numel(missing_months));
        for i = 1:numel(missing_months)
            m0 = missing_months(i);
            monthstart = datenum([ ...
                data_year_beg + floor((m0-1)/12), ...
                mod(m0-1, 12) + 1, 1]);
            monthend = datenum([ ...
                data_year_beg + floor(m0/12), ...
                mod(m0, 12) + 1, 1]) - 1;
            missing_dates(i) = (monthstart + monthend)/2;
        end
        Slept_dates = [data_dates, missing_dates];
        Slept_months = [data_months; missing_months'];
        fill_dates   = zeros(1, numel(fill_months));
        fill_dates(use_months)     = use_dates;
        fill_dates(missing_months) = missing_dates;
    else
        Slept_dates  = data_dates;
        Slept_months = data_months;
        fill_dates   = data_dates;
    end

    DateInfo = struct();
    DateInfo.fill_nmonths     = fill_nmonths;
    DateInfo.fill_months      = fill_months;
    DateInfo.fill_dates       = fill_dates;
    DateInfo.use_months       = use_months;
    DateInfo.use_nmonths      = use_nmonths;
    DateInfo.use_dates        = use_dates;
    DateInfo.missing_months   = missing_months;
    DateInfo.missing_dates    = missing_dates;
    DateInfo.data_months      = data_months;
    DateInfo.data_dates       = data_dates;
    DateInfo.data_nmonths     = data_nmonths;
    DateInfo.data2use_monIndex = data2use_monIndex;
    DateInfo.data_year_beg    = data_year_beg;
    DateInfo.data_year_end    = data_year_end;
    DateInfo.Slept_months     = Slept_months;
    DateInfo.Slept_dates      = Slept_dates;

end

function [timeRes, fillRes, ArtiRes] = fit_slepian_timeseries( ...
    data_slepcoffs, CC, TH, S, S_sig, Radius, DateInfo, fitwhat)
% FIT_SLEPIAN_TIMESERIES  Fit temporal model to Slepian coefficients.
%
% This wraps slept2resid_m and handles:
%   * signal + residual decomposition for Slepian coefficients
%   * total mass time series from functionintegrals
%   * construction of filled series (signal, residual, delta)

    Slept_dates      = DateInfo.Slept_dates;
    data_dates       = DateInfo.Slept_dates(1:numel(DateInfo.data_months));
    data_months      = DateInfo.data_months;
    data_nmonths     = DateInfo.data_nmonths;
    use_months       = DateInfo.use_months;
    use_nmonths      = DateInfo.use_nmonths;
    use_dates        = DateInfo.use_dates;
    fill_months      = DateInfo.fill_months;
    fill_nmonths     = DateInfo.fill_nmonths;
    missing_dates    = DateInfo.missing_dates;
    missing_months   = DateInfo.missing_months;
    data2use_monIndex = DateInfo.data2use_monIndex;

    %----- Call slept2resid_m to get signal/residual and total mass -----
    [data_signalcoffs, data_residcoffs, ~, data_extravalues, ...
     data_MASS, ~, data_MASSparams, data_MASSparamerrors, ...
     data_MASSfit, functionintegrals, ~] = ...
        slept2resid_m(data_slepcoffs, Slept_dates, fitwhat, [], [], ...
                      CC, TH, S_sig, Radius);

    % Subset to unique months
    use_signalcoffs = data_signalcoffs(data2use_monIndex, :);
    use_residcoffs  = data_residcoffs(data2use_monIndex, :);
    use_MASS        = data_MASS(data2use_monIndex);
    use_MASSfit     = data_MASSfit(data2use_monIndex, :);

    % Covariance of residuals
    Cab = slepresid2cov(use_residcoffs);

    % Mean-referenced coefficients (delta) for unique months
    use_slepcoffs  = data_slepcoffs(data2use_monIndex, :);
    use_sleptdelta = use_slepcoffs(1:use_nmonths, :) - ...
        repmat(mean(data_slepcoffs(1:data_nmonths, :), 1), use_nmonths, 1);

    % Integrate to total mass (code version)
    [use_MASS_CC, ~, ~, ~, use_MASSfit_CC, ~] = ...
        integral_fit(S, use_dates, functionintegrals, use_sleptdelta, Cab);

    %----- Fill gaps using fitted values in data_extravalues -----
    if any(missing_months)
        fill_signalcoffs_com = zeros(fill_nmonths, size(use_signalcoffs, 2));
        fill_residcoffs_com  = zeros(fill_nmonths, size(use_residcoffs,  2));

        % Insert existing months
        for i = 1:use_nmonths
            fill_signalcoffs_com(use_months(i), :) = use_signalcoffs(i, :);
            fill_residcoffs_com(use_months(i), :)  = use_residcoffs(i, :);
        end
        % Insert gap estimates
        for i = 1:numel(missing_months)
            fill_signalcoffs_com(missing_months(i), :) = data_extravalues(i, :);
            fill_residcoffs_com(missing_months(i), :)  = 0;
        end

        fill_signaldelta = fill_signalcoffs_com(1:fill_nmonths, :) - ...
            repmat(mean(use_signalcoffs(1:use_nmonths, :), 1), fill_nmonths, 1);
        fill_residdelta = fill_residcoffs_com(1:fill_nmonths, :) - ...
            repmat(mean(use_residcoffs(1:use_nmonths, :), 1), fill_nmonths, 1);
    else
        fill_signaldelta = use_signalcoffs(1:use_nmonths, :) - ...
            repmat(mean(use_signalcoffs(1:use_nmonths, :), 1), use_nmonths, 1);
        fill_residdelta = use_residcoffs(1:use_nmonths, :) - ...
            repmat(mean(use_residcoffs(1:use_nmonths, :), 1), use_nmonths, 1);
    end

    fill_deltacoffs = fill_signaldelta + fill_residdelta;

    % Integrate filled series using function 'integratebasis'
    [fill_MASSsig_CC, ~, ~, ~, fill_MASSsigfit_CC, ~] = ...
        integral_fit(S, DateInfo.fill_dates, functionintegrals, fill_signaldelta, Cab);
    [fill_MASSres_CC, ~, ~, ~, fill_MASSresfit_CC, ~] = ...
        integral_fit(S, DateInfo.fill_dates, functionintegrals, fill_residdelta, Cab);

    fill_MASSsig = sum(fill_MASSsig_CC);
    fill_MASSres = sum(fill_MASSres_CC);
    fill_MASS    = fill_MASSsig + fill_MASSres;

    %----- Pack outputs -----
    timeRes = struct();
    timeRes.use_sleptdelta  = use_sleptdelta;
    timeRes.use_signalcoffs = use_signalcoffs;
    timeRes.use_residcoffs  = use_residcoffs;
    timeRes.use_MASS        = use_MASS;
    timeRes.use_MASSfit     = use_MASSfit;
    timeRes.use_MASS_CC     = use_MASS_CC;
    timeRes.use_MASSfit_CC  = use_MASSfit_CC;
    timeRes.functionintegrals = functionintegrals;
    timeRes.Cab             = Cab;
    timeRes.data_MASSparams = data_MASSparams;
    timeRes.data_MASSparamerrors = data_MASSparamerrors;

    fillRes = struct();
    fillRes.fill_signaldelta    = fill_signaldelta;
    fillRes.fill_residdelta     = fill_residdelta;
    fillRes.fill_deltacoffs     = fill_deltacoffs;
    fillRes.fill_MASS           = fill_MASS;
    fillRes.fill_MASSsig        = fill_MASSsig;
    fillRes.fill_MASSres        = fill_MASSres;
    fillRes.fill_MASSsig_CC     = fill_MASSsig_CC;
    fillRes.fill_MASSres_CC     = fill_MASSres_CC;
    fillRes.fill_MASSsigfit_CC  = fill_MASSsigfit_CC;
    fillRes.fill_MASSresfit_CC  = fill_MASSresfit_CC;

    ArtiRes = struct();

end

function ArtInfo = handle_artificial_months(data_slepcoffs, DateInfo, ...
    artificial_months, CC, TH, S_sig, Radius, fitwhat)
% HANDLE_MISSING_MONTHS  Build time axes, detect gaps, and define dates.
%
% This function:
%   * Converts calendar dates to "month index" from the start of dataset
%   * Removes duplicate months
%   * Detects missing months (gaps)
%   * Builds Slept_dates (including gaps) and fill_dates

    data_dates       = DateInfo.data_dates;
    data_months      = DateInfo.data_months;
    use_dates        = DateInfo.use_dates;
    use_months       = DateInfo.use_months;
    fill_dates       = DateInfo.fill_dates;
    fill_nmonths     = DateInfo.fill_nmonths;
    data_year_beg    = DateInfo.data_year_beg;
    missing_months   = DateInfo.missing_months;
    missing_dates    = DateInfo.missing_dates;

    Art_data_months=data_months;Art_data_dates=data_dates;
    Art_data_slepcoffs=data_slepcoffs;
    [ia,~] = ismember(Art_data_months,artificial_months);
    Art_data_months(ia)=[];Art_data_dates(ia)=[];
    Art_data_slepcoffs(ia,:)=[];

    Art_use_months=use_months;Art_use_dates=use_dates;
    [ia,~] = ismember(Art_use_months,artificial_months);
    Art_use_months(ia)=[];Art_use_dates(ia)=[];
    Art_use_nmonths = length(Art_use_months);

    % raw data month (may include repeat months)
    Art_data_months=(str2num(datestr(Art_data_dates,'yyyy'))-data_year_beg)*12+ ...
        str2num(datestr(Art_data_dates,'mm'));
    Art_data_nmonths=numel(Art_data_months);

    % find repeat month (if repeat, use the first values)
    % transfer from data_mon to use_mon (without repeat months)
    [Art_use_months,Art_data2use_monIndex] = unique(Art_data_months,'stable'); % IndexTimeA返回索引

    for i=1:numel(artificial_months)
        artificial_monthstart = datenum([data_year_beg+floor((artificial_months(i)-1)/12) ...
            mod(artificial_months(i)-1,12)+1 1]);
        artificial_monthend = datenum([data_year_beg+floor((artificial_months(i))/12) ...
            mod(artificial_months(i),12)+1 1])-1;
        artificial_monthmid = (artificial_monthstart+artificial_monthend)/2;
        artificial_dates(i) = artificial_monthmid;
    end
    Art_Slept_dates=[Art_data_dates missing_dates artificial_dates];


    [Art_data_signalcoffs,Art_data_residcoffs,Art_data_ftests,...
        Art_data_extravalues]=slept2resid_m(Art_data_slepcoffs,...
        Art_Slept_dates,fitwhat,[],[],CC,TH,S_sig,Radius);

    % try my best to be consistent with previous definitions that each month
    % should have one and only have one value.
    Art_ESTsignal=Art_data_signalcoffs(Art_data2use_monIndex,:);
    Art_ESTresid=Art_data_residcoffs(Art_data2use_monIndex,:);

    % complete estiamted signal with filled leakage values
    Art_fill_signalcoffs_com=zeros(length(fill_dates),size(Art_ESTsignal,2));
    Art_fill_residcoffs_com=zeros(length(fill_dates),size(Art_ESTresid,2));
    for i=1:size(Art_ESTsignal,1)
        Art_fill_signalcoffs_com(Art_use_months(i),:)=Art_ESTsignal(i,:);
        Art_fill_residcoffs_com(Art_use_months(i),:)=Art_ESTresid(i,:);
    end

    if any(missing_months)

        for i=1:size(missing_months,2)
            Art_fill_signalcoffs_com(missing_months(i),:)=Art_data_extravalues(i,:);
            Art_fill_residcoffs_com(missing_months(i),:)=0;
        end

        for i=1:size(artificial_months,2)
            Art_fill_signalcoffs_com(artificial_months(i),:)=Art_data_extravalues(numel(missing_months)+i,:);
            %             ESTresid_com(leakage(i),:)=extravalues_resid(i,:);
            Art_fill_residcoffs_com(artificial_months(i),:)=0;
        end

    else

        for i=1:size(artificial_months,2)
            Art_fill_signalcoffs_com(artificial_months(i),:)=Art_data_extravalues(i,:);
            %             ESTresid_com(leakage(i),:)=extravalues_resid(i,:);
            Art_fill_residcoffs_com(artificial_months(i),:)=0;
        end

    end

    Art_fill_signaldelta=Art_fill_signalcoffs_com(1:fill_nmonths,:) - repmat(mean(Art_ESTsignal(1:Art_use_nmonths,:),1),fill_nmonths,1);
    Art_fill_residdelta=Art_fill_residcoffs_com(1:fill_nmonths,:) - repmat(mean(Art_ESTresid(1:Art_use_nmonths,:),1),fill_nmonths,1);

    Art_fill_delta=Art_fill_signaldelta+Art_fill_residdelta; % unit kg/m^2

    ArtInfo = struct();
    ArtInfo.Art_fill_delta     = Art_fill_delta;
    ArtInfo.artificial_months  = artificial_months;

end

function gridRes = compute_grid_EWH_and_mass( ...
    CC, S, S_sig, Radius, Lwindow, c11cmn, ...
    XY, gapInfo, fillRes)
% COMPUTE_GRID_EWH_AND_MASS  Build gridded EWH and basin-averaged mass.
%
% This function:
%   * Converts Slepian coefficients (kg/m^2) to gridded EWH (cm)
%   * Computes area-weighted basin-average mass time series (Gt)

    fill_dates      = gapInfo.fill_dates;
    fill_nmonths    = gapInfo.fill_nmonths;
    fill_signaldelta = fillRes.fill_signaldelta;
    fill_residdelta  = fillRes.fill_residdelta;

    % Example grid from first Slepian function
    [r_example, lon, lat, Plm, degres] = plm2xyz(CC{1}, 1, c11cmn, Lwindow); %#ok<ASGLU>

    if any(lon < 0)
        lon = lon + 360;
    end

    % Area weighting for each grid cell
    [lonlon, latlat] = meshgrid(lon, lat);
    [in, on] = check_polygon_in(XY, lonlon, latlat);
    c11cmn_area = zeros(numel(lat), numel(lon));
    for i = 1:numel(lat)
        for j = 1:numel(lon)
            c11cmn_area(i,j) = areaquad(lat(i)-0.5, lon(j)-0.5, ...
                                        lat(i)+0.5, lon(j)+0.5);
        end
    end

    % Precompute kernels for each Slepian function
    if numel(S_sig) > 1
        S1 = S_sig(1);
        S2 = S_sig(2);
        r_record = zeros(S2, size(r_example,1), size(r_example,2));
        for j = 1:S1
            [r_record(j,:,:), lon, lat, ~, ~] = ...
                plm2xyz(CC{j}, 1, c11cmn, Lwindow); 
        end
        for j = S1+1:S2
            CC_smo = plm_Gausssmooth(CC{j}, Radius);
            [r_record(j,:,:), lon, lat, ~, ~] = ...
                plm2xyz(CC_smo, 1, c11cmn, Lwindow); 
        end
    else
        r_record = zeros(S, size(r_example,1), size(r_example,2));
        for j = 1:S
            [r_record(j,:,:), lon, lat, ~, ~] = ...
                plm2xyz(CC{j}, 1, c11cmn, Lwindow); 
        end
    end

    fill_Grid_EWHsig = zeros(fill_nmonths, size(r_example,1), size(r_example,2));
    fill_Grid_EWHres = zeros(fill_nmonths, size(r_example,1), size(r_example,2));
    fill_awMASS      = zeros(1, fill_nmonths);

    for i = 1:fill_nmonths
        sp_ewh_signal   = 0;
        sp_ewh_residual = 0;
        for j = 1:S
            r = squeeze(r_record(j,:,:));
            sp_ewh_signal   = sp_ewh_signal   + r * fill_signaldelta(i,j)' / 1000 * 100;
            sp_ewh_residual = sp_ewh_residual + r * fill_residdelta(i,j)' / 1000 * 100;
        end
        fill_Grid_EWHsig(i,:,:) = sp_ewh_signal;
        fill_Grid_EWHres(i,:,:) = sp_ewh_residual;

        fill_awMASS(i) = sum( ...
            (sp_ewh_signal(in) + sp_ewh_residual(in)) / 100 .* c11cmn_area(in) ...
            ) * 4 * pi * 6370000^2 * 1e3 / 1e3 / 1e9;
    end

    fill_Grid_EWH = fill_Grid_EWHsig + fill_Grid_EWHres;

    gridRes = struct();
    gridRes.lon = lon;
    gridRes.lat = lat;
    gridRes.fill_Grid_EWH    = fill_Grid_EWH;
    gridRes.fill_Grid_EWHsig = fill_Grid_EWHsig;
    gridRes.fill_Grid_EWHres = fill_Grid_EWHres;
    gridRes.fill_awMASS      = fill_awMASS;
    gridRes.r_record         = r_record;
end

function oceanRes = apply_IB_correction( ...
    fnpl_IB, fitwhat, DateInfo, gridRes, timeRes, fillRes, S)
% APPLY_IB_CORRECTION  Apply inverted barometer (IB) correction for ocean basins.
%
% This function:
%   * Fits IB time series with the same model
%   * Constructs IB-corrected MSL fields on a grid (cm of salty water)
%   * Returns corrected fields and IB-related time series in oceanRes

    fill_nmonths    = DateInfo.fill_nmonths;
    use_nmonths     = DateInfo.use_nmonths;
    data2use_monIndex = DateInfo.data2use_monIndex;
    missing_months  = DateInfo.missing_months;

    use_sleptdelta   = timeRes.use_sleptdelta;
    fill_signaldelta = fillRes.fill_signaldelta;
    fill_residdelta  = fillRes.fill_residdelta;

%     fill_signaldelta = gridRes.fill_Grid_EWHsig;
%     fill_residdelta  = gridRes.fill_Grid_EWHres;

    r_record = gridRes.r_record;

    % Load IB time series
    load(fnpl_IB, 'IB_C', 'data_dates');  % IB_C(:,2): IB in mmH2O or kg/m^2
%     save(fnpl_IB, 'IB_C', 'data_dates', 'use_months');  % keep updated metadata

    data_IB = smooth(IB_C(:,2), 1);  % no smoothing here (can be increased)
    % unit kg/m^2 or mmH2O

    % Fit IB with same temporal model
    Slept_dates = DateInfo.Slept_dates;
    [data_IBsignal, data_IBresid, ~, data_IB_extravalues] = ...
        slept2resid_m(data_IB, Slept_dates, fitwhat);

    use_IBsignal = data_IBsignal(data2use_monIndex, :);
    use_IBresid  = data_IBresid(data2use_monIndex, :);

    % Fill gaps in IB
    if any(missing_months)
        fill_IBsignal_com = zeros(fill_nmonths, size(use_IBsignal, 2));
        fill_IBresid_com  = zeros(fill_nmonths, size(use_IBresid,  2));
        for i = 1:use_nmonths
            fill_IBsignal_com(DateInfo.use_months(i), :) = use_IBsignal(i, :);
            fill_IBresid_com(DateInfo.use_months(i), :)  = use_IBresid(i, :);
        end
        for i = 1:numel(missing_months)
            fill_IBsignal_com(missing_months(i), :) = data_IB_extravalues(i, :);
            fill_IBresid_com(missing_months(i), :)  = 0;
        end
        fill_IBsignaldelta = fill_IBsignal_com(1:fill_nmonths, :) - ...
            repmat(mean(use_IBsignal(1:use_nmonths, :), 1), fill_nmonths, 1);
        fill_IBresiddelta = fill_IBresid_com(1:fill_nmonths, :) - ...
            repmat(mean(use_IBresid(1:use_nmonths, :), 1), fill_nmonths, 1);
    else
        fill_IBsignaldelta = use_IBsignal(1:use_nmonths, :) - ...
            repmat(mean(use_IBsignal(1:use_nmonths, :), 1), use_nmonths, 1);
        fill_IBresiddelta = use_IBresid(1:use_nmonths, :) - ...
            repmat(mean(use_IBresid(1:use_nmonths, :), 1), use_nmonths, 1);
    end

    fill_IBdeltacoffs = fill_IBsignaldelta + fill_IBresiddelta;

    % Build MSL grids (cm of salty water) from OBP (kg/m^2) and IB series
    r_example = squeeze(r_record(1,:,:));
    fill_Grid_MSLsig = zeros(fill_nmonths, size(r_example,1), size(r_example,2));
    fill_Grid_MSLres = zeros(fill_nmonths, size(r_example,1), size(r_example,2));

    for i=1:fill_nmonths
        sp_obp_signal=0;sp_obp_residual=0;
        for j=1:S
            r=squeeze(r_record(j,:,:));

            sp_obp_signal=sp_obp_signal+r*fill_signaldelta(i,j)' ... % (you can manually adjust this period)
                /1000*100; % change from kg/m^2 to cm (/ 1000 kg/m*3 * 100)
            sp_obp_residual=sp_obp_residual+r*fill_residdelta(i,j)' ... % (you can manually adjust this period)
                /1000*100; % change from kg/m^2 to cm (/ 1000 kg/m*3 * 100)
        end
        % substract IB
        fill_Grid_MSLsig(i,:,:)=sp_obp_signal*1000/1028-fill_IBsignaldelta(i)'*1000/1028/10; % change unit from mmH2O to cm(salty H2O)
        fill_Grid_MSLres(i,:,:)=sp_obp_residual*1000/1028-fill_IBresiddelta(i)'*1000/1028/10; % change unit from mmH2O to cm(salty H2O)
    end

    fill_Grid_MSL=fill_Grid_MSLsig+fill_Grid_MSLres;

    % signal + residual (raw data)
    use_Grid_MSL=zeros(use_nmonths,size(r_example,1),size(r_example,2));
    for i=1:use_nmonths
        sp_obp=0;
        for j=1:S
            r=squeeze(r_record(j,:,:));

            sp_obp=sp_obp+r*use_sleptdelta(i,j)' ... % (you can manually adjust this period)
                /1000*100; % change from kg/m^2 to cm (/ 1028 kg/m*3 * 100)
        end
        % add IB
        use_Grid_MSL(i,:,:)=sp_obp*1000/1028-data_IB(data2use_monIndex(i))'*1000/1028/10; % change unit from mmH2O to cm(salty H2O)
    end

    oceanRes = struct();
    oceanRes.fill_Grid_MSLsig   = fill_Grid_MSLsig;
    oceanRes.fill_Grid_MSLres   = fill_Grid_MSLres;
    oceanRes.fill_Grid_MSL      = fill_Grid_MSL;
    oceanRes.fill_IBsignaldelta = fill_IBsignaldelta;
    oceanRes.fill_IBresiddelta  = fill_IBresiddelta;
    oceanRes.fill_IBdeltacoffs  = fill_IBdeltacoffs;
    oceanRes.use_Grid_MSL       = use_Grid_MSL;
    oceanRes.data_IB            = data_IB;
end
