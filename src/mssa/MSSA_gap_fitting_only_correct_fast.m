function varargout = MSSA_gap_fitting_only_correct_fast(A5_CJGI_Sort, MSLA1_CJGI, M, N, opts, fig_note)

if nargin < 5 || isempty(opts)
    opts = struct();
end
if nargin >= 6 && ~isempty(fig_note)
    opts.fig = true;
    opts.fig_note = fig_note;
end

if ~isfield(opts, 'maxIter'), opts.maxIter = 100; end
if ~isfield(opts, 'chiThreshold'), opts.chiThreshold = 0.1; end
if ~isfield(opts, 'doTest'), opts.doTest = false; end
if ~isfield(opts, 'verbose'), opts.verbose = false; end
if ~isfield(opts, 'useSvds'), opts.useSvds = true; end
if ~isfield(opts, 'fig'), opts.fig = false; end
if ~isfield(opts, 'fig_note'), opts.fig_note = 'MSSA_gap_fitting_fast'; end
if ~isfield(opts, 'returnFullEvalues'), opts.returnFullEvalues = false; end

fig = opts.fig;
fig_note = opts.fig_note;

persistent pathAdded
if isempty(pathAdded)
    addpath(fileparts(mfilename('fullpath')));
    pathAdded = true;
end

Method = 2;
intit_num = numel(A5_CJGI_Sort);
time_num = size(MSLA1_CJGI, 1);
fill_months = (1:time_num).';

if isfield(A5_CJGI_Sort, 'data_year_beg')
    data_year_beg = A5_CJGI_Sort(1).data_year_beg;
else
    data_year_beg = 0;
end

% Plot settings
font_Size = 18;
F_gap = [.07 .03];
F_marg_h = [.09 .04];
F_marg_w = [.04 .02];
F_Position = [100,50,1500,900];

col_gap = [0 205 205] / 255;
col_est = [153 204 255] / 255;
col_out = [128 0 128] / 255;
abc = 'abcdefghijklmnopqrstuvwxyz';

tt_fil = floor(fill_months/12) + data_year_beg + mod(fill_months,12)/12 - 1/24;

% ===== first MSSA =====
[~,~,~,MSSA1_RC] = MSSA_fast(fill_months, MSLA1_CJGI, M, Method, N, opts);

MSSA1_CJGIave = sum(sum(MSSA1_RC, 3), 2) / intit_num;
std_MSSA1 = std(detrend(MSSA1_CJGIave));

Total_CGJI_fill_reconst_MSSA1 = MSLA1_CJGI;

if fig
    fl2 = figure;
    ha = tight_subplot(intit_num, 1, F_gap, F_marg_h, F_marg_w);
end

for i = 1:intit_num

    MSSA_TS = MSLA1_CJGI(:,i);
    missing_months = A5_CJGI_Sort(i).missing_months(:);

    outlier_months = find( ...
        MSSA_TS < MSSA1_CJGIave - 3 * std_MSSA1 | ...
        MSSA_TS > MSSA1_CJGIave + 3 * std_MSSA1);

    Total_fill_outfree = MSSA_TS;
    Total_fill_outfree(outlier_months) = MSSA1_CJGIave(outlier_months);
    Total_CGJI_fill_reconst_MSSA1(:,i) = Total_fill_outfree;

    if fig
        ybou1 = [-20, 20];

        if isfield(A5_CJGI_Sort, 'use_months')
            use_months = A5_CJGI_Sort(i).use_months(:);
        else
            use_months = fill_months;
        end

        tt_use = floor(use_months/12) + data_year_beg + mod(use_months,12)/12 - 1/24;
        tt_lea = floor(missing_months/12) + data_year_beg + mod(missing_months,12)/12 - 1/24;

        axes(ha(i));

        ybou1(2) = max(ybou1(2), max(MSSA_TS));
        ybou1(1) = min(ybou1(1), min(MSSA_TS));

        if isempty(tt_use)
            Cap_posix1 = tt_fil(1);
        else
            Cap_posix1 = tt_use(min(5, numel(tt_use)));
        end
        Cap_posiy1 = ybou1(2) - (ybou1(2)-ybou1(1))/10;

        baseline = ybou1(2);

        he_std = area(tt_fil, [MSSA1_CJGIave - 3*std_MSSA1, 6*std_MSSA1*ones(numel(tt_fil),1)], baseline);
        he_std(1).FaceColor = [1 1 1];
        he_std(1).EdgeColor = col_est;
        he_std(2).FaceColor = col_est;
        he_std(2).EdgeColor = col_est;
        he_std(2).FaceAlpha = 0.4;
        he_std(2).EdgeAlpha = 0.7;

        hold on

        p1 = plot(tt_fil, MSSA_TS, 'color', 'black', 'linewidth', 2);

        if ~isempty(outlier_months)
            tt_out = tt_fil(outlier_months);

            for j = 1:numel(tt_out)
                if outlier_months(j) == 1
                    outlier_tt_neibour = [tt_out(j), tt_out(j)+1/24];
                elseif outlier_months(j) == numel(tt_fil)
                    outlier_tt_neibour = [tt_out(j)-1/24, tt_out(j)];
                else
                    outlier_tt_neibour = [tt_out(j)-1/24, tt_out(j), tt_out(j)+1/24];
                end

                he_out = area(outlier_tt_neibour, ones(numel(outlier_tt_neibour),1)*ybou1(2), ybou1(1));
                he_out.FaceColor = col_out;
                he_out.EdgeColor = col_out;
                he_out.FaceAlpha = 0.4;
                he_out.EdgeAlpha = 0.7;
            end
        end

        for j = 1:numel(tt_lea)
            if missing_months(j) <= 1 || missing_months(j) >= numel(MSSA_TS)
                continue
            end

            leak_tt_neibour = [tt_lea(j)-1/24, tt_lea(j), tt_lea(j)+1/24];
            leak_MSL_neibour = [ ...
                (MSSA_TS(missing_months(j)-1) + MSSA_TS(missing_months(j))) / 2, ...
                MSSA_TS(missing_months(j)), ...
                (MSSA_TS(missing_months(j)) + MSSA_TS(missing_months(j)+1)) / 2];

            p2 = plot(leak_tt_neibour, leak_MSL_neibour, 'r', 'linewidth', 2);

            he_gap = area(leak_tt_neibour([1,3]), [ybou1(2); ybou1(2)], ybou1(1));
            he_gap.FaceColor = col_gap;
            he_gap.EdgeColor = col_gap;
            he_gap.FaceAlpha = 0.4;
            he_gap.EdgeAlpha = 0.7;
        end

        p3 = plot(tt_fil, MSSA1_CJGIave, 'blue', 'linewidth', 1);

        local_title(A5_CJGI_Sort, i);
        text(Cap_posix1, Cap_posiy1, ['(' abc(i) ')'], 'FontSize', font_Size-1);

        xlim([tt_fil(1)-0.5, tt_fil(end)+0.5])
        ylim(ybou1)
        xticks(floor(tt_fil(1)):1:floor(tt_fil(end)))
        set(gca, 'xticklabel', get(gca, 'xtick'), 'FontSize', font_Size-2)
        set(gca, 'FontName', 'Times New Roman', 'TickLength', [0.004,0.035]);
    end
end

if fig
    if exist('he_out', 'var') && exist('he_gap', 'var') && exist('p2', 'var')
        l1 = legend([p1,p2,p3,he_gap,he_out], ...
            {'Raw data','Fit data','1st MSSA average','GRACE gaps','Outliers'}, ...
            'NumColumns', 5);
    elseif exist('he_gap', 'var') && exist('p2', 'var')
        l1 = legend([p1,p2,p3,he_gap], ...
            {'Raw data','Fit data','1st M-SSA data','GRACE gaps'}, ...
            'NumColumns', 4);
    else
        l1 = legend([p1,p3], {'Raw data','1st M-SSA data'}, 'NumColumns', 2);
    end

    set(l1, 'Position', [0.25,0.02,0.5,0.03], 'box', 'off', 'fontsize', font_Size)
    set(gcf, 'Position', F_Position)
    print(fl2, '-dtiff', '-r125', [fig_note '_MSSA1st.tif']);
end

% ===== iterative MSSA =====
MSSAn_Total_CGJI_fill_reconst = cell(opts.maxIter + 1, 1);
MSSAn_Total_CGJI_fill_reconst{1} = Total_CGJI_fill_reconst_MSSA1;

MSSAn_std = nan(opts.maxIter + 1, intit_num);
MSSAn_std(1,:) = std(Total_CGJI_fill_reconst_MSSA1, 0, 1);

chi = nan(opts.maxIter + 1, intit_num);
flag = zeros(intit_num, 1);
iter_num = zeros(intit_num, 1);

MSSAn_RC = [];

for ii = 1:opts.maxIter

    MSLAn_CJGI = MSSAn_Total_CGJI_fill_reconst{ii};

    [~,~,~,MSSAn_RC] = MSSA_fast(fill_months, MSLAn_CJGI, M, Method, N, opts);

    MSLAn_CJGI_reconst = squeeze(sum(MSSAn_RC, 2));
    if intit_num == 1
        MSLAn_CJGI_reconst = MSLAn_CJGI_reconst(:);
    end

    for i = find(flag == 0).'

        missing_months = A5_CJGI_Sort(i).missing_months(:);

        if isempty(missing_months)
            flag(i) = 1;
            continue
        end

        MSLAn_c = MSLAn_CJGI(:,i);
        old_gap_values = MSLAn_CJGI(missing_months,i);

        MSLAn_c(missing_months) = MSLAn_CJGI_reconst(missing_months,i);

        MSSAn_std(ii+1,i) = std(MSLAn_c(missing_months));

        denom = MSSAn_std(ii+1,i) * MSSAn_std(ii,i);
        if denom <= eps || isnan(denom)
            chi(ii+1,i) = 0;
        else
            chi(ii+1,i) = sqrt(sum((old_gap_values - MSLAn_c(missing_months)).^2) / denom);
        end

        if chi(ii+1,i) < opts.chiThreshold
            flag(i) = 1;
        end

        MSLAn_CJGI(:,i) = MSLAn_c;
        iter_num(i) = iter_num(i) + 1;
    end

    MSSAn_Total_CGJI_fill_reconst{ii+1} = MSLAn_CJGI;

    if sum(flag) == intit_num
        break
    end
end

MSSAn_Total_CGJI_fill_reconst(max(iter_num)+2:end) = [];

% ===== plot iterative process =====
if fig
    fl3 = figure;
    ha = tight_subplot(intit_num, 1, F_gap, F_marg_h, F_marg_w);

    max_iter_plot = max(max(iter_num), 1);
    map_color = jet(max_iter_plot);

    for i = 1:intit_num

        ybou1 = [-20, 20];

        if isfield(A5_CJGI_Sort, 'use_months')
            use_months = A5_CJGI_Sort(i).use_months(:);
        else
            use_months = fill_months;
        end

        missing_months = A5_CJGI_Sort(i).missing_months(:);
        MSSA_TS = MSLA1_CJGI(:,i);

        tt_use = floor(use_months/12) + data_year_beg + mod(use_months,12)/12 - 1/24;
        tt_lea = floor(missing_months/12) + data_year_beg + mod(missing_months,12)/12 - 1/24;

        axes(ha(i));

        iter_plot_i = max(iter_num(i), 1);
        map_colori = map_color(round(linspace(1, max_iter_plot, iter_plot_i)), :);
        TickLabels_c = cell(iter_plot_i, 1);

        for ii = 1:iter_plot_i
            MSSAi_Total_CGJI_EST_reconst = MSSAn_Total_CGJI_fill_reconst{ii};

            p1 = plot(tt_fil, MSSAi_Total_CGJI_EST_reconst(:,i), ...
                'color', map_colori(ii,:), 'linewidth', 2);
            hold on

            TickLabels_c{ii} = num2str(ii-1);

            ybou1(2) = max(ybou1(2), max(MSSAi_Total_CGJI_EST_reconst(:,i)));
            ybou1(1) = min(ybou1(1), min(MSSAi_Total_CGJI_EST_reconst(:,i)));
        end

        if isempty(tt_use)
            Cap_posix1 = tt_fil(1);
        else
            Cap_posix1 = tt_use(min(5, numel(tt_use)));
        end
        Cap_posiy1 = ybou1(2) - (ybou1(2)-ybou1(1))/10;

        colormap(ha(i), map_colori)
        colorbar('Ticks', 1/iter_plot_i/2:1/iter_plot_i:1, ...
            'TickLabels', TickLabels_c)

        for j = 1:numel(tt_lea)
            leak_tt_neibour = [tt_lea(j)-1/24, tt_lea(j), tt_lea(j)+1/24];

            he_gap = area(leak_tt_neibour([1,3]), [ybou1(2); ybou1(2)], ybou1(1));
            he_gap.FaceColor = col_gap;
            he_gap.EdgeColor = col_gap;
            he_gap.FaceAlpha = 0.4;
            he_gap.EdgeAlpha = 0.7;
            hold on
        end

        text(Cap_posix1, Cap_posiy1, ...
            ['(' abc(i) ') The number of iteration: ' num2str(iter_num(i)-1)], ...
            'FontSize', font_Size-1);

        xlim([tt_fil(1)-0.5, tt_fil(end)+0.5])
        ylim(ybou1)
        xticks(floor(tt_fil(1)):1:floor(tt_fil(end)))
        set(gca, 'xticklabel', get(gca, 'xtick'), 'FontSize', font_Size-2)
        ylabel('cm')

        local_title(A5_CJGI_Sort, i);

        set(gca, 'FontName', 'Times New Roman', 'TickLength', [0.004,0.035]);
    end

    set(gcf, 'Position', F_Position)

    if exist('he_gap', 'var')
        l1 = legend([he_gap], {'GRACE Gap'}, 'NumColumns', 4, 'box', 'off');
        set(l1, 'Position', [0.25,0.02,0.5,0.03], 'box', 'off', 'fontsize', font_Size)
    end

    print(fl3, '-dtiff', '-r125', [fig_note '_MSSAN.tif']);
end

RCtest = [];
if opts.doTest
    RCtest = local_rc_test(MSSAn_RC, N, intit_num);
end

varns = {MSSAn_Total_CGJI_fill_reconst, MSSAn_RC, iter_num, chi, RCtest};
varargout = varns(1:nargout);

end


function local_title(A5_CJGI_Sort, i)

if isfield(A5_CJGI_Sort, 'name')
    name_i = A5_CJGI_Sort(i).name;
else
    name_i = ['Institution ' num2str(i)];
end

if numel(name_i) > 4
    name_i = name_i(1:end-4);
end

if isfield(A5_CJGI_Sort, 'MSSA_TS_order')
    order_i = A5_CJGI_Sort(i).MSSA_TS_order;

    if ischar(order_i)
        title([name_i ' inverted barometer'], 'FontWeight', 'bold')
    else
        title([name_i ' spherical Slepian coefficient (\alpha = ' num2str(order_i) ')'], ...
            'FontWeight', 'bold')
    end
else
    title(name_i, 'FontWeight', 'bold')
end

end
