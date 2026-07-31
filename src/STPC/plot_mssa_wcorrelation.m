function plot_mssa_wcorrelation(wcorrDiagnostics, turn_MSSA_four, Attach)
%PLOT_MSSA_WCORRELATION Plot compact w-correlation and score diagnostics.

coefficients = wcorrDiagnostics.coefficients;
coefficientCount = numel(coefficients);
showCount = min(coefficientCount, 12);
if showCount < 1
    return
end

fontSize = 9;
matrixFigure = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100, 100, 1200, 850]);
columns = min(4, showCount);
rows = ceil(showCount / columns);
tiledlayout(rows, columns, 'TileSpacing', 'compact', ...
    'Padding', 'compact');

for ss = 1:showCount
    nexttile;
    matrix = coefficients{ss}.wcorr;
    imagesc(abs(matrix), [0, 1]);
    axis image
    set(gca, 'YDir', 'normal', 'FontSize', fontSize);
    hold on
    candidates = unique(turn_MSSA_four{ss}(:, 1));
    for i = 1:numel(candidates)
        xline(candidates(i) + 0.5, 'b:', 'LineWidth', 0.8);
        yline(candidates(i) + 0.5, 'b:', 'LineWidth', 0.8);
    end
    title(sprintf('SSF %d; N=%s', ss, ...
        strjoin(string(candidates'), ',')), 'FontWeight', 'normal');
    xlabel('RC index');
    ylabel('RC index');
end
colormap(matrixFigure, flipud(gray(256)));
matrixName = [Attach.Attach_ALL 'WCorrelation_MSSA.tif'];
print(matrixFigure, '-dtiff', '-r200', ...
    fullfile(Attach.fig_path_ALL, matrixName));
close(matrixFigure);

scoreFigure = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100, 100, 1050, 650]);
hold on
colors = parula(max(coefficientCount, 2));
for ss = 1:coefficientCount
    diagnostic = coefficients{ss};
    plot(diagnostic.domain, diagnostic.score, ...
        'Color', colors(ss, :), 'LineWidth', 0.8);
    candidates = unique(turn_MSSA_four{ss}(:, 1));
    [present, locations] = ismember(candidates, diagnostic.domain);
    plot(candidates(present), diagnostic.score(locations(present)), ...
        'o', 'Color', colors(ss, :), 'MarkerSize', 3, ...
        'MarkerFaceColor', colors(ss, :));
end
yline(0, 'k:');
xlabel('Retained MSSA reconstructed components, N');
ylabel('Tail RMS(|\rho_w|) - cross-boundary RMS(|\rho_w|)');
title(sprintf(['W-correlation boundary scores (%d SSFs; ' ...
    'block width %d)'], coefficientCount, ...
    wcorrDiagnostics.requested_block_width));
set(gca, 'FontSize', 11, 'Box', 'on');
scoreName = [Attach.Attach_ALL 'WCorrelation_Score_MSSA.tif'];
print(scoreFigure, '-dtiff', '-r250', ...
    fullfile(Attach.fig_path_ALL, scoreName));
close(scoreFigure);
end
