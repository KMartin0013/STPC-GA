function plot_stpc_spatial_temporal_RMS( ...
    STPC_p_SN_results, Turning_number, Area, XY_ori, XY_buffer, ...
    lon, lat, Noise_SigLev, use_rms, use_rms_title, Attach)

% PLOT_SPATIAL_TEMPORAL_STD
% Plot:
%   1) spatial RMS maps for EWH and EWH(noise-removal) in two panels
%   2) spatio-temporal RMS curves for different S and N
%   3) heatmap of RMS matrix to find optimal (S,N)
%
% Returns:
%   Square_Seq_SIG : the RMS matrix (for disp_sn) used as selection criterion
%

font_Size = 16;
abc       = 'abcdefghijklmnopqrstuvwxyz';

lon_matl=lon-0.5;
lat_matl=lat+0.5;
[lonlon_matl,latlat_matl]=meshgrid(lon_matl,lat_matl);

% Sub-sample S and N
use_ss = 1:2:Turning_number;
use_nn = 1:2:Turning_number;

% Display types to be plotted in curves
% 1: EWH, 3: EWH(noise-removal), 2: EWH(noise)
use_disp = [1,3,2];

%% 1) First spatial panel (EWH, disp_data = 1)

disp_data = 1;

F_Position = [2100,-350,1300,1100];
f = figure;

F_gap0    = [0.05 0.125];
F_marg_h0 = [0.325 0.04];
F_marg_w0 = [0.04 0.36];

% Colorbar bounds depending on Area
if strcmp(Area,'greenland')
    ybou1 = [0 120];
    ybou2 = [0 80];
    ybou3 = [50 65];
elseif strcmp(Area,'SCSpTH')
    ybou1 = [0 25];
    ybou2 = [0 25];
    ybou3 = [0 25];
elseif strcmp(Area,'yangtze')
    ybou1 = [0 20];
    ybou2 = [0 15];
    ybou3 = [0 15];
elseif strcmp(Area,'mekong')
    ybou1 = [0 40];
    ybou2 = [0 40];
    ybou3 = [0 20];
else
    ybou1 = [0 40];
    ybou2 = [0 40]; 
    ybou3 = [0 15];
end

h0 = tight_subplot(numel(use_ss), numel(use_nn), F_gap0, F_marg_h0, F_marg_w0);
set(f, 'Position', F_Position);

ff = 0;

for ss = use_ss
    for nn = use_nn

        ff = ff + 1;
        axes(h0(ff));

        EWH_sn = STPC_p_SN_results{ss,nn}.EWH;

        switch disp_data
            case 1
                r_res   = squeeze(EWH_sn.fillM_Grid_MSSA);
                fig_name = 'fillM_EWH_MSSA';
            case 2
                r_res   = squeeze(EWH_sn.fillM_Grid_MSSA - EWH_sn.fillM_Grid_STPC);
                fig_name = 'fillM_EWH_MSSA_noise';
            case 3
                r_res   = squeeze(EWH_sn.fillM_Grid_STPC);
                fig_name = 'fillM_EWH_STPC';
            case 4
                r_res   = squeeze(EWH_sn.fillM_Grid_sig_STPC);
                fig_name = 'fillM_EWH_signal_STPC';
            case 5
                r_res   = squeeze(EWH_sn.fillM_Grid_res_STPC);
                fig_name = 'fillM_EWH_resid_STPC';
        end

        s_rms = squeeze(rms(r_res, 1));

        % Map projection
        if strcmp(Area,'greenland')
            m_proj('stereographic','long',(lon_matl(1)+lon_matl(end))/2, ...
                'lat',(lat_matl(1)+lat_matl(end))/2,'rad',20,'rec','on');
        else
            m_proj('mercator','long',[lon_matl(1), lon_matl(end)], ...
                'lat',[lat_matl(end), lat_matl(1)]);
        end

        m_pcolor(lonlon_matl, latlat_matl, s_rms);
        shading flat
        m_coast

        m_line(XY_ori(:,1), XY_ori(:,2),'Linestyle','-','color','k','linewidth',1);
        hold on
        m_line(XY_buffer(:,1), XY_buffer(:,2),'Linestyle','--','color','r','linewidth',1);

        if nn==use_nn(1) && ss==use_ss(end)
            m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','xtick',lon_matl(6:10:end));
        elseif nn==use_nn(1)
            m_grid('tickdir','in','yaxislocation','left','xticklabels',[]);
        elseif ss==use_ss(end)
            m_grid('tickdir','in','xaxislocation','bottom','yticklabels',[],'xtick',lon_matl(6:10:end));
        else
            m_grid('tickdir','in','xticklabels',[],'yticklabels',[])
        end

%         m_text(lon_matl(end)+55, lat_matl(1)-9, ['(' abc(ff) ')'], ...
%             'color','black','fontsize',font_Size,'FontName','Times New Roman');

%         m_text(lon_matl(end), lat_matl(1)+5, ['(' abc(ff) ')'], ...
%             'color','black','fontsize',font_Size,'FontName','Times New Roman');

        m_text(mean(lon_matl), lat_matl(1)+5, ['(' abc(ff) ') S(' num2str(ss) ') and N(' num2str(nn) ')'], ...
            'color','black','fontsize',font_Size,'FontName','Times New Roman');

        caxis(ybou1);
        set(gca,'FontName','Times New Roman','FontSize',font_Size-2, ...
            'TickLength',[0.012,0.055]);

        if disp_data == 1
            title('EWH','FontName','Times New Roman','FontSize',font_Size-4);
        elseif disp_data == 3
            title(['EWH (p=' num2str(Noise_SigLev) ')'],'FontName','Times New Roman','FontSize',font_Size-4);
        end
    end
    colormap('jet');
end

%% 2) Second spatial panel (EWH noise-removal, disp_data = 3)

disp_data = 3;

F_gap0e    = F_gap0;
F_marg_h0e = F_marg_h0;
F_marg_w0e = [F_marg_w0(1)+0.115, F_marg_w0(2)-0.115];

h0e = tight_subplot(numel(use_ss), numel(use_nn), F_gap0e, F_marg_h0e, F_marg_w0e);

ff = 0;

for ss = use_ss
    for nn = use_nn

        ff = ff + 1;
        axes(h0e(ff));

        EWH_sn = STPC_p_SN_results{ss,nn}.EWH;
        switch disp_data
            case 1
                r_res   = squeeze(EWH_sn.fillM_Grid_MSSA);
                fig_name = 'fillM_EWH_MSSA';
            case 2
                r_res   = squeeze(EWH_sn.fillM_Grid_MSSA - EWH_sn.fillM_Grid_STPC);
                fig_name = 'fillM_EWH_MSSA_noise';
            case 3
                r_res   = squeeze(EWH_sn.fillM_Grid_STPC);
                fig_name = 'fillM_EWH_STPC';
            case 4
                r_res   = squeeze(EWH_sn.fillM_Grid_sig_STPC);
                fig_name = 'fillM_EWH_signal_STPC';
            case 5
                r_res   = squeeze(EWH_sn.fillM_Grid_res_STPC);
                fig_name = 'fillM_EWH_resid_STPC';
        end

        s_rms = squeeze(rms(r_res, 1));

        if strcmp(Area,'greenland')
            m_proj('stereographic','long',(lon_matl(1)+lon_matl(end))/2, ...
                'lat',(lat_matl(1)+lat_matl(end))/2,'rad',20,'rec','on');
        else
            m_proj('mercator','long',[lon_matl(1), lon_matl(end)], ...
                'lat',[lat_matl(end), lat_matl(1)]);
        end

        m_pcolor(lonlon_matl, latlat_matl, s_rms);
        shading flat
        m_coast

        m_line(XY_ori(:,1), XY_ori(:,2),'Linestyle','-','color','k','linewidth',1);
        hold on
        m_line(XY_buffer(:,1), XY_buffer(:,2),'Linestyle','--','color','r','linewidth',1);

        if ss==use_ss(end)
            m_grid('tickdir','in','yticklabels',[],'yaxislocation','left','xaxislocation','bottom','xtick',lon_matl(6:10:end));
        else
            m_grid('tickdir','in','xticklabels',[],'yticklabels',[])
        end

        caxis(ybou1);
        set(gca,'FontName','Times New Roman','FontSize',font_Size-2, ...
            'TickLength',[0.012,0.055]);

        if disp_data == 1
            title('EWH','FontName','Times New Roman','FontSize',font_Size-4);
        elseif disp_data == 3
            title(['EWH (p=' num2str(Noise_SigLev) ')'],'FontName','Times New Roman','FontSize',font_Size-4);
        end
    end
    colormap('jet');
end

cb = colorbar;
cb.Location = 'South';
cb.Position = [0.10 0.275 0.6 0.015];
title(cb,'RMS','fontsize',14);

%% 3) Spatio-temporal RMS curves + heatmap

disp_sn = 3;  % which disp_data is interpreted as "signal"

% Square matrix for selected S,N subset
use_ss_nn_square = squeeze(use_rms(disp_sn, :, :));

% Colors (color-blind friendly)
colorr_blind = [0,114,178; ...
                204,121,167; ...
                0,158,115; ...
                230,159,0; ...
                213,94,0; ...
                240,228,66; ...
                86,180,233]/255;

% Tick labels
xticklabelS     = cell(1,length(use_ss));
xticklabelN     = cell(1,length(use_nn));
xticklabelSN    = cell(1,length(use_nn));
xticklabelS_gap = cell(1,max(use_ss));
xticklabelN_gap = cell(1,max(use_nn));

for i = 1:length(use_ss)
    xticklabelS{i}          = ['S(' num2str(use_ss(i)) ')'];
    xticklabelS_gap{use_ss(i)} = xticklabelS{i};
end
for i = 1:length(use_nn)
    xticklabelN{i}          = ['N(' num2str(use_nn(i)) ')'];
    xticklabelN_gap{use_nn(i)} = xticklabelN{i};
    xticklabelSN{i}         = ['S(' num2str(use_ss(i)) '), N(' num2str(use_nn(i)) ')'];
end

% 3.1: RMS vs S for fixed N
F_gap1    = [0.05 0.06];
F_marg_h1 = [0.06 0.77];
F_marg_w1 = [F_marg_w0(1)+0.02, F_marg_w0e(2)+0.02];

h1 = tight_subplot(1, length(use_nn), F_gap1, F_marg_h1, F_marg_w1);

xbou1 = [0.5, Turning_number+0.5];

for i = 1:length(use_nn)

    ff = ff + 1;
    axes(h1(i));

    p1 = gobjects(1, max(use_disp));
    for j = use_disp
        p1(j) = plot(1:Turning_number, squeeze(use_rms(j,:,use_nn(i))), ...
            'Marker','none','LineStyle','-', ...
            'LineWidth',2,'color',colorr_blind(j,:));
        hold on
    end

    if i == 1
        ylabel(use_rms_title);
    end

    tbox = text(1, ybou2(2)*4.5/5, ...
        ['N=N(' num2str(use_nn(i)) ')'], ...
        'FontName','Times New Roman','FontSize',font_Size-3);
    tbox.EdgeColor = 'k';
    tbox.LineStyle = '-';
    tbox.LineWidth = 1;

    ylim(ybou2);
    xlim(xbou1);
    xticks(use_ss);
    xticklabels(xticklabelS);

    title(['(' abc(ff) ')'], 'color','black','fontsize',font_Size, ...
        'FontName','Times New Roman');

    set(gca,'FontName','Times New Roman','FontSize',font_Size-2, ...
        'TickLength',[0.012,0.055]);
end

% 3.2: RMS vs N for fixed S
F_gap2    = [0.06 0.03];
F_marg_h2 = [0.34 0.03];
F_marg_w2 = [0.81 0.02];

h2 = tight_subplot(length(use_ss), 1, F_gap2, F_marg_h2, F_marg_w2);

for i = 1:length(use_ss)

    ff = ff + 1;
    axes(h2(i));

    p1 = gobjects(1, max(use_disp));
    for j = use_disp
        p1(j) = plot(1:Turning_number, squeeze(use_rms(j, use_ss(i), :)), ...
            'Marker','none','LineStyle','-', ...
            'LineWidth',2,'color',colorr_blind(j,:));
        hold on
    end

    ylabel(use_rms_title);

    tbox = text(1, ybou2(2)*4.5/5, ...
        ['S=S(' num2str(use_ss(i)) ')'], ...
        'FontName','Times New Roman','FontSize',font_Size-3);
    tbox.EdgeColor = 'k';
    tbox.LineStyle = '-';
    tbox.LineWidth = 1;

    ylim(ybou2);
    xlim(xbou1);
    xticks(use_nn);
    xticklabels(xticklabelN);

    title(['(' abc(ff) ')'], 'color','black','fontsize',font_Size, ...
        'FontName','Times New Roman');

    set(gca,'FontName','Times New Roman','FontSize',font_Size-2, ...
        'TickLength',[0.012,0.055]);
end

% 3.3: Heatmap of RMS matrix
F_gap3    = F_gap0;
F_marg_h3 = F_marg_h1;
F_marg_w3 = [F_marg_w2(1)-0.01, F_marg_w2(2)];

h3 = tight_subplot(1, 1, F_gap3, F_marg_h3, F_marg_w3);
axes(h3(1));

ff = ff + 1;

SHM = SHeatmap(use_ss_nn_square, 'Format', 'sq');
SHM = SHM.draw();
SHM.setText('FontSize', font_Size+2);

% Mark global maximum with "x"
max_val = max(use_ss_nn_square, [], 'all');
for i = 1:size(use_ss_nn_square,1)
    for j = 1:size(use_ss_nn_square,2)
        if use_ss_nn_square(i,j) == max_val
            SHM.setTextMN(i,j,'String','x');
        else
            SHM.setTextMN(i,j,'String','');
        end
    end
end

title(['(' abc(ff) ')'], 'color','black','fontsize',font_Size, ...
    'FontName','Times New Roman');

ax = gca;
ax.XTickLabel = xticklabelN_gap;
ax.YTickLabel = xticklabelS_gap;
ax.FontSize   = 14;

cb = colorbar;
cb.Location = 'eastoutside';
cb.Position = [0.965 0.06 0.012 0.17];
title(cb, {use_rms_title}, 'fontsize',14,'Interpreter','latex');

caxis(ybou3);

set(gca,'FontName','Times New Roman','FontSize',font_Size-2, ...
    'TickLength',[0.012,0.055]);

% Legend
legend_text = {'EWH',['noise EWH (p=' num2str(Noise_SigLev) ')'], ...
    ['EWH (p=' num2str(Noise_SigLev) ')'], ...
    ['Fitted signal of EWH (p=' num2str(Noise_SigLev) ')'], ...
    ['Unfitted signal of EWH (p=' num2str(Noise_SigLev) ')']};

% Use only the last set of handles 'p1' from previous loops (non-empty indices use_disp)
pn = p1(use_disp);

legend(pn, legend_text(use_disp), ...
    'box','off','Location',[0.15,0.013,0.7,0.02], ...
    'FontSize',font_Size,'NumColumns',3);

% Save figure
tif_name = [fig_name, '_Spatio_temporal_SIG' num2str(Noise_SigLev) '.tif'];
print(f, '-dtiff', '-r500', fullfile(Attach.fig_path_ALL, tif_name));

end
