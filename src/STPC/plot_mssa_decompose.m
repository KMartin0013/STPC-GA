function plot_mssa_decompose(mssa_Sort_sn, mssa_Sort_coff_sn, ...
            r_record, lonlon, c11cmn, Area, XY_ori, XY_buffer, ...
            V, MSSA_evalues, S_rec, S_ol, intit_num, Attach)

fill_months     = mssa_Sort_sn(1).fill_months;
use_months      = mssa_Sort_sn(1).use_months;
data_year_beg   = mssa_Sort_sn(1).data_year_beg;
missing_months  = mssa_Sort_sn(1).missing_months;

lon_matl = c11cmn(1)-0.5 : c11cmn(3)-0.5;
lat_matl = c11cmn(2)+0.5 : -1 : c11cmn(4)+0.5;
[lonlon_matl,latlat_matl] = meshgrid(lon_matl,lat_matl);

select_S = 1:S_rec;      % how many coefficients (plots)
Nc       = 6;            % the number of RCs for the plot

col_gap  = [127 127 127]/255;

font_Size = 16;
MarkSize  = 2;
F_Position = [2100,-250,1500,1100];
abc = 'abcdefghijklmnopqrstuvwxyz';

% Blind-friendly color
colorr_blind = [0,114,178;204,121,167;0,158,115;230,159,0;213,94,0;240,228,66;86,180,233]/255;

% time epoch
tt_use = floor(use_months/12)+data_year_beg+mod(use_months,12)/12-1/24;
tt_fil = floor(fill_months/12)+data_year_beg+mod(fill_months,12)/12-1/24;
tt_lea = floor(missing_months/12)+data_year_beg+mod(missing_months,12)/12-1/24;

%--------------------- Main loops ---------------------%
for i = 1:numel(select_S)

    Nc_plot = min(Nc, size(mssa_Sort_coff_sn(intit_num+1).fillM_reconcoffs_RC{i},2));

    f2 = figure;
    set(f2,'Position',F_Position)

    SSF = select_S(i);

    %% (a) Normalized slepian basis
    F_gap    = [.02 .03];
    F_marg_h = [.29 .24];
    F_marg_w = [.03 .78];
    h0 = tight_subplot(1,1,F_gap,F_marg_h,F_marg_w);

    SSF_flat = reshape(squeeze(r_record(SSF,:,:)),[],1);
    SSF_flat = normalize(SSF_flat,'norm',Inf);
    SSF_nor  = reshape(SSF_flat,size(lonlon,1),size(lonlon,2));

    axes(h0(1));

    if strcmp(Area,'greenland')
        m_proj('stereographic','long',(c11cmn(1)+c11cmn(3))/2,'lat',(c11cmn(4)+c11cmn(2))/2, ...
               'rad',20,'rec','on');
        tt_lon = (c11cmn(3)-c11cmn(1))/5 + c11cmn(1);
        tt_lat = c11cmn(4) + (c11cmn(2)-c11cmn(4))/12;
    elseif strcmp(Area,'yangtze')
        tt_lon = (c11cmn(3)-c11cmn(1))/20 + c11cmn(1);
        tt_lat = c11cmn(4) + (c11cmn(2)-c11cmn(4))/6;
        m_proj('mercator','long',[c11cmn(1)-0.5, c11cmn(3)-0.5],'lat',[c11cmn(4)+0.5, c11cmn(2)+0.5]);
    else
        tt_lon = (c11cmn(3)-c11cmn(1))/12 + c11cmn(1);
        tt_lat = c11cmn(2) - (c11cmn(2)-c11cmn(4))/12;
        m_proj('mercator','long',[c11cmn(1)-0.5, c11cmn(3)-0.5],'lat',[c11cmn(4)+0.5, c11cmn(2)+0.5]);
    end

    m_pcolor(lonlon_matl,latlat_matl,SSF_nor);
    shading flat;
    colormap(mymap("coolwarm"));
    h = colorbar('h');

    m_coast
    m_grid('box','fancy','tickdir','in');
    m_line(XY_ori(:,1),   XY_ori(:,2),   'LineStyle','-','color','k','linewidth',1);
    hold on
    m_line(XY_buffer(:,1),XY_buffer(:,2),'LineStyle','--','color','r','linewidth',1);

    m_text(tt_lon,tt_lat, ...
        ['\lambda = ' num2str(V(i),'%.2f')],'FontName','Times New Roman', ...
        'fontsize',13,'color','b','verticalalignment','middle','horizontalalignment','left', ...
        'backgroundColor','white')

    caxis([-1 1]);

    if strcmp(Area,'greenland')
        m_text(tt_lon,90, ...
            ['\bf (a) Normalized Slepian basis ($g_{' num2str(SSF) '}\left(\hat{\mathbf{r}}\right)$)'], ...
            'FontWeight','bold','FontSize',font_Size-1,'horizontalalignment','center','Interpreter','latex')
    else
        title(['\bf (a) Normalized Slepian basis ($g_{' num2str(SSF) '}\left(\hat{\mathbf{r}}\right)$)'], ...
            'FontWeight','bold','FontSize',font_Size-1,'Interpreter','latex')
    end
    set(gca,'FontName','Times New Roman','FontSize',font_Size-2,'TickLength',[0.004,0.035]);

    %% (b) original time series of slepian coefficient
    F_gap    = [.02 .03];
    F_marg_h = [.82 .04];
    F_marg_w = [.3 .3];
    ha = tight_subplot(1,1,F_gap,F_marg_h,F_marg_w);

    ybou1 = [-1,1];
    ybou2 = [-20,20];

    axes(ha(1));

    pn = [];
    pn_name = {};
    line_ins = [];

    for ii = 1:intit_num
        line_ins(ii,:) = squeeze(mssa_Sort_coff_sn(ii).fillM_deltacoffs(:,i));
        p1(ii) = plot(tt_fil,line_ins(ii,:), ...
            'color',colorr_blind(ii,:),'linewidth',2);
        hold on

        pn      = [pn p1(ii)];
        pn_name{ii} = mssa_Sort_sn(ii).name(1:end-4);
    end

    Maxx = max(line_ins,[],'all');
    Minn = min(line_ins,[],'all');
    ybou1(1) = min(ybou1(1),Minn);
    ybou1(2) = max(ybou1(2),Maxx);

    Cap_posix1 = tt_fil(5);
    Cap_posiy1 = ybou1(2)-(ybou1(2)-ybou1(1))/10;

    for k = 1:numel(tt_lea)
        leak_tt_neibour = [tt_lea(k)-1/24,tt_lea(k),tt_lea(k)+1/24];
        he_gap = area(leak_tt_neibour([1,3]),[ybou1(2);ybou1(2)],ybou1(1));
        he_gap.FaceColor = col_gap;
        he_gap.EdgeColor = col_gap;
        he_gap.FaceAlpha = 0.4;
        he_gap.EdgeAlpha = 0.7;
        hold on
    end

    text(Cap_posix1,Cap_posiy1,['(' abc(2) ')'],'FontSize',font_Size-1);

    ylabel('EWH (cm)','FontWeight','bold','FontSize',font_Size)

    xlim([tt_fil(1)-0.5, tt_fil(end)+0.5])
    ylim(ybou1)

    pn_name{intit_num+1} = 'GRACE gap';
    xticks(floor(tt_fil(1)):2:floor(tt_fil(end)))
    xticklabels([]);

    title(['\bf Spherical Slepian coefficient ($d_{' num2str(SSF) '}$)'],'Interpreter','latex')

    set(gca,'FontName','Times New Roman','FontSize',font_Size-2,'TickLength',[0.004,0.035]);

    %% (c) Power Spectrum of the Coefficient
    F_gap    = [.02 .03];
    F_marg_h = [.82 .04];
    F_marg_w = [.74 .02];
    ha = tight_subplot(1,1,F_gap,F_marg_h,F_marg_w);

    Fs   = 12;
    T    = 1/Fs;
    L    = length(tt_fil);
    freq = Fs*(0:(L/2))/L;

    ybou1 = [0,0.5];

    axes(ha(1));

    Yss = zeros(size(line_ins(1,:)));
    for ii = 1:intit_num
        Yss = Yss + squeeze(line_ins(ii,:));

        Y  = fft(squeeze(line_ins(ii,:)));
        P2 = abs(Y/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);

        semilogx(1./freq,P1,'color',colorr_blind(ii,:),'LineWidth',1.2, ...
                 'Marker','o','MarkerSize',MarkSize);
        hold on

        ybou1(2) = max(ybou1(2),max(P1));
    end

    xticks([3/12,6/12,1,4,10]);
    xticklabels({'3','6','12','48','120'})
    xticklabels([]);

    ylim(ybou1)

    Cap_posiy1 = ybou1(2)-(ybou1(2)-ybou1(1))/10;

    text(0.12,Cap_posiy1,['(' abc(2+Nc_plot+1) ')'],'FontSize',font_Size-1);

    title('Power spectrum','FontWeight','bold','FontSize',font_Size-1)
    set(gca,'FontName','Times New Roman','FontSize',font_Size-2,'TickLength',[0.007,0.035]);

    %% (d) Nc_plot RC time series
    F_gap    = [.02 .03];
    F_marg_h = [.06 .22];
    F_marg_w = [.3 .3];
    ha = tight_subplot(Nc_plot,1,F_gap,F_marg_h,F_marg_w);

    for j = 1:Nc_plot

        ybou1 = [-1,1];

        axes(ha(j));

        pn = [];
        pn_name = {};
        line_ins = [];
        for ii = 1:intit_num
            line_ins(ii,:) = squeeze(mssa_Sort_coff_sn(ii).fillM_reconcoffs_RC{i}(:,j));
            p1(ii) = plot(tt_fil,squeeze(mssa_Sort_coff_sn(ii).fillM_reconcoffs_RC{i}(:,j)), ...
                          'color',colorr_blind(ii,:),'linewidth',2);
            hold on
            pn = [pn p1(ii)];
            pn_name{ii} = mssa_Sort_sn(ii).name(1:end-4);
        end

        Maxx = max(line_ins,[],'all');
        Minn = min(line_ins,[],'all');
        ybou1(1) = min(ybou1(1),Minn);
        ybou1(2) = max(ybou1(2),Maxx);

        Cap_posix1 = tt_fil(5);
        Cap_posiy1 = ybou1(2)-(ybou1(2)-ybou1(1))/7;

        for k = 1:numel(tt_lea)
            leak_tt_neibour = [tt_lea(k)-1/24,tt_lea(k),tt_lea(k)+1/24];
            he_gap = area(leak_tt_neibour([1,3]),[ybou1(2);ybou1(2)],ybou1(1));
            he_gap.FaceColor = col_gap;
            he_gap.EdgeColor = col_gap;
            he_gap.FaceAlpha = 0.4;
            he_gap.EdgeAlpha = 0.7;
            hold on
        end

        text(Cap_posix1,Cap_posiy1, ...
            ['\bf (' abc(2+j) ') $\widetilde{\lambda}=' num2str(MSSA_evalues(j,i),'%.3f') '$'], ...
            'FontName','Times New Roman','FontSize',font_Size-1,'Interpreter','latex');

        ylabel(['RC-' num2str(j)],'FontWeight','bold','FontSize',font_Size)

        xlim([tt_fil(1)-0.5, tt_fil(end)+0.5])
        ylim(ybou1)

        pn_name{intit_num+1} = 'GRACE gap';
        xticks(floor(tt_fil(1)):2:floor(tt_fil(end)))

        if j == 1
            if (SSF == S_ol) && (S_ol > S_rec)
                title('\bf Leading reconstructed components (RCs) for inverted barometer','FontWeight','bold')
            else
                title(['\bf Leading reconstructed components (RCs) for $d_{' num2str(SSF) '}$'], ...
                      'Interpreter','latex')
            end
        end

        if j == Nc_plot
            xlabel('Time (year)')
        else
            xticklabels([]);
        end

        set(gca,'FontName','Times New Roman','FontSize',font_Size-2,'TickLength',[0.004,0.035]);
    end

    %% (e) Power spectrum for each RC + normality test information
    F_gap    = [.02 .03];
    F_marg_h = [.06 .22];
    F_marg_w = [.74 .02];
    ha = tight_subplot(Nc_plot,1,F_gap,F_marg_h,F_marg_w);

    Fs   = 12;
    T    = 1/Fs;
    L    = length(tt_fil);
    freq = Fs*(0:(L/2))/L;

    RCTest = mssa_Sort_coff_sn(intit_num+1).fillM_reconcoffs_RCTest{i};

    for j = 1:Nc_plot

        axes(ha(j));

        Yss = zeros(size(mssa_Sort_coff_sn(1).fillM_reconcoffs_RC{i}(:,j)))';

        ybou1 = [0,0.5];

        for ii = 1:intit_num
            line_ins(ii,:) = squeeze(mssa_Sort_coff_sn(ii).fillM_reconcoffs_RC{i}(:,j));

            Yss = Yss + squeeze(line_ins(ii,:));

            Y  = fft(squeeze(line_ins(ii,:)));
            P2 = abs(Y/L);
            P1 = P2(1:L/2+1);
            P1(2:end-1) = 2*P1(2:end-1);

            semilogx(1./freq,P1,'color',colorr_blind(ii,:),'LineWidth',1.2, ...
                     'Marker','o','MarkerSize',MarkSize);
            hold on

            ybou1(2) = max(ybou1(2),max(P1));
        end

        xticks([3/12,6/12,1,4,10]);
        xticklabels({'3','6','12','48','120'})
        if j < Nc_plot
            xticklabels([]);
        else
            xlabel('Period (month)')
        end

        ylim(ybou1);

        Cap_posiy1 = ybou1(2)-(ybou1(2)-ybou1(1))/10;
        Cap_posiy2 = ybou1(2)-(ybou1(2)-ybou1(1))/3;

        text(0.12,Cap_posiy1,['(' abc(2+j+Nc_plot+1) ')'],'FontSize',font_Size-1);

        hcdf_fre = RCTest(4,:);
        pks      = RCTest(8,:);
        pli      = RCTest(9,:);

        test_x = 2;
        if hcdf_fre(j) < 1
            test_x = 0.3;
        end

        if     pks(j) < 0.05
            text(test_x,Cap_posiy1,'p(K-S)<0.05','FontSize',font_Size-3);
        elseif pks(j) < 0.1
            text(test_x,Cap_posiy1,'p(K-S)<0.1','FontSize',font_Size-3);
        elseif pks(j) < 0.3
            text(test_x,Cap_posiy1,'p(K-S)<0.3','FontSize',font_Size-3,'color',[0.4940 0.1840 0.5560]);
        else
            text(test_x,Cap_posiy1,'p(K-S)>0.3','FontSize',font_Size-3,'color','r');
        end

        if     pli(j) < 0.05
            text(test_x,Cap_posiy2,'p(Lilliefors)<0.05','FontSize',font_Size-3);
        elseif pli(j) < 0.1
            text(test_x,Cap_posiy2,'p(Lilliefors)<0.1','FontSize',font_Size-3);
        elseif pli(j) < 0.3
            text(test_x,Cap_posiy2,'p(Lilliefors)<0.3','FontSize',font_Size-3,'color',[0.4940 0.1840 0.5560]);
        else
            text(test_x,Cap_posiy2,'p(Lilliefors)>0.3','FontSize',font_Size-3,'color','r');
        end

        set(gca,'FontName','Times New Roman','FontSize',font_Size-2,'TickLength',[0.007,0.035]);
    end

    %% Legend and Output
    if any(missing_months)
        pn_name{intit_num+1} = 'GRACE gap';
        l1 = legend([pn,he_gap],pn_name,'NumColumns',2,'box','off','Location','northwest');
    else
        l1 = legend(pn,pn_name,'NumColumns',2,'box','off','Location','northwest');
    end
    set(l1,'Position',[0.01,0.2,0.2,0.01],'box','off','fontsize',font_Size+1)

    tif_name2 = ['MSSA_Decompose_CC' num2str(SSF) '.tif'];
    if i == 4
        print(f2,'-dtiff','-r500',fullfile(Attach.fig_path_ALL,tif_name2));
    else
        print(f2,'-dtiff','-r300',fullfile(Attach.fig_path_ALL,tif_name2));
    end

end

end
