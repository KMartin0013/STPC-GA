function M_rec = det_opt_M(mssaSort, mssaInfo, slepcoff_gapfilling, ...
    S_ol, S, Attach, doPlot)

intit_num   = mssaInfo.intit_num;

fill_dates  = mssaSort(1).fill_dates;
fill_nomths = numel(fill_dates);
ddr1        = Attach.fig_path_ALL;
% doPlot      = true;

%% determine the suitable M by the separation number
% three years and up to the maximum
test_M_rec  = [36:12:numel(fill_dates)/2];
test_M_num  = numel(test_M_rec);

if ~any(test_M_rec)
    warning('Notice: This code is for monthly datasets now.')
    error('Your data lengths had better exceeds 72 months.')
end

RC_Seq=[];RC_RMS_Seq=[];
for i=1:test_M_num

    disp(['M = ' num2str(test_M_rec(i))]);

    for ss=1:S_ol

        mssaSort_optM=mssaSort(1:intit_num);

        for ins=1:intit_num
            if ss<S+1
                mssaSort_optM(ins).MSSA_TS_order=ss;

            else
                disp('we need to consider IB for ocean study.')
                mssaSort_optM(ins).MSSA_TS_order='IB';

            end
        end

        if intit_num*test_M_rec(i)<50
            warning('You used too small time series.')
            test_nssa=intit_num*test_M_rec(i);
        else
            test_nssa=50;
        end

        [~, MSSAn_RC, MSSAf_evalues, RCTest]=MSSA_final_noCDF_freqsort(mssaSort_optM, ...
            slepcoff_gapfilling(:,:,ss),test_M_rec(i),test_nssa);

        RC_Seq(i,ss,:)=mean(RCTest(5:7,:)');
        RC_RMS_Seq(i,ss,:)=mean(RCTest(5:7,:)'.*rms(mean(MSSAn_RC,3))');
        RC_eva_Seq(i,ss,:)=mean(RCTest(5:7,:)'.*MSSAf_evalues(1:test_nssa));
    end

end


%%
Seq12_13_23=zeros(test_M_num,3);
WSeq12_13_23=zeros(test_M_num,3);
ESeq12_13_23=zeros(test_M_num,3);

for i=1:test_M_num
    lg_str{i}=num2str(test_M_rec(i));
end

for i=1:test_M_num

    Seq12_13_23(i,1)=mean(squeeze(RC_Seq(i,1:S,1)));
    Seq12_13_23(i,2)=mean(squeeze(RC_Seq(i,1:S,2)));
    Seq12_13_23(i,3)=mean(squeeze(RC_Seq(i,1:S,3)));

    WSeq12_13_23(i,1)=mean(squeeze(RC_RMS_Seq(i,1:S,1)));
    WSeq12_13_23(i,2)=mean(squeeze(RC_RMS_Seq(i,1:S,2)));
    WSeq12_13_23(i,3)=mean(squeeze(RC_RMS_Seq(i,1:S,3)));

    ESeq12_13_23(i,1)=mean(squeeze(RC_eva_Seq(i,1:S,1)));
    ESeq12_13_23(i,2)=mean(squeeze(RC_eva_Seq(i,1:S,2)));
    ESeq12_13_23(i,3)=mean(squeeze(RC_eva_Seq(i,1:S,3)));
end

% simple version, use the maximum to determine the M
[WSeq_max,WSeq_num]=max(WSeq12_13_23);
test_M_rec_use=WSeq_num(1);
M_rec=test_M_rec(WSeq_num(1));

[ESeq_max,ESeq_num]=max(ESeq12_13_23);
test_M_rec_use_E=ESeq_num(1);
M_rec_E=test_M_rec(ESeq_num(1));

if doPlot

    figure
    subplot(3,3,1)
    for i=1:test_M_num
        plot(1:S,squeeze(RC_Seq(i,1:S,1)));
        hold on

    end
    title('First Separation number (1/2)')
    legend(lg_str)

    subplot(3,3,2)
    for i=1:test_M_num
        plot(1:S,squeeze(RC_Seq(i,1:S,2)));
        hold on

    end
    title('Second Separation number (1/3)')
    legend(lg_str)

    subplot(3,3,3)
    for i=1:test_M_num
        plot(1:S,squeeze(RC_Seq(i,1:S,3)));
        hold on

    end
    title('Third Separation number (2/3)')
    legend(lg_str)

    subplot(3,3,4)
    for i=1:test_M_num

        plot(1:S,squeeze(RC_RMS_Seq(i,1:S,1)));
        hold on

    end
    title('First Separation number (1/2)')
    legend(lg_str)

    subplot(3,3,5)
    for i=1:test_M_num

        plot(1:S,squeeze(RC_RMS_Seq(i,1:S,2)));
        hold on

    end
    title('Second Separation number (1/3)')
    legend(lg_str)

    subplot(3,3,6)
    for i=1:test_M_num

        plot(1:S,squeeze(RC_RMS_Seq(i,1:S,3)));
        hold on

    end
    title('Third Separation number (2/3)')
    legend(lg_str)

        subplot(3,3,7)
    for i=1:test_M_num

        plot(1:S,squeeze(RC_eva_Seq(i,1:S,1)));
        hold on

    end
    title('First Separation number (1/2)')
    legend(lg_str)

    subplot(3,3,8)
    for i=1:test_M_num

        plot(1:S,squeeze(RC_eva_Seq(i,1:S,2)));
        hold on

    end
    title('Second Separation number (1/3)')
    legend(lg_str)

    subplot(3,3,9)
    for i=1:test_M_num

        plot(1:S,squeeze(RC_eva_Seq(i,1:S,3)));
        hold on

    end
    title('Third Separation number (2/3)')
    legend(lg_str)

end

%% plot for the optimal M

if doPlot

    font_Size=14;F_Position=[100,250,1000,450];
    F_gap=[.03 .08];F_marg_h=[.13 .05];F_marg_w=[.07 .02];

    fm=figure

    ha=tight_subplot(1,2,F_gap,F_marg_h,F_marg_w);

    set(gcf,'Position',F_Position)

    axes(ha(1))
    for i=1:test_M_num
        %     plot(1:N_ol,squeeze(RC_Seq(i,:,1))./N_rec);
        %     hold on
        %     WSeq12_13_23(i,1)=mean(squeeze(RC_Seq(i,:,1))./N_rec);

        if i<8
            plot(1:S,squeeze(RC_RMS_Seq(i,1:S,1)),...
                'LineWidth',1.2);
        else
            plot(1:S,squeeze(RC_RMS_Seq(i,1:S,1)),'--',...
                'LineWidth',1.2);
        end
        hold on

    end

    xlim([0,S+1])
    legend(lg_str,'box','off','NumColumns',2,'FontSize',font_Size-2)
    ylabel('Separation index','FontWeight','bold')
    xlabel('Spherical Slepian coefficients','FontWeight','bold')
    title('(a)')
    % title(Area)

    set(gca,'FontName','Times New Roman','FontSize',font_Size,'TickLength',[0.004,0.035]);

    axes(ha(2))
    plot(test_M_rec,WSeq12_13_23(:,1),'LineWidth',2,'color','k');
    hold on
    % p1=plot(M_rec,WSeq12_13_23(turn_WM_four(:,used_sta_M),1),'o','MarkerSize',8,...
    %     'MarkerFaceColor','r','MarkerEdgeColor','k');
    p1=plot(M_rec,WSeq12_13_23(test_M_rec_use,1),'o','MarkerSize',8,...
        'MarkerFaceColor','r','MarkerEdgeColor','k');

    xticks(test_M_rec)

    ylabel('Mean separation index','FontWeight','bold')
    xlabel('Width of Sliding window M (month)','FontWeight','bold')
    title('(b)')

    legend(p1,'Optimal M','box','off','Location','northwest','FontSize',font_Size)
    xlim([test_M_rec(1)-5,test_M_rec(end)+5])
    set(gca,'FontName','Times New Roman','FontSize',font_Size,'TickLength',[0.004,0.035]);

    tif_namefm=[Attach.Attach_ALL 'M_rec_Determine.tif'];
    print(fm,'-dtiff','-r300',fullfile(Attach.fig_path_ALL,tif_namefm));

%% plot for the optimal M (eigenvalues)

    font_Size=14;F_Position=[100,250,1000,450];
    F_gap=[.03 .08];F_marg_h=[.13 .05];F_marg_w=[.07 .02];

    fm=figure

    ha=tight_subplot(1,2,F_gap,F_marg_h,F_marg_w);

    set(gcf,'Position',F_Position)

    axes(ha(1))
    for i=1:test_M_num
        %     plot(1:N_ol,squeeze(RC_Seq(i,:,1))./N_rec);
        %     hold on
        %     WSeq12_13_23(i,1)=mean(squeeze(RC_Seq(i,:,1))./N_rec);

        if i<8
            plot(1:S,squeeze(RC_eva_Seq(i,1:S,1)),...
                'LineWidth',1.2);
        else
            plot(1:S,squeeze(RC_eva_Seq(i,1:S,1)),'--',...
                'LineWidth',1.2);
        end
        hold on

    end

    xlim([0,S+1])
    legend(lg_str,'box','off','NumColumns',2,'FontSize',font_Size-2)
    ylabel('Separation index','FontWeight','bold')
    xlabel('Spherical Slepian coefficients','FontWeight','bold')
    title('(a)')
    % title(Area)

    set(gca,'FontName','Times New Roman','FontSize',font_Size,'TickLength',[0.004,0.035]);

    axes(ha(2))
    plot(test_M_rec,ESeq12_13_23(:,1),'LineWidth',2,'color','k');
    hold on
    % p1=plot(M_rec,WSeq12_13_23(turn_WM_four(:,used_sta_M),1),'o','MarkerSize',8,...
    %     'MarkerFaceColor','r','MarkerEdgeColor','k');
    p1=plot(M_rec_E,ESeq12_13_23(test_M_rec_use_E,1),'o','MarkerSize',8,...
        'MarkerFaceColor','r','MarkerEdgeColor','k');

    xticks(test_M_rec)

    ylabel('Mean separation index','FontWeight','bold')
    xlabel('Width of Sliding window M (month)','FontWeight','bold')
    title('(b)')

    legend(p1,'Optimal M','box','off','Location','northwest','FontSize',font_Size)
    xlim([test_M_rec(1)-5,test_M_rec(end)+5])
    set(gca,'FontName','Times New Roman','FontSize',font_Size,'TickLength',[0.004,0.035]);

    tif_namefm=[Attach.Attach_ALL 'M_rec_Eigen_Determine.tif'];
    print(fm,'-dtiff','-r300',fullfile(Attach.fig_path_ALL,tif_namefm));

end

disp(['So the optimal M for reconstruction is determined to be : ' num2str(M_rec) ])

end