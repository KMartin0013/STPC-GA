function [turn_MSSA_four, MSSA_evalues_sumup, MSSA_evalues] = det_N_TP(mssa_results, ...
    N_bou, Turning_number, used_sta_SN, Attach, plotProcess);

intit_num           = mssa_results.intit_num;
mssa_Sort           = mssa_results.mssa_Sort;
fillM_deltacoffs    = mssa_results.fillM_deltacoffs;
S_ol                = mssa_results.S_ol;
S                   = mssa_results.S;
N_gap               = mssa_results.N_gap;
M_rec               = mssa_results.M_rec;

threshold_N=[N_bou, 1-N_bou];

% for MSSA
% Lilliefors can legitimately return a p-value outside the tabulated
% range.  The value is still usable for the later binary significance
% decision, so silence only these two known boundary notices locally and
% restore the caller's warning state automatically.
warningStateLow = warning('query', 'stats:lillietest:OutOfRangePLow');
warningStateHigh = warning('query', 'stats:lillietest:OutOfRangePHigh');
warning('off', 'stats:lillietest:OutOfRangePLow');
warning('off', 'stats:lillietest:OutOfRangePHigh');
warningCleanup = onCleanup(@() local_restore_lillietest_warnings( ...
    warningStateLow, warningStateHigh)); %#ok<NASGU>

mssa_Sort_N=mssa_Sort(1:intit_num);

MSSA_evalues=zeros(intit_num*M_rec,S_ol);
for ss=1:S_ol
    [~,~,MSSA_evalue,~]=MSSA_final_noCDF_freqsort(mssa_Sort_N, ...
        fillM_deltacoffs(:,:,ss),M_rec,N_gap); %error Mssa1 should be Mssa2 to be consistent with the calculation of M*4

    MSSA_evalues(:,ss)=MSSA_evalue;
end

MSSA_evalues_sumup=zeros(size(MSSA_evalues));
for i=1:size(MSSA_evalues,1)
    MSSA_evalues_sumup(i,:)=sum(MSSA_evalues(1:i,:),1);
end

%% turning point of the sum of eigenvalues of MSSA (MSSA part)
option_sta=["mean","rms","std","linear"];

used_sta_MSSA=used_sta_SN; % what kind of turning point do you want? (default is RSS)
MaxNumChanges_MSSA=Turning_number; % how many turning point do you want?

turn_MSSA_four={};
for j=1:size(MSSA_evalues_sumup,2)
    turn_MSSA_four{j}=nan(MaxNumChanges_MSSA,4);
    for i=1:4
        used_sta=char(option_sta(i));
        %         turn_MSSA_four{j}=nan(max_turning,4);

        use_N_TP=find(MSSA_evalues_sumup(:,j)>threshold_N(1) & ...
            MSSA_evalues_sumup(:,j)<threshold_N(2));
        if isempty(use_N_TP)
            error('V3:NoNTurningPointDomain', ...
                ['No cumulative MSSA eigenvalues for coefficient %d fall ', ...
                 'inside the configured N_bound.'], j);
        end

        % whole SCS
        [turn_MSSA] = findchangepts(MSSA_evalues_sumup(use_N_TP,j), ...
            'Statistic', used_sta, 'MaxNumChanges', MaxNumChanges_MSSA);

        leng_MSSA=length(turn_MSSA);

        if leng_MSSA == 0
            error('V3:NoNTurningPoints', ...
                ['No N turning point was found for coefficient %d. ', ...
                 'Reduce config.turningNumber or adjust config.N_bound.'], j);
        end
        if leng_MSSA<MaxNumChanges_MSSA
            turn_MSSA(leng_MSSA+1:MaxNumChanges_MSSA)=turn_MSSA(leng_MSSA);
            warning('V3:FewerNTurningPoints', ...
                ['Only %d distinct N turning points were found for ', ...
                 'coefficient %d; the final candidate is repeated to ', ...
                 'keep turningNumber=%d.'], ...
                leng_MSSA, j, MaxNumChanges_MSSA)
        end

        turn_MSSA_four{j}(1:length(turn_MSSA),i)=use_N_TP(turn_MSSA);

    end
end

%% turning point MSSA
if plotProcess

    font_Size=14;F_Position=[100,250,700,500];

    f1b=figure

    max_mode=min(150,intit_num*M_rec);

    if N_bou~=0
        disp(['You choose the boundary restriction for N: ' num2str(N_bou)])
        he_gap1=area([0,max_mode],[N_bou 1-2*N_bou N_bou;N_bou 1-2*N_bou N_bou],0);
        he_gap1(2).FaceColor = [1,1,1];
        he_gap1(2).EdgeColor = [1,1,1];

        he_gap1(1).FaceColor = [0.4,0.4,0.4];
        he_gap1(1).EdgeColor = [0.4,0.4,0.4];
        he_gap1(1).FaceAlpha = 0.4;
        he_gap1(1).EdgeAlpha = 0.7;
        he_gap1(3).FaceColor = [0.4,0.4,0.4];
        he_gap1(3).EdgeColor = [0.4,0.4,0.4];
        he_gap1(3).FaceAlpha = 0.4;
        he_gap1(3).EdgeAlpha = 0.7;

        hold on
    end

    colorr=colormap('winter');
    dd=linspace(1,256,S);
    map_colori=colorr(floor(flipud(dd')),:);TickLabels_c={};
    for j=1:S
        %     color_line=repmat(0.2+0.02*j,1,3);
        color_line=map_colori(j,:);
        plot(1:max_mode,MSSA_evalues_sumup(1:max_mode,j),'color',color_line,'LineWidth',1);
        hold on

        p3=plot(turn_MSSA_four{j}(:,used_sta_MSSA), ...
            MSSA_evalues_sumup(turn_MSSA_four{j}(:,used_sta_MSSA),j), ...
            'o','MarkerSize',5,'MarkerFaceColor','y',...
            'MarkerEdgeColor','k');
        hold on

        xlim([0,max_mode])

        for i=1:MaxNumChanges_MSSA
            turning_MSSA(j,i)=turn_MSSA_four{j}(i,used_sta_MSSA);
            turning_MSSA_evalues(j,i)=MSSA_evalues_sumup( ...
                turn_MSSA_four{j}(i,used_sta_MSSA),j);
        end

        TickLabels_c{j}=num2str(j);
    end

    colormap(gca,map_colori)

    cb=colorbar('Ticks',[1/S/2:10/S:1,1],...
        'TickLabels',TickLabels_c([1:10:end,end]))

    for i=1:MaxNumChanges_MSSA
        plot(turning_MSSA(1:S,:),turning_MSSA_evalues(1:S,:),'k')
    end

    % cb = colorbar;
    cb.Position = [0.8 0.2 0.04,0.3];
    title(cb,'SSF','fontsize',14);

    ylim([0,1])

    ylabel({'\bf Sum of normalized eigenvalues ($\sum\widetilde{\lambda}$)','\bf of M-SSA'},...
        'Interpreter','latex','FontWeight','bold')
    xlabel('Reconstructed components (RC) of each SSF','FontWeight','bold')

    legend([p3],{'TP(N)'}, ...
        'NumColumns',2,'Position',[0.55,0.7,0.4,0.03],'FontSize',font_Size+2,'box','off')

    set(gca,'FontName','Times New Roman','FontSize',font_Size,'TickLength',[0.004,0.035]);

    xlim([0,max_mode])

    % set(gca, 'xticklabel', get(gca, 'xtick'))
    % xlim([0,Max_SSF])

    set(gcf,'Position',F_Position)

    tif_name1b=[Attach.Attach_ALL 'TurningPoint_MSSA.tif'];
    print(f1b,'-dtiff','-r300',fullfile(Attach.fig_path_ALL,tif_name1b));

end

end

function local_restore_lillietest_warnings(warningStateLow, warningStateHigh)
warning(warningStateLow);
warning(warningStateHigh);
end
