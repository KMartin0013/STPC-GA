function [turn_V_four] = det_S_TP(V, Max_S, S_bou, ...
    Turning_number, used_sta_SN, Attach, plotProcess);

threshold_S=[S_bou, 1-S_bou]; % restriction on extremely high or low lambda

% turning point of V (Slepian part)

disp('Now we determine different S by searching for breaking points of the distribution of lambda.')

% so we choose 'mean' as the chosen one
option_sta=["mean","rms","std","linear"];

% use_S_TP=find(V>0.1 & V<0.9);
use_S_TP=find(V>threshold_S(1) & V<threshold_S(2));
if isempty(use_S_TP)
    error('V3:NoSTurningPointDomain', ...
        'No Slepian eigenvalues fall inside the configured S_bound.');
end

used_sta_V=used_sta_SN; % what kind of turning point do you want? (default is RSS)
MaxNumChanges_V=Turning_number; % how many turning point do you want?

turn_V_four=nan(MaxNumChanges_V,4);
for i=1:4
    used_sta=char(option_sta(i));
    % whole SCS
    [turn_V] = findchangepts(V(use_S_TP), 'Statistic', used_sta, ...
        'MaxNumChanges', MaxNumChanges_V);

    leng_V=length(turn_V);

    if leng_V == 0
        error('V3:NoSTurningPoints', ...
            ['No S turning point was found. Reduce config.turningNumber ', ...
             'or adjust config.S_bound.']);
    end
    if leng_V<MaxNumChanges_V
        turn_V(leng_V+1:MaxNumChanges_V)=turn_V(leng_V);
        warning('V3:FewerSTurningPoints', ...
            ['Only %d distinct S turning points were found; the final ', ...
             'candidate is repeated to keep turningNumber=%d.'], ...
            leng_V, MaxNumChanges_V)
    end

    turn_V_four(1:length(turn_V),i)=use_S_TP(turn_V);

end

%%

if plotProcess
    font_Size=14;F_Position=[100,250,700,500];

    f1a=figure

    if S_bou~=0
        disp(['You choose the boundary restriction for S: ' num2str(S_bou)])
        he_gap1=area([0,Max_S],[S_bou 1-2*S_bou S_bou;S_bou 1-2*S_bou S_bou],0);
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

    p1=plot(1:Max_S,V(1:Max_S),'b','LineWidth',2);
    hold on

    p2=plot([0 turn_V_four(1,used_sta_V)],[mean(V(1:turn_V_four(1,used_sta_V))) mean(V(1:turn_V_four(1,used_sta_V)))],...
        'Color','k','LineWidth',1.5,'LineStyle','--');
    hold on
    for i=1:MaxNumChanges_V-1
        plot([turn_V_four(i,used_sta_V) turn_V_four(i,used_sta_V)],[0 1],...
            'Color',[0.8 0.8 0.8],'LineWidth',1)
        hold on
        plot([turn_V_four(i,used_sta_V) turn_V_four(i+1,used_sta_V)],[mean(V(turn_V_four(i,used_sta_V):turn_V_four(i+1,used_sta_V))) mean(V(turn_V_four(i,used_sta_V):turn_V_four(i+1,used_sta_V)))],...
            'Color','k','LineWidth',1.5,'LineStyle','--');
        hold on
    end
    plot([turn_V_four(MaxNumChanges_V,used_sta_V) turn_V_four(MaxNumChanges_V,used_sta_V)],[0 1],...
        'Color',[0.8 0.8 0.8],'LineWidth',1)
    hold on
    plot([turn_V_four(MaxNumChanges_V,used_sta_V) length(V)],[mean(V(turn_V_four(MaxNumChanges_V,used_sta_V):length(V))) mean(V(turn_V_four(MaxNumChanges_V,used_sta_V):length(V)))],...
        'Color','k','LineWidth',1.5,'LineStyle','--');
    hold on

    p3=plot(turn_V_four(:,used_sta_V),V(turn_V_four(:,used_sta_V)),'o','MarkerSize',5,...
        'MarkerFaceColor','r','MarkerEdgeColor','k');
    hold on

    xlim([0,Max_S])

    ylabel('Concentration ratio (\lambda)','FontWeight','bold')
    xlabel('Spherical Slepian functions (SSF)','FontWeight','bold')
    legend([p3],{'TP(S)'}, ...
        'NumColumns',2,'Position',[0.5,0.84,0.4,0.03],'FontSize',font_Size+2,'box','off')

    set(gca,'FontName','Times New Roman','FontSize',font_Size,'TickLength',[0.004,0.035]);

    % set(gca, 'xticklabel', get(gca, 'xtick'))
    xlim([0,Max_S])

    set(gcf,'Position',F_Position)

    tif_name1a=[Attach.Attach_ALL 'TurningPoint_V.tif'];
    print(f1a,'-dtiff','-r300',fullfile(Attach.fig_path_ALL,tif_name1a));

end

end
