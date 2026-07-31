% This code is to reconstruct the time series based on
% Multiple Singular Spectrum Analysis (M-SSA) with additional noise detection procedure for the p values.
% Code version: 1
% Code Note:
%   Main function:
%       (1) decompose the temporal reconstructed component (RC) of each 
%       spherical Slepian coefficients, 
%       (2) test the significance of each RC using K-S test and Lilliefors
%       test
%   Main features:
%       (1) Sampling figures were also plotted if fig_note is not empty
%       (2) Default plotting RCs is 8.
%
% Tested on 9.13.0.2049777 (R2022b)
% Last modified by zhongtian.ma@connect.polyu.hk, 12/18/2024

function varargout = MSSA_final_noCDF_freqsort(mssa_Sort,mssa_TS,Mnum,Nc,fig_note,Sig)
% MSSA usage for multiple institutions
%[mssa_TS_reconst,MSSAf_RC,MSSAf_evalues,RCTest]=MSSA_gap_fitting(A5_CJGI_Sort,mssa_TS,M,Mnum,fig_note);

% M=numel(tt_EST)/2;
Method=2;
% Mnum=8;                    % we want 8 RC
% smo_data=1;
emperi=0.9;

if exist('fig_note','var') && ~isempty(fig_note)
    fig=1;
else
    fig=0;
end

if exist('Sig','var')
    Significance=Sig;
else
    Significance=0.05;
end


intit_num=size(mssa_TS,2);
time_num=size(mssa_TS,1);

%% for the total time series

fill_months=1:time_num;
data_year_beg=mssa_Sort.data_year_beg;

tt_fil=floor(fill_months/12)+data_year_beg+mod(fill_months,12)/12-1/24;

% [MSSAf_evalues,MSSAf_st_eofs,MSSAf_st_pcs,MSSAf_RC]=MSSA(tt_fil,mssa_TS,Mnum,Method,Nc);

opts = struct();
opts.useSvds = true;
opts.returnFullEvalues = true;
opts.verbose = false;

[MSSAf_evalues,MSSAf_st_eofs,MSSAf_st_pcs,MSSAf_RC] = ...
    MSSA_fast(tt_fil,mssa_TS,Mnum,Method,Nc,opts);

mssa_TS_reconst=squeeze(sum(MSSAf_RC,2));

Fs = 12;                    % Sampling frequency
T = 1/Fs;                     % Sampling period
L = length(tt_fil);                     % Length of signal
t = (0:L-1)*T;                % Time vector
freq = Fs*(0:(L/2))/L;

for i=1:Nc

    Yss1=zeros(size(MSSAf_RC(:,i,1)));
    Yss2=zeros(1,L/2+1);
    Yss3=zeros(1,L/2+1);
    for ii=1:intit_num

        Yss1=Yss1+MSSAf_RC(:,i,ii);

        Y = fft(MSSAf_RC(:,i,ii));
        P2 = abs(Y/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);

        for jj=1:L/2+1
            cdf_P1(jj) = sum(P1(1:jj))/sum(P1);
        end

        Yss2=Yss2+cdf_P1;
        Yss3=Yss3+P1';
    end

    [hks(i),pks(i)] = kstest((Yss1-mean(Yss1))/std(Yss1),'Alpha',Significance);
%     hli(i) = lillietest((Yss1-mean(Yss1))/std(Yss1),'Alpha',0.05);
%     hks(i) = kstest(Yss1,'Alpha',0.05);
    [hli(i),pli(i)] = lillietest(Yss1,'Alpha',Significance);

    [a,b]=min(abs(1./freq-4/12));
    hcdf(i)=Yss2(b)>emperi*intit_num;

    [a,b]=max(Yss3);
    hcdf_fre(i)=freq(b);

    [Yss3_sort,I] = sort(Yss3/intit_num,'descend');
    Sep12(i)=Yss3_sort(1)/Yss3_sort(2);
    Sep13(i)=Yss3_sort(1)/Yss3_sort(3);
    Sep23(i)=Yss3_sort(2)/Yss3_sort(3);
end

% additional test for the sum of annual, long-term trend, semi-annual and
% high-frequency noise.

%% PLOT
col_gap=[0 205 205]/255;col_est=[153 204 255]/255;col_out=[128 0 128]/255;
col_art=[238 173 14]/255;
% map_color=mymap("rainbow");
font_Size=16;
MarkSize=2;
F_Position=[100,50,1500,900];
abc='abcdefghijklmnopqrstuvwxyz';

colorr_blind=[0,114,178;204,121,167;0,158,115;230,159,0;213,94,0;240,228,66;86,180,233]/255; % Blind-friendly colors

if fig

    % if Nc is too large for plot, we only main the first 8th of them
    if Nc>8
        Nc=8;
    end
    
    fl4=figure
    F_gap=[.02 .03];F_marg_h=[.09 .04];F_marg_w=[.04 .3];
    ha=tight_subplot(Nc,1,F_gap,F_marg_h,F_marg_w);
    map_colori=jet(4);

    % Total_CGJI_EST_reconst_MSSA1=zeros(numel(tt_EST),intit_num);
    for i=1:Nc

        ybou1=[-1,1];ybou2=[-20,20];
        tt_use=floor(mssa_Sort(1).use_months/12)+data_year_beg+mod(mssa_Sort(1).use_months,12)/12-1/24;
        tt_fil=floor(fill_months/12)+data_year_beg+mod(fill_months,12)/12-1/24;
        tt_lea=floor(mssa_Sort(1).missing_months/12)+data_year_beg+mod(mssa_Sort(1).missing_months,12)/12-1/24;

        axes(ha(i));

        pn=[];pn_name={};
        for ii=1:intit_num

%             p1(ii)=plot(tt_fil,MSSAf_RC(:,i,ii),'color',map_colori(ii,:),'linewidth',2);
            p1(ii)=plot(tt_fil,MSSAf_RC(:,i,ii),'color',colorr_blind(ii,:),'linewidth',2);
            hold on

            pn=[pn p1(ii)];pn_name{ii}=mssa_Sort(ii).name(1:end-4);
        end

        Maxx=max(MSSAf_RC(:,i,1:intit_num),[],'all');
        Minn=min(MSSAf_RC(:,i,1:intit_num),[],'all');

        ybou1(1)=min(ybou1(1),Minn);
        ybou1(2)=max(ybou1(2),Maxx);

        Cap_posix1=tt_fil(5);Cap_posiy1=ybou1(2)-(ybou1(2)-ybou1(1))/10;
        Cap_posix2=tt_fil(5);Cap_posiy2=ybou2(2)-(ybou2(2)-ybou2(1))/10;

        for j=1:numel(tt_lea)
            leak_tt_neibour=[tt_lea(j)-1/24,tt_lea(j),tt_lea(j)+1/24];
            %         leak_MSL_neibour=[Total_MSL_EST_reconst(missing_months(j)-1),...
            %             Total_MSL_EST_reconst(missing_months(j)), ...
            %             Total_MSL_EST_reconst(missing_months(j)+1)];
            %         p2=plot(leak_tt_neibour,leak_MSL_neibour,'r','linewidth',2);

            he_gap=area(leak_tt_neibour([1,3]),[ybou1(2);ybou1(2)],ybou1(1));
            he_gap.FaceColor = col_gap;
            he_gap.EdgeColor = col_gap;
            he_gap.FaceAlpha = 0.4;
            he_gap.EdgeAlpha = 0.7;

            hold on
        end

        text(Cap_posix1,Cap_posiy1,['(' abc(i) ')'],'FontSize',font_Size-1);

        %     p3=plot(tt_EST,MSSA1_MSLA_CJGIave,'black','linewidth',1);

        %     if i==1
        ylabel(['RC-' num2str(i)],'FontWeight','bold','FontSize',font_Size)
        %     end

        xlim([tt_fil(1)-0.5, tt_fil(end)+0.5])
        ylim(ybou1)

        pn_name{intit_num+1}='GRACE gap';
        if i==Nc
            xticks(floor(tt_fil(1)):1:floor(tt_fil(end)))
            l1=legend([pn,he_gap],pn_name,'NumColumns',5,'box','off','Location','northwest')
        else
            xticks([])
        end

        if i==1
            if ischar(mssa_Sort(i).MSSA_TS_order)
                title(['M-SSA Reconstruction for inverted barometer'],'FontWeight','bold')
            else
                title(['M-SSA Reconstruction for spherical slepian coefficient (\alpha = ' num2str(mssa_Sort(i).MSSA_TS_order) ')'],'FontWeight','bold')
            end
%             title(['M-SSA Reconstruction'],'FontWeight','bold','FontSize',font_Size-1)
        end

        set(gca,'FontName','Times New Roman','FontSize',font_Size-2,'TickLength',[0.004,0.035]);

        set(gca, 'xticklabel', get(gca, 'xtick'))

%         title(['GRACE SSF MSSA-' num2str(i)],'FontWeight','bold')

    end

    set(l1,'Position',[0.1,0.01,0.5,0.03],'box','off','fontsize',font_Size+1)

    F_gap=[.02 .03];F_marg_h=[.09 .04];F_marg_w=[.74 .02];
    ha=tight_subplot(Nc,1,F_gap,F_marg_h,F_marg_w);
    map_colori=jet(4);

    Fs = 12;                    % Sampling frequency
    T = 1/Fs;                     % Sampling period
    L = length(tt_fil);                     % Length of signal
    t = (0:L-1)*T;                % Time vector
    freq = Fs*(0:(L/2))/L;

    ybou1=[0,0.5];ybou2=[0,1];
    % Total_CGJI_EST_reconst_MSSA1=zeros(numel(tt_EST),intit_num);
    for i=1:Nc
        ybou1=[0,0.5];ybou2=[0,1];

        axes(ha(i));

        Yss=zeros(size(MSSAf_RC(:,i,1)));
        for ii=1:intit_num

            Yss=Yss+MSSAf_RC(:,i,ii);

            Y = fft(MSSAf_RC(:,i,ii));
            P2 = abs(Y/L);
            P1 = P2(1:L/2+1);
            P1(2:end-1) = 2*P1(2:end-1);

%             semilogx(1./freq,P1,'color',map_colori(ii,:),'LineWidth',1.2,'Marker','o','MarkerSize',MarkSize);
            semilogx(1./freq,P1,'color',colorr_blind(ii,:),'LineWidth',1.2,'Marker','o','MarkerSize',MarkSize);
            hold on
            
            ybou1(2)=max(ybou1(2),max(P1));
        end
% 
%         hks(i) = kstest((Yss-mean(Yss))/std(Yss),'Alpha',0.05);
%         hli(i) = lillietest((Yss-mean(Yss))/std(Yss),'Alpha',0.05);

        xticks([3/12,6/12,1,4,10]);
        xticklabels({'3','6','12','48','120'})
        if i<Nc
            xticklabels([]);
        else
            xlabel('Period (month)')
        end

        ylim(ybou1)

        Cap_posiy1=ybou1(2)-(ybou1(2)-ybou1(1))/10;
        Cap_posiy2=ybou1(2)-(ybou1(2)-ybou1(1))/3;

        text(0.12,Cap_posiy1,['(' abc(i+Nc) ')'],'FontSize',font_Size-1);

        test_x=2;
        if hcdf_fre(i)<1
            test_x=0.3;
        end

        
        %         
%         if hks(i)
%             text(test_x,Cap_posiy1,['K-S: pass (5%)'],'FontSize',font_Size-3);
%         else
%             text(test_x,Cap_posiy1,['K-S: fail (5%)'],'FontSize',font_Size-3);
%         end
%         if hli(i)
%             text(test_x,Cap_posiy2,['Lilliefors: pass (5%)'],'FontSize',font_Size-3);
%         else
%             text(test_x,Cap_posiy2,['Lilliefors: fail (5%)'],'FontSize',font_Size-3);
%         end

        if pks(i)<0.05
            text(test_x,Cap_posiy1,['p(K-S)<0.05'],'FontSize',font_Size-3);
        elseif pks(i)<0.1
            text(test_x,Cap_posiy1,['p(K-S)<0.1'],'FontSize',font_Size-3);
        elseif pks(i)<0.3
            text(test_x,Cap_posiy1,['p(K-S)<0.3'],'FontSize',font_Size-3);
        else
            text(test_x,Cap_posiy1,['p(K-S)>0.3'],'FontSize',font_Size-3);
        end
        if pli(i)<0.05
            text(test_x,Cap_posiy2,['p(Lilliefors)<0.05'],'FontSize',font_Size-3);
        elseif pli(i)<0.1
            text(test_x,Cap_posiy2,['p(Lilliefors)<0.1'],'FontSize',font_Size-3);
        elseif pli(i)<0.3
            text(test_x,Cap_posiy2,['p(Lilliefors)<0.3'],'FontSize',font_Size-3);
        else
            text(test_x,Cap_posiy2,['p(Lilliefors)>0.3'],'FontSize',font_Size-3);
        end

        if i==1
        title(['Power Spectrum'],'FontWeight','bold','FontSize',font_Size-1)
        end

        set(gca,'FontName','Times New Roman','FontSize',font_Size-2,'TickLength',[0.007,0.035]);

    end

%     F_gap=[.04 .03];F_marg_h=[.09 .04];F_marg_w=[.80 .02];
%     ha=tight_subplot(Nc,1,F_gap,F_marg_h,F_marg_w);
%     map_colori=jet(4);
% 
%     Fs = 12;                    % Sampling frequency
%     T = 1/Fs;                     % Sampling period
%     L = length(tt_fil);                     % Length of signal
%     t = (0:L-1)*T;                % Time vector
%     freq = Fs*(0:(L/2))/L;
% 
%     ybou1=[0,0.5];ybou2=[0,1];
%     Cap_posiy2=ybou2(2)-(ybou2(2)-ybou2(1))/10;
%     % Total_CGJI_EST_reconst_MSSA1=zeros(numel(tt_EST),intit_num);
%     for i=1:Nc
% 
%         axes(ha(i));
% 
%         Yss=zeros(1,L/2+1);
%         for ii=1:intit_num
% 
%             Y = fft(MSSAf_RC(:,i,ii));
%             P2 = abs(Y/L);
%             P1 = P2(1:L/2+1);
%             P1(2:end-1) = 2*P1(2:end-1);
% 
%             for jj=1:L/2+1
%                 cdf_P1(jj) = sum(P1(1:jj))/sum(P1);
%             end
% 
%             Yss=Yss+cdf_P1;
% 
%             semilogx(1./freq,cdf_P1,'color',map_colori(ii,:),'LineWidth',1.2,'Marker','o','MarkerSize',MarkSize);
%             hold on
%         end
% 
% %         [a,b]=min(abs(1./freq-4/12));
% %         hcdf(i)=Yss(b)>emperi*4;
% 
%         semilogx(4/12,emperi,'Marker','+','MarkerSize',MarkSize+4,'MarkerEdgeColor','k');
%         hold on
% 
%         if i<Nc
%             xticklabels([]);
%         else
%             xlabel('Period (year)')
%         end
% 
%         if hcdf(i)
%             text(0.8,Cap_posiy2,['CDF(f<0.33)>' num2str(emperi*100) '%'],'FontSize',font_Size-2);
%         else
%             text(0.8,Cap_posiy2,['CDF(f<0.33)<' num2str(emperi*100) '%'],'FontSize',font_Size-2);
%         end
% 
%         ylim(ybou2)
%         set(gca,'Fontname','Times New Roman','FontSize',font_Size-2)
% %         text(2004.8,Cap_posiy1,['Number of iteration: ' num2str(iter_num(i)-1)],'FontSize',font_Size-2);
%         title(['CDF Test'],'FontWeight','bold','FontSize',font_Size-2)

%     end

    set(gcf,'Position',F_Position)
    % if exist('he_out','var')
    %     l1=legend([p1,he_gap,he_out],{'Raw data','Interpolated data','1st M-SSA average', ...
    %         'Outliers Correction','Missing gap','Outliers'},'NumColumns',6)
    % else
    %     l1=legend([p1,p2,p3,he_gap],{'Raw data','Interpolated data','1st M-SSA data','Missing gap'},'NumColumns',4)
    % end

    %     set(l1,'Position',[0.25,0.02,0.5,0.03],'box','off','fontsize',font_Size)

    %     tif_name1=['AFig3v1_',FIG_Attach_ALL,'_Time_MSSAN.tif'];
    print(fl4,'-dtiff','-r125',[fig_note '_MSSAN_RC.tif']);

end

% Collect output
varns={mssa_TS_reconst,MSSAf_RC,MSSAf_evalues,[hks;hli;hcdf;hcdf_fre;Sep12;Sep13;Sep23;pks;pli]};
varargout=varns(1:nargout);


