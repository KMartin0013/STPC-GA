% This code is to fill the gap of the time series based on
% Multiple Singular Spectrum Analysis (M-SSA).
% Code version: 1
% Code Note:
%   Main function:
%       (1) fill the gap of required time series based an iterative method.
%       (2) Sampling figures of the comparison between different iterative
%       steps.
%   Main features:
%       (1) Default maximum iterative steps is 100.
%
% Tested on 9.13.0.2049777 (R2022b)
% Last modified by zhongtian.ma@connect.polyu.hk, 12/18/2024

% Main Reference:
% Ref1: Gauer, L. M., Chanard, K., & Fleitout, L. (2023). Data‐driven gap filling and spatio‐temporal filtering of the GRACE and GRACE‐FO records. Journal of Geophysical Research: Solid Earth, 128(5), e2022JB025561.

function varargout = MSSA_gap_fitting_correct(mssa_Sort,M,N,fig_note)
% MSSA usage for multiple institutions
%[MSSAn_Total_CGJI_EST_reconst,MSSAn_RC,iter_num,chi,RCtest]=MSSA_gap_fitting(A5_CJGI_Sort,M,Mnum,fig_note);

addpath(fileparts(mfilename('fullpath')))
% M=numel(tt_EST)/2;
Method=2;
% Mnum=8;                    % we want 8 RC
% smo_data=1;
emperi=0.9;

if exist('fig_note','var')
    fig=1;
else
    fig=0;
end

intit_num=numel(mssa_Sort);

%% for the total time series
% first M-SSA
for i=1:intit_num
    MSLA1_CJGI(:,i)=mssa_Sort(i).MSSA_TS';
end

time_num=size(MSLA1_CJGI,1);
fill_months=1:time_num;
data_year_beg=mssa_Sort.data_year_beg;

% [MSSA1_evalues,MSSA1_st_eofs,MSSA1_st_pcs,MSSA1_RC]=MSSA(1:time_num,MSLA1_CJGI,M,Method,N);

opts = struct();
opts.useSvds = true;
opts.returnFullEvalues = false;
opts.verbose = false;

[MSSA1_evalues,MSSA1_st_eofs,MSSA1_st_pcs,MSSA1_RC]=MSSA_fast(1:time_num,MSLA1_CJGI,M,Method,N, opts);

MSSA1_CJGIave=sum(sum(MSSA1_RC,3),2)/intit_num;

% std_MSSA1=std(MSSA1_CJGIave);
std_MSSA1=std(detrend(MSSA1_CJGIave));

% PLOT
font_Size=18;
F_gap=[.07 .03];F_marg_h=[.09 .04];F_marg_w=[.04 .02];
F_Position=[100,50,1500,900];
MarkSize=2;
col_gap=[150 150 150]/255;
col_est=[153 204 255]/255;col_out=[128 0 128]/255;
col_art=[238 173 14]/255;
abc='abcdefghijklmnopqrstuvwxyz';

if fig
    fl2=figure
    ha=tight_subplot(intit_num,1,F_gap,F_marg_h,F_marg_w);
end

Total_CGJI_fill_reconst_MSSA1=zeros(time_num,intit_num);
for i=1:intit_num

    tt_use=floor(mssa_Sort(i).use_months/12)+data_year_beg+mod(mssa_Sort(i).use_months,12)/12-1/24;
    tt_fil=floor(fill_months/12)+data_year_beg+mod(fill_months,12)/12-1/24;
    tt_lea=floor(mssa_Sort(i).missing_months/12)+data_year_beg+mod(mssa_Sort(i).missing_months,12)/12-1/24;

    missing_months=mssa_Sort(i).missing_months;
    MSSA_TS=mssa_Sort(i).MSSA_TS;

    % find the outliers
    outlier_down=find(MSSA_TS'<(MSSA1_CJGIave-3*std_MSSA1)==1);
    outlier_up=find(MSSA_TS'>(MSSA1_CJGIave+3*std_MSSA1)==1);
    outlier_months=[outlier_down;outlier_up];

    %     replace the outliers for the next M-SSA
    Total_fill_outfree=MSSA_TS';
    Total_fill_outfree(outlier_months)=MSSA1_CJGIave(outlier_months);
    Total_CGJI_fill_reconst_MSSA1(:,i)=Total_fill_outfree;

    if fig
        ybou1=[-20,20];ybou2=[-20,20];

        axes(ha(i));

        ybou1(2)=max(ybou1(2),max(MSSA_TS));
        ybou1(1)=min(ybou1(1),min(MSSA_TS));

        Cap_posix1=tt_use(5);Cap_posiy1=ybou1(2)-(ybou1(2)-ybou1(1))/10;
        Cap_posix2=tt_use(5);Cap_posiy2=ybou2(2)-(ybou2(2)-ybou2(1))/10;

        % plot the standard deviations of mean values for the active centers
        baseline=ybou1(2);

        he_std=area(tt_fil,[MSSA1_CJGIave-3*std_MSSA1 6*std_MSSA1*ones(numel(tt_fil),1)],baseline);
        he_std(1).FaceColor = [1 1 1];
        he_std(1).EdgeColor = col_est;

        he_std(2).FaceColor = col_est;
        he_std(2).EdgeColor = col_est;
        he_std(2).FaceAlpha = 0.4;
        he_std(2).EdgeAlpha = 0.7;

        hold on

        % plot the MSL before the replacement of outliers
        p1=plot(tt_fil,MSSA_TS,'color','black','linewidth',2);
        hold on

        if ~isempty(outlier_months)

            tt_out=tt_fil(outlier_months);

            for j=1:numel(tt_out)

                if outlier_months(j)==1
                    outlier_tt_neibour=[tt_out(j),tt_out(j)+1/24];
                    outlier_data_neibour=[Total_fill_outfree(outlier_months(j)), ...
                        (Total_fill_outfree(outlier_months(j))+Total_fill_outfree(outlier_months(j)+1))/2];

                elseif outlier_months(j)==numel(tt_fil)
                    outlier_tt_neibour=[tt_out(j)-1/24,tt_out(j)];
                    outlier_data_neibour=[(Total_fill_outfree(outlier_months(j)-1)+Total_fill_outfree(outlier_months(j)))/2,...
                        Total_fill_outfree(outlier_months(j))];

                else
                    outlier_tt_neibour=[tt_out(j)-1/24,tt_out(j),tt_out(j)+1/24];
                    outlier_data_neibour=[(Total_fill_outfree(outlier_months(j)-1)+Total_fill_outfree(outlier_months(j)))/2,...
                        Total_fill_outfree(outlier_months(j)), ...
                        (Total_fill_outfree(outlier_months(j))+Total_fill_outfree(outlier_months(j)+1))/2];

                end

%                 p4=plot(outlier_tt_neibour,outlier_MSL_neibour,'color',col_out,'linewidth',2);

                he_out=area(outlier_tt_neibour,ones(numel(outlier_tt_neibour),1)*ybou1(2),ybou1(1));
                he_out.FaceColor = col_out;
                he_out.EdgeColor = col_out;
                he_out.FaceAlpha = 0.4;
                he_out.EdgeAlpha = 0.7;
                hold on
            end

        end

        for j=1:numel(tt_lea)
            leak_tt_neibour=[tt_lea(j)-1/24,tt_lea(j),tt_lea(j)+1/24];
            leak_MSL_neibour=[(MSSA_TS(missing_months(j)-1)+MSSA_TS(missing_months(j)))/2,...
                MSSA_TS(missing_months(j)), ...
                (MSSA_TS(missing_months(j))+MSSA_TS(missing_months(j)+1))/2];
            p2=plot(leak_tt_neibour,leak_MSL_neibour,'r','linewidth',2);

            he_gap=area(leak_tt_neibour([1,3]),[ybou1(2);ybou1(2)],ybou1(1));
            he_gap.FaceColor = col_gap;
            he_gap.EdgeColor = col_gap;
            he_gap.FaceAlpha = 0.4;
            he_gap.EdgeAlpha = 0.7;

            hold on
        end

        p3=plot(tt_fil,MSSA1_CJGIave,'blue','linewidth',1);

        if ischar(mssa_Sort(i).MSSA_TS_order)
            title([mssa_Sort(i).name(1:end-4) ' inverted barometer'],'FontWeight','bold')
        else
            title([mssa_Sort(i).name(1:end-4) ' spherical Slepian coefficient (\alpha = ' num2str(mssa_Sort(i).MSSA_TS_order) ')'],'FontWeight','bold')
        end
%         title(['Slepian expansion coefficient (' A5_CJGI_Sort(i).name(1:end-4) ')'],'FontWeight','bold')

        text(Cap_posix1,Cap_posiy1,['(' abc(i) ')'],'FontSize',font_Size-1);

        %     if i==1
%         ylabel('Mass sea level (cm)','FontWeight','bold','FontSize',font_Size)
        %     end

        ylabel('(cm)')
        xlim([tt_fil(1)-0.5, tt_fil(end)+0.5])
        ylim(ybou1)
        xticks(floor(tt_fil(1)):1:floor(tt_fil(end)))
        set(gca, 'xticklabel', get(gca, 'xtick'),'FontSize',font_Size-2)

        set(gca,'FontName','Times New Roman','TickLength',[0.004,0.035]);

    end

end

if fig

    if exist('he_out','var')
        l1=legend([p1,p2,p3,he_gap,he_out],{'Raw data','Fit data','1st MSSA average', ...
            'GRACE gaps','Outliers'},'NumColumns',5)
    else
        l1=legend([p1,p2,p3,he_gap],{'Raw data','Fit data','1st MSSA data','GRACE gaps',},'NumColumns',4)
    end

    set(l1,'Position',[0.25,0.02,0.5,0.03],'box','off','fontsize',font_Size)

    set(gcf,'Position',F_Position)

    print(fl2,'-dtiff','-r125',[fig_note '_MSSA1st.tif']);

end

%% second step (iterative process)

% M=numel(tt_EST)/2;
Method=2;
% Mnum=8;                    % we want 8 RCs
% smo_data=1;

% maximum to iterate 100 times
MSSAn_Total_CGJI_fill_reconst=cell(100,1);

MSSAn_Total_CGJI_fill_reconst{1}=Total_CGJI_fill_reconst_MSSA1;
MSSAn_MSLA_CJGIave{1}=MSSA1_CJGIave;
MSSAn_std(1,:)=std(Total_CGJI_fill_reconst_MSSA1);

% the maximum number of iteration is 100
flag=zeros(intit_num,1);iter_num=zeros(intit_num,1);
for ii=1:100
    % iterative M-SSA

    MSLAn_CJGI=MSSAn_Total_CGJI_fill_reconst{ii};

%     [MSSAn_evalues,MSSAn_st_eofs,MSSAn_st_pcs,MSSAn_RC]=MSSA(tt_fil,MSLAn_CJGI,M,Method,N);

    opts = struct();
    opts.useSvds = true;
    opts.returnFullEvalues = false;
    opts.verbose = false;

    [MSSAn_evalues,MSSAn_st_eofs,MSSAn_st_pcs,MSSAn_RC]=MSSA_fast(tt_fil,MSLAn_CJGI,M,Method,N, opts);

    MSLAn_CJGI_reconst=squeeze(sum(MSSAn_RC,2));

    % ONLY institutions haven't meet the requirement of X will be iterated
    for i=find(flag==0)'

        missing_months=mssa_Sort(i).missing_months;

        % replace the gap
        MSLAn_c=MSLAn_CJGI(:,i);
        MSLAn_c(missing_months)=MSLAn_CJGI_reconst(missing_months,i);

        % Here is a correction for only considering the missing values
        % this is the old version
%         MSSAn_std(ii+1,i)=std(MSLAn_c);
% 
%         chi(ii+1,i)=sqrt( sum((MSLAn_CJGI(:,i)-MSLAn_c).^2)...
%             /MSSAn_std(ii+1,i)/MSSAn_std(ii,i) );

        % this is the corrected version
        MSSAn_std(ii+1,i)=std(MSLAn_c(missing_months));

        chi(ii+1,i)=sqrt( sum((MSLAn_CJGI(missing_months,i)-MSLAn_c(missing_months)).^2)...
            /MSSAn_std(ii+1,i)/MSSAn_std(ii,i) );


        if chi(ii+1,i)<0.1
            flag(i)=1;
        end

        MSLAn_CJGI(:,i)=MSLAn_c;

        iter_num(i)=iter_num(i)+1;
    end

    MSSAn_Total_CGJI_fill_reconst{ii+1}=MSLAn_CJGI;

    if sum(flag)==intit_num
        break
    end

end

% remove null cell
MSSAn_Total_CGJI_fill_reconst(max(iter_num)+2:end)=[];


Fs = 12;                    % Sampling frequency
T = 1/Fs;                     % Sampling period
L = length(tt_fil);                     % Length of signal
t = (0:L-1)*T;                % Time vector
freq = Fs*(0:(L/2))/L;

for i=1:N

    Yss1=zeros(size(MSSAn_RC(:,i,1)));
    Yss2=zeros(1,L/2+1);
    Yss3=zeros(1,L/2+1);
    for ii=1:intit_num

        Yss1=Yss1+MSSAn_RC(:,i,ii);

        Y = fft(MSSAn_RC(:,i,ii));
        P2 = abs(Y/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);

        for jj=1:L/2+1
            cdf_P1(jj) = sum(P1(1:jj))/sum(P1);
        end

        Yss2=Yss2+cdf_P1;
        Yss3=Yss3+P1';
    end

    hks(i) = kstest((Yss1-mean(Yss1))/std(Yss1),'Alpha',0.05);
    hli(i) = lillietest((Yss1-mean(Yss1))/std(Yss1),'Alpha',0.05);

    [a,b]=min(abs(1./freq-4/12));
    hcdf(i)=Yss2(b)>emperi*intit_num;

    [a,b]=max(Yss3);
    hcdf_fre(i)=freq(b);

end

%% PLOT
% col_gap=[255 99 71]/255;col_est=[153 204 255]/255;col_out=[128 0 128]/255;
col_gap=[150 150 150]/255;
col_est=[153 204 255]/255;col_out=[128 0 128]/255;
col_art=[238 173 14]/255;

% map_color=mymap("rainbow");
% font_Size=15;

if fig

    fl3=figure
    ha=tight_subplot(intit_num,1,F_gap,F_marg_h,F_marg_w);
    map_color=jet(max(iter_num));

    % Total_CGJI_EST_reconst_MSSA1=zeros(numel(tt_EST),intit_num);
    for i=1:intit_num

        ybou1=[-20,20];ybou2=[-20,20];

        tt_use=floor(mssa_Sort(i).use_months/12)+data_year_beg+mod(mssa_Sort(i).use_months,12)/12-1/24;
        tt_fil=floor(fill_months/12)+data_year_beg+mod(fill_months,12)/12-1/24;
        tt_lea=floor(mssa_Sort(i).missing_months/12)+data_year_beg+mod(mssa_Sort(i).missing_months,12)/12-1/24;

        missing_months=mssa_Sort(i).missing_months;
        MSSA_TS=mssa_Sort(i).MSSA_TS;

        axes(ha(i));
        map_colori=map_color(round(linspace(1,max(iter_num),iter_num(i))),:);

        for ii=1:iter_num(i)

            MSSAi_Total_CGJI_EST_reconst=MSSAn_Total_CGJI_fill_reconst{ii};

            p1=plot(tt_fil,MSSAi_Total_CGJI_EST_reconst(:,i),'color',map_colori(ii,:),'linewidth',2);
            hold on

            TickLabels_c{ii}=num2str(ii-1);

            ybou1(2)=max(ybou1(2),max(MSSAi_Total_CGJI_EST_reconst(:,i)));
            ybou1(1)=min(ybou1(1),min(MSSAi_Total_CGJI_EST_reconst(:,i)));
        end

        Cap_posix1=tt_use(5);Cap_posiy1=ybou1(2)-(ybou1(2)-ybou1(1))/10;
        Cap_posix2=tt_use(5);Cap_posiy2=ybou2(2)-(ybou2(2)-ybou2(1))/10;

        colormap(ha(i),map_colori)
        if iter_num(i)<=5
            colorbar('Ticks',linspace(0,1,iter_num(i)),...
                'TickLabels',TickLabels_c)
        else
            colorbar('Ticks',1/iter_num(i)/2:1/iter_num(i)*5:1,...
                'TickLabels',TickLabels_c(1:5:end))
        end

        for j=1:numel(tt_lea)
            leak_tt_neibour=[tt_lea(j)-1/24,tt_lea(j),tt_lea(j)+1/24];
            %         leak_MSL_neibour=[MSSA_TS(missing_months(j)-1),...
            %             MSSA_TS(missing_months(j)), ...
            %             MSSA_TS(missing_months(j)+1)];
            %         p2=plot(leak_tt_neibour,leak_MSL_neibour,'r','linewidth',2);

            he_gap=area(leak_tt_neibour([1,3]),[ybou1(2);ybou1(2)],ybou1(1));
            he_gap.FaceColor = col_gap;
            he_gap.EdgeColor = col_gap;
            he_gap.FaceAlpha = 0.4;
            he_gap.EdgeAlpha = 0.7;

            hold on
        end

        %     p3=plot(tt_EST,MSSA1_MSLA_CJGIave,'black','linewidth',1);

        %     if i==1
%         ylabel('Mass sea level (cm)','FontWeight','bold','FontSize',font_Size)
        %     end

        text(Cap_posix1,Cap_posiy1,['(' abc(i) ') The number of iteration: ' num2str(iter_num(i)-1)],'FontSize',font_Size-1);

        ylabel('(cm)')
        xlim([tt_fil(1)-0.5, tt_fil(end)+0.5])
        ylim(ybou1)
        xticks(floor(tt_fil(1)):1:floor(tt_fil(end)))
        set(gca, 'xticklabel', get(gca, 'xtick'),'FontSize',font_Size-2)

        if ischar(mssa_Sort(i).MSSA_TS_order)
            title([mssa_Sort(i).name(1:end-4) ' inverted barometer'],'FontWeight','bold')
        else
            title([mssa_Sort(i).name(1:end-4) ' spherical Slepian coefficient (\alpha = ' num2str(mssa_Sort(i).MSSA_TS_order) ')'],'FontWeight','bold')
        end
%         title(['Slepian expansion coefficient (' A5_CJGI_Sort(i).name(1:end-4) ')'],'FontWeight','bold')

        set(gca,'FontName','Times New Roman','TickLength',[0.004,0.035]);

    end

    set(gcf,'Position',F_Position)
    % if exist('he_out','var')
    %     l1=legend([p1,he_gap,he_out],{'Raw data','Interpolated data','1st M-SSA average', ...
    %         'Outliers Correction','Missing gap','Outliers'},'NumColumns',6)
    % else
        l1=legend([he_gap],{'GRACE Gap'},'NumColumns',4,'box','off')
    % end

    set(l1,'Position',[0.25,0.02,0.5,0.03],'box','off','fontsize',font_Size)

    set(gcf,'Position',F_Position)

    %     tif_name1=['AFig3v1_',FIG_Attach_ALL,'_Time_MSSAN.tif'];
    print(fl3,'-dtiff','-r125',[fig_note '_MSSAN.tif']);

end

% Collect output
varns={MSSAn_Total_CGJI_fill_reconst,MSSAn_RC,iter_num,chi,[hks;hli;hcdf;hcdf_fre]};
varargout=varns(1:nargout);


