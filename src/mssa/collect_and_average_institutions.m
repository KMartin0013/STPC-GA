function [mssaSort, S, S_ol, SHC_degord] = collect_and_average_institutions(basicInfo, ...
    slepian_results, Attach)

Code_Version = basicInfo.Code_Version;
Area         = basicInfo.Area;
Lwindow      = basicInfo.Lwindow;
c11cmn       = basicInfo.c11cmn;
ddir1        = basicInfo.ddir1;      % output results directory (Results_*)
ddir2        = basicInfo.ddir2;      % output figure directory (Figures_*)
ddir3        = basicInfo.ddir3;      % output figure directory (Results_*/AddData)
ifilesRoot   = basicInfo.ifilesRoot; % root directory of original IFILES
saveAddData  = basicInfo.saveAddData;
redo         = basicInfo.redo;

use_institu  = basicInfo.use_institu;

intit_num    = length(slepian_results);


if ~redo && exist(fullfile(ddir3,['MainMSSA_ADD_' Attach.Attach_ALL '.mat']),'file')

    load(fullfile(ddir3,['MainMSSA_ADD_' Attach.Attach_ALL '.mat']),...
        'mssaSort','S','S_ol','SHC_degord')

    return

end

mssaSort=struct();
for ins=1:intit_num

    Dataproduct  = slepian_results(ins).Dataproduct;

    S            = slepian_results(ins).S;
    XY           = slepian_results(ins).XY;
    XY_ori       = slepian_results(ins).XY_ori;
    XY_buffer    = slepian_results(ins).XY_buffer;

    Slept_dates  = slepian_results(ins).date.Slept_dates;
    fill_dates   = slepian_results(ins).date.fill_dates;
    fill_months  = slepian_results(ins).date.fill_months;
    fill_nmonths = slepian_results(ins).date.fill_nmonths;
    use_nmonths  = slepian_results(ins).date.use_nmonths;
    use_months   = slepian_results(ins).date.use_months;
    data2use_monIndex   = slepian_results(ins).date.data2use_monIndex;
    missing_months      = slepian_results(ins).date.missing_months;

    functionintegrals   = slepian_results(ins).time.functionintegrals;

    %% Boundary output
    fid=fopen(fullfile(char(Attach.txt_path_each(ins)),'boundary_XY.txt'),'wt');
    fprintf(fid,'>\n');
    for i=1:size(XY,1)
        fprintf(fid,'%f\t%f\n',XY(i,1),XY(i,2));
    end
    fclose(fid);
    % fprintf(fid,'%f\t%f\n',XY(1,1),XY(1,2));

    figure
    plot(XY(:,1),XY(:,2))
    hold on
    plot(XY_ori(:,1),XY_ori(:,2))

    close all

    %% Calcualte GRACE

    %     save(fullfile(getenv('IFILES'),['Results_' parrent_path],['Main_' FIG_Attach '.mat']),'Table1','Table1_un')

    Dataprocess=Dataproduct;
    Dataprocess(1)=[];

    mssaSort(ins).name=char(use_institu(ins));
    mssaSort(ins).data_year_beg=slepian_results(ins).date.data_year_beg;
    mssaSort(ins).process=Dataprocess;
    mssaSort(ins).data_dates=slepian_results(ins).date.data_dates;
    mssaSort(ins).fill_dates=slepian_results(ins).date.fill_dates;
    mssaSort(ins).fill_months=slepian_results(ins).date.fill_months;
    mssaSort(ins).use_dates=slepian_results(ins).date.use_dates;
    mssaSort(ins).use_months=slepian_results(ins).date.use_months';
    mssaSort(ins).missing_dates=slepian_results(ins).date.missing_dates;
    mssaSort(ins).missing_months=slepian_results(ins).date.missing_months;

    %         mssa_Sort(ins).V=V;

    % basin-all
    mssaSort(ins).fill_deltacoffs=slepian_results(ins).fill.fill_deltacoffs;
    mssaSort(ins).fill_residdelta=slepian_results(ins).fill.fill_residdelta;
    mssaSort(ins).fill_signaldelta=slepian_results(ins).fill.fill_signaldelta;
    mssaSort(ins).fill_MASS=slepian_results(ins).fill.fill_MASS; % Gt
    mssaSort(ins).fill_MASSsig=slepian_results(ins).fill.fill_MASSsig; % Gt
    mssaSort(ins).fill_MASSres=slepian_results(ins).fill.fill_MASSres; % Gt

    % SSF coefficients
    mssaSort(ins).use_residcoffs=slepian_results(ins).time.use_residcoffs;
    mssaSort(ins).use_sleptdelta=slepian_results(ins).time.use_sleptdelta;

    mssaSort(ins).Fit.use_MASSfit=slepian_results(ins).time.use_MASSfit; % Gt
    mssaSort(ins).Fit.data_MASSslope=slepian_results(ins).time.data_MASSparams(2)*365.25; % Gt
    mssaSort(ins).Fit.data_MASSparamerrors = slepian_results(ins).time.data_MASSparamerrors(2); % Gt

    % this is the correct version of old A5 one.
    %         S=N;S_sig=N_sig; % we want to use S to replace the old 'N' to denote the number of SSF
    if strcmp(Dataproduct{4},'ocean')

        S_ol=S+1;

        mssaSort(ins).fill_IBdeltacoffs=slepian_results(ins).ocean.fill_IBdeltacoffs; % mmH2O (~kg/m2)
        mssaSort(ins).fill_IBresiddelta=slepian_results(ins).ocean.fill_IBresiddelta; % mmH2O
        mssaSort(ins).fill_IBsignaldelta=slepian_results(ins).ocean.fill_IBsignaldelta; % mmH2O

    else

        S_ol=S;

    end

    mssaSort(ins).fill_Grid_EWH=slepian_results(ins).grid.fill_Grid_EWH; % cm
    mssaSort(ins).fill_Grid_EWHsig=slepian_results(ins).grid.fill_Grid_EWHsig; % cm
    mssaSort(ins).fill_Grid_EWHres=slepian_results(ins).grid.fill_Grid_EWHres; % cm

    %     for SHC (actually it is not the main purpose here)
    disp('Calculate the SHC, it may take a while. please wait.')
    [potcoffs,~,~]=grace2plmt_m(Dataproduct{1},Dataproduct{2},Dataproduct{3},'SD', ...
        0,Dataproduct{4},Dataproduct{5});

    SHCpotcoffsdelta=potcoffs(1:size(potcoffs,1),:,:) - [repmat(mean(potcoffs(1:size(potcoffs,1),:,:),1),size(potcoffs,1),1)];

    SHCpotcoffsdelta(:,:,1:2)=potcoffs(:,:,1:2);

    [SHC_time,~]=size(SHCpotcoffsdelta);
    %         sample_SHC=squeeze(SHCpotcoffsdelta(1,:,:));

    % reoder for first fit
    SHC_reoderfit=[];
    for i=1:SHC_time
        SHC_reoderfit(i,:)=reshape(squeeze(SHCpotcoffsdelta(i,:,3:4)),[],1);
    end

    fitwhat=[3 365.25 365.25/2]; %linear trend and an annual and a semiannual harmonic term

    % only use the fit function in this code to get the residual of IB
    [data_SHCsignalcoffs,data_SHCresidcoffs,data_SHCftests,data_SHCextravalues]=slept2resid_m(SHC_reoderfit,...
        Slept_dates,fitwhat); % mm (values equivalent to kg/m^2)
    use_SHCsignalcoffs=data_SHCsignalcoffs(data2use_monIndex,:);
    use_SHCresidcoffs=data_SHCresidcoffs(data2use_monIndex,:);

    if any(missing_months)
        % complete estiamted signal with filled leakage values
        fill_SHCsignalcoffs_com=zeros(fill_nmonths,size(use_SHCsignalcoffs,2));
        fill_SHCresidcoffs_com=zeros(fill_nmonths,size(use_SHCresidcoffs,2));
        for i=1:size(use_SHCsignalcoffs,1)
            fill_SHCsignalcoffs_com(use_months(i),:)=use_SHCsignalcoffs(i,:);
            fill_SHCresidcoffs_com(use_months(i),:)=use_SHCresidcoffs(i,:);
        end

        for i=1:size(missing_months,2)
            fill_SHCsignalcoffs_com(missing_months(i),:)=data_SHCextravalues(i,:);
            %             SHCresid_com(leakage(i),:)=extravalues_resid(i,:);
            fill_SHCresidcoffs_com(missing_months(i),:)=0;
        end
        fill_SHCsignaldelta=fill_SHCsignalcoffs_com(1:fill_nmonths,:) - repmat(mean(use_SHCsignalcoffs(1:use_nmonths,:),1),fill_nmonths,1);
        fill_SHCresiddelta=fill_SHCresidcoffs_com(1:fill_nmonths,:) - repmat(mean(use_SHCresidcoffs(1:use_nmonths,:),1),fill_nmonths,1);

    else
        fill_SHCsignaldelta=use_SHCsignalcoffs(1:use_nmonths,:) - repmat(mean(use_SHCsignalcoffs(1:use_nmonths,:),1),use_nmonths,1);
        fill_SHCresiddelta=use_SHCresidcoffs(1:use_nmonths,:) - repmat(mean(use_SHCresidcoffs(1:use_nmonths,:),1),use_nmonths,1);
    end
    fill_SHCfilldelta=fill_SHCsignaldelta+fill_SHCresiddelta; % unit kg/m^2

    mssaSort(ins).fill_SHCfilldelta=fill_SHCfilldelta;
    mssaSort(ins).fill_SHCresiddelta=fill_SHCresiddelta;
    mssaSort(ins).fill_SHCsignaldelta=fill_SHCsignaldelta;

    if strcmp(Dataproduct{4},'ocean')
        % for GAD SHC (actually it is not the main purpose here)
        % this calculation is for subtracting from Gauer et al.,
        disp('Calculate the SHC for GAD, it may take a while. please wait.')
        [potcoffs_gsm,~,~]=grace2plmt_m(Dataproduct{1},Dataproduct{2},Dataproduct{3},'SD', ...
            0,'land',Dataproduct{5});

        % recover the GAD by removing the so-called EWH(GSM) from OBP
        % (GSM+GAD)
        SHCpotcoffsdelta_gad=potcoffs(1:size(potcoffs,1),:,:) - potcoffs_gsm(1:size(potcoffs_gsm,1),:,:) - ...
            [repmat(mean(potcoffs(1:size(potcoffs,1),:,:) - potcoffs_gsm(1:size(potcoffs_gsm,1),:,:),1),size(potcoffs_gsm,1),1)];

        SHCpotcoffsdelta_gad(:,:,1:2)=potcoffs_gsm(:,:,1:2);

        [SHC_time,~]=size(SHCpotcoffsdelta_gad);
        sample_SHC=squeeze(SHCpotcoffsdelta_gad(1,:,:));

        %             [r_GAD_Sample,lon,lat,Plm,degres]=plm2xyz(sample_SHC,1,c11cmn,60);

        %             figure
        %             imagesc(r_GAD_Sample)

        % reoder for first fit
        SHC_reoderfit=[];
        for i=1:SHC_time
            SHC_reoderfit(i,:)=reshape(squeeze(SHCpotcoffsdelta_gad(i,:,3:4)),[],1);
        end

        fitwhat=[3 365.25 365.25/2]; %linear trend and an annual and a semiannual harmonic term

        % only use the fit function in this code to get the residual of IB
        [data_SHCsignalcoffs,data_SHCresidcoffs,data_SHCftests,data_SHCextravalues]=slept2resid_m(SHC_reoderfit,...
            Slept_dates,fitwhat); % mm (values equivalent to kg/m^2)
        use_SHCsignalcoffs=data_SHCsignalcoffs(data2use_monIndex,:);
        use_SHCresidcoffs=data_SHCresidcoffs(data2use_monIndex,:);

        if any(missing_months)
            % complete estiamted signal with filled leakage values
            fill_SHCsignalcoffs_com=zeros(fill_nmonths,size(use_SHCsignalcoffs,2));
            fill_SHCresidcoffs_com=zeros(fill_nmonths,size(use_SHCresidcoffs,2));
            for i=1:size(use_SHCsignalcoffs,1)
                fill_SHCsignalcoffs_com(use_months(i),:)=use_SHCsignalcoffs(i,:);
                fill_SHCresidcoffs_com(use_months(i),:)=use_SHCresidcoffs(i,:);
            end

            for i=1:size(missing_months,2)
                fill_SHCsignalcoffs_com(missing_months(i),:)=data_SHCextravalues(i,:);
                %             SHCresid_com(leakage(i),:)=extravalues_resid(i,:);
                fill_SHCresidcoffs_com(missing_months(i),:)=0;
            end
            fill_SHCsignaldelta=fill_SHCsignalcoffs_com(1:fill_nmonths,:) - repmat(mean(use_SHCsignalcoffs(1:use_nmonths,:),1),fill_nmonths,1);
            fill_SHCresiddelta=fill_SHCresidcoffs_com(1:fill_nmonths,:) - repmat(mean(use_SHCresidcoffs(1:use_nmonths,:),1),fill_nmonths,1);

        else
            fill_SHCsignaldelta=use_SHCsignalcoffs(1:use_nmonths,:) - repmat(mean(use_SHCsignalcoffs(1:use_nmonths,:),1),use_nmonths,1);
            fill_SHCresiddelta=use_SHCresidcoffs(1:use_nmonths,:) - repmat(mean(use_SHCresidcoffs(1:use_nmonths,:),1),use_nmonths,1);
        end
        fill_SHCfilldelta=fill_SHCsignaldelta+fill_SHCresiddelta; % unit kg/m^2

        mssaSort(ins).fill_GADfilldelta=fill_SHCfilldelta;
        mssaSort(ins).fill_GADresiddelta=fill_SHCresiddelta;
        mssaSort(ins).fill_GADsignaldelta=fill_SHCsignaldelta;
    end

end

% search for the common months with real GRACE measurements among products
common_use_months=1:fill_nmonths;
for ins=1:intit_num
    common_use_months_index=ismember(mssaSort(ins).use_months,common_use_months);
    common_use_months=mssaSort(ins).use_months(common_use_months_index);
end

%
mssaSort(intit_num+1).name=char(use_institu(intit_num+1));
mssaSort(intit_num+1).data_year_beg=mssaSort(1).data_year_beg;
mssaSort(intit_num+1).process=Dataprocess;
mssaSort(intit_num+1).data_dates=[];
mssaSort(intit_num+1).fill_dates=fill_dates;
mssaSort(intit_num+1).fill_months=fill_months;
mssaSort(intit_num+1).use_dates=fill_dates(common_use_months); % default using common months % bug 1 fix
mssaSort(intit_num+1).use_months=common_use_months;  % default using common months % bug 1 fix
mssaSort(intit_num+1).missing_dates=fill_dates(setdiff(fill_months,fill_months(common_use_months))); % default using common months
mssaSort(intit_num+1).missing_months=fill_months(setdiff(fill_months,fill_months(common_use_months)))'; % default using common months


fill_deltacoffs=0;fill_residdelta=0;fill_signaldelta=0;
fill_MASS=0;fill_MASSsig=0;fill_MASSres=0;
use_sleptdelta=0;use_residcoffs=0;

fill_SHCfilldelta=0;fill_SHCresiddelta=0;fill_SHCsignaldelta=0;

fill_Grid_EWH=0;fill_Grid_EWHres=0;fill_Grid_EWHsig=0;
for ins=1:intit_num
    fill_deltacoffs=fill_deltacoffs+mssaSort(ins).fill_deltacoffs;
    fill_residdelta=fill_residdelta+mssaSort(ins).fill_residdelta;
    fill_signaldelta=fill_signaldelta+mssaSort(ins).fill_signaldelta;

    fill_MASS=fill_MASS+mssaSort(ins).fill_MASS;
    fill_MASSsig=fill_MASSsig+mssaSort(ins).fill_MASSsig;
    fill_MASSres=fill_MASSres+mssaSort(ins).fill_MASSres;

    fill_Grid_EWH=fill_Grid_EWH+mssaSort(ins).fill_Grid_EWH;
    fill_Grid_EWHres=fill_Grid_EWHres+mssaSort(ins).fill_Grid_EWHres;
    fill_Grid_EWHsig=fill_Grid_EWHsig+mssaSort(ins).fill_Grid_EWHsig;

    fill_SHCfilldelta=fill_SHCfilldelta+mssaSort(ins).fill_SHCfilldelta;
    fill_SHCresiddelta=fill_SHCresiddelta+mssaSort(ins).fill_SHCresiddelta;
    fill_SHCsignaldelta=fill_SHCsignaldelta+mssaSort(ins).fill_SHCsignaldelta;

    % transform the use month to be consistent with the first instution
    % (common months here)
    common_index=ismember(mssaSort(ins).use_months, common_use_months);

    use_sleptdelta=use_sleptdelta+mssaSort(ins).use_sleptdelta(common_index,:);
    use_residcoffs=use_residcoffs+mssaSort(ins).use_residcoffs(common_index,:);
end
mssaSort(intit_num+1).fill_deltacoffs=fill_deltacoffs/intit_num;
mssaSort(intit_num+1).fill_residdelta=fill_residdelta/intit_num;
mssaSort(intit_num+1).fill_signaldelta=fill_signaldelta/intit_num;

mssaSort(intit_num+1).fill_MASS=fill_MASS/intit_num; % Gt
mssaSort(intit_num+1).fill_MASSsig=fill_MASSsig/intit_num; % Gt
mssaSort(intit_num+1).fill_MASSres=fill_MASSres/intit_num; % Gt
%     mssa_Sort(intit_num+1).V=V;

mssaSort(intit_num+1).fill_Grid_EWH=fill_Grid_EWH/intit_num; % cm
mssaSort(intit_num+1).fill_Grid_EWHres=fill_Grid_EWHres/intit_num; % cm
mssaSort(intit_num+1).fill_Grid_EWHsig=fill_Grid_EWHsig/intit_num; % cm

mssaSort(intit_num+1).fill_SHCfilldelta=fill_SHCfilldelta/intit_num;
mssaSort(intit_num+1).fill_SHCresiddelta=fill_SHCresiddelta/intit_num;
mssaSort(intit_num+1).fill_SHCsignaldelta=fill_SHCsignaldelta/intit_num;

mssaSort(intit_num+1).use_sleptdelta=use_sleptdelta/intit_num;
mssaSort(intit_num+1).use_residcoffs=use_residcoffs/intit_num;

if strcmp(Dataproduct{4},'ocean')
    fill_IBdeltacoffs=0;fill_IBresiddelta=0;fill_IBsignaldelta=0;
    for ins=1:intit_num
        fill_IBdeltacoffs=fill_IBdeltacoffs+mssaSort(ins).fill_IBdeltacoffs;
        fill_IBresiddelta=fill_IBresiddelta+mssaSort(ins).fill_IBresiddelta;
        fill_IBsignaldelta=fill_IBsignaldelta+mssaSort(ins).fill_IBsignaldelta;
    end
    mssaSort(intit_num+1).fill_IBdeltacoffs=fill_IBdeltacoffs/intit_num;
    mssaSort(intit_num+1).fill_IBresiddelta=fill_IBresiddelta/intit_num;
    mssaSort(intit_num+1).fill_IBsignaldelta=fill_IBsignaldelta/intit_num;

    fill_GADfilldelta=0;fill_GADresiddelta=0;fill_GADsignaldelta=0;
    for ins=1:intit_num
        fill_GADfilldelta=fill_GADfilldelta+mssaSort(ins).fill_GADfilldelta;
        fill_GADresiddelta=fill_GADresiddelta+mssaSort(ins).fill_GADresiddelta;
        fill_GADsignaldelta=fill_GADsignaldelta+mssaSort(ins).fill_GADsignaldelta;
    end
    mssaSort(intit_num+1).fill_GADfilldelta=fill_GADfilldelta/intit_num;
    mssaSort(intit_num+1).fill_GADresiddelta=fill_GADresiddelta/intit_num;
    mssaSort(intit_num+1).fill_GADsignaldelta=fill_GADsignaldelta/intit_num;
end

% calculate the mean trend
[Cab_ave] = slepresid2cov(mssaSort(intit_num+1).use_residcoffs);
alphavarall_ave = functionintegrals*Cab_ave(1:S,1:S)*functionintegrals';

[use_MASSave_CC,~,~,~,~,~]=integral_fit(S,mssaSort(intit_num+1).use_dates,functionintegrals, ...
    mssaSort(intit_num+1).use_sleptdelta,Cab_ave);

use_MASSave=sum(use_MASSave_CC);

[use_ave_fit,use_ave_delta,use_ave_params,use_ave_paramerrors] = timeseriesfit([mssaSort(intit_num+1).use_dates' use_MASSave'],...
    alphavarall_ave,1,1);
% Make a matrix for the line, and 95% confidence in the fit
use_MASSavefit = [mssaSort(intit_num+1).use_dates' use_ave_fit use_ave_delta];
data_MASSaveparamerrors = use_ave_paramerrors*365.25;

mssaSort(intit_num+1).Fit.use_MASSfit=use_MASSavefit; % Gt
mssaSort(intit_num+1).Fit.data_MASSslope=use_ave_params(2)*365.25;
mssaSort(intit_num+1).Fit.data_MASSparamerrors = data_MASSaveparamerrors(2);

% degree and order of SHC
SHC_degord = SHCpotcoffsdelta(1,:,1:2);

if saveAddData
    save(fullfile(ddir3,['MainMSSA_ADD_' Attach.Attach_ALL '.mat']), ...
        'mssaSort','S','S_ol','SHC_degord')
end


end