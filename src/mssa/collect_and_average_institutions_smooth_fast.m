function [mssaSort_Me] = collect_and_average_institutions_smooth_fast(basicInfo, ...
    slepian_results, method_order, Attach)

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

fill_dates   = slepian_results(1).date.fill_dates;
fill_months  = slepian_results(1).date.fill_months;
fill_nmonths = slepian_results(1).date.fill_nmonths;

CC           = slepian_results(1).CC;
XY           = slepian_results(1).XY;
XY_ori       = slepian_results(1).XY_ori;
BasinArea    = slepian_results(1).BasinArea;

Default_Method = ["None","Gaussian300km","Gaussian500km","DDK3","DDK4",...
    "DDK5","DDK6","DDK7"];

matFile1=fullfile(ddir1,['MainMSSA_SMO1_' Attach.Attach_ALL '.mat']);

matFile2=fullfile(ddir1,['MainMSSA_SMO2_' Attach.Attach_ALL '.mat']);

if ~redo && exist(matFile2,'file')

    disp(['Loading: ' matFile2 ])
    load(matFile2,'mssaSort_Me')

    return

end

% matFile2=fullfile(ddir1,['MainMSSA_SMO2_' Attach.Attach_ALL '.mat']);

% redo=false;
% if ~redo && exist(matFile2,'file')
%
%     disp(['Loading: ' matFile2 ])
%     load(matFile2,'mssaSort_Me')
%
%     return
%
% end

% redo=true;

if ~redo && exist(matFile1,'file')

    load(matFile1,'mssaSort','SHC_degord')

else
    mssaSort=struct();
    SHC_degord=[];
    cacheDirSMO1=fullfile(ddir3,'CacheSMO1');
    if saveAddData && ~exist(cacheDirSMO1,'dir')
        mkdir(cacheDirSMO1);
    end

    for ins=1:intit_num

        Dataproduct  = slepian_results(ins).Dataproduct;
        Slept_dates  = slepian_results(ins).date.Slept_dates;
        data_dates   = slepian_results(ins).date.data_dates;
        use_nmonths  = slepian_results(ins).date.use_nmonths;
        use_months   = slepian_results(ins).date.use_months;
        data2use_monIndex   = slepian_results(ins).date.data2use_monIndex;
        missing_months      = slepian_results(ins).date.missing_months;

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

        % this is the correct version of old A5 one.
        %         S=N;S_sig=N_sig; % we want to use S to replace the old 'N' to denote the number of SSF
        if strcmp(Dataproduct{4},'ocean')

            mssaSort(ins).fill_IBdeltacoffs=slepian_results(ins).ocean.fill_IBdeltacoffs; % mmH2O (~kg/m2)
            mssaSort(ins).fill_IBresiddelta=slepian_results(ins).ocean.fill_IBresiddelta; % mmH2O
            mssaSort(ins).fill_IBsignaldelta=slepian_results(ins).ocean.fill_IBsignaldelta; % mmH2O

        end

        %     mssaSort(ins).fill_Grid_EWH=slepian_results(ins).grid.fill_Grid_EWH; % cm
        %     mssaSort(ins).fill_Grid_EWHsig=slepian_results(ins).grid.fill_Grid_EWHsig; % cm
        %     mssaSort(ins).fill_Grid_EWHres=slepian_results(ins).grid.fill_Grid_EWHres; % cm

        instCacheTag=char(Attach.Attach_each(ins));
        if isfield(basicInfo,'note') && ~isempty(basicInfo.note)
            instCacheTag=erase(instCacheTag, char(basicInfo.note));
        end
        matFile1_ins=fullfile(cacheDirSMO1,['SMO1INST_' instCacheTag '.mat']);

        if ~redo && exist(matFile1_ins,'file')
            cacheData=load(matFile1_ins,'mssaSort_ins','SHC_degord');
            if isfield(cacheData,'mssaSort_ins')
                cacheFields=fieldnames(cacheData.mssaSort_ins);
                for cacheFieldSeq=1:numel(cacheFields)
                    mssaSort(ins).(cacheFields{cacheFieldSeq})= ...
                        cacheData.mssaSort_ins.(cacheFields{cacheFieldSeq});
                end
                if isempty(SHC_degord) && isfield(cacheData,'SHC_degord')
                    SHC_degord=cacheData.SHC_degord;
                end
                if cache_has_smooth_fields(cacheData.mssaSort_ins, method_order)
                    disp(['Loading institution SMO1 cache: ' matFile1_ins])
                    continue
                end
            end
        end

        %%%   for SHC   %%%
        disp('Calculate the SHC, it may take a while. Please wait.')
        %         [potcoffs,~,~]=grace2plmt_m(Dataproduct{1},Dataproduct{2},Dataproduct{3},'SD', ...
        %             0,Dataproduct{4},Dataproduct{5});

        % the DDK should be applied to geopotential (POT) but not surface
        % density (SD), so we get POT here
        % notice that we only process the GSM with GAC subtracted
        [potcoffs_NoGIA,~,~]=grace2plmt_m(Dataproduct{1},Dataproduct{2},Dataproduct{3},'POT', ...
            0,'land');

        % if you want to restore GAD, first do the smoothing, then add it
        % back.
        if strcmp(Dataproduct{4},'ocean')

            [potcoffs_obp_NoGIA,~,~]=grace2plmt_m(Dataproduct{1},Dataproduct{2},Dataproduct{3},'POT', ...
                0,'ocean');

            potcoffs_gad_NoGIA=potcoffs_NoGIA;

            potcoffs_gad_NoGIA(:,:,3:4)=potcoffs_obp_NoGIA(:,:,3:4)-potcoffs_NoGIA(:,:,3:4);

        end

        % The GIA data is only in the SD form, so we first smooth the POT
        % data, and then transform to SD and do GIA correction.
        % GIA removal (SD)
        [data_dates,GIAt]=correct4gia(data_dates,Dataproduct{5});

        % Process the active center list.
        defval('a',fralmanac('a_EGM96','Earth'))

        for me=method_order
            %         for me=3

            sdcoffs_data=zeros(size(potcoffs_NoGIA,1),size(potcoffs_NoGIA,2),4);

            for i=1:length(data_dates)
                lmcosi_data_raw(:,1:4)=squeeze(potcoffs_NoGIA(i,:,1:4));

                switch me-1
                    case 0
                        % this is original
                        lmcosi_pot_data=lmcosi_data_raw;
                    case 1
                        % this is Gaussian smooth (300 km)
                        lmcosi_pot_data=plm_Gausssmooth(lmcosi_data_raw,300);
                    case 2
                        % this is Gaussian smooth (500 km)
                        lmcosi_pot_data=plm_Gausssmooth(lmcosi_data_raw,500);
                    case 3
                        % this is DDK3
                        lmcosi_pot_data=plm_DDK(lmcosi_data_raw,3);
                    case 4
                        % this is DDK4
                        lmcosi_pot_data=plm_DDK(lmcosi_data_raw,4);
                    case 5
                        % this is DDK5
                        lmcosi_pot_data=plm_DDK(lmcosi_data_raw,5);
                    case 6
                        % this is DDK4
                        lmcosi_pot_data=plm_DDK(lmcosi_data_raw,6);
                    case 7
                        % this is DDK5
                        lmcosi_pot_data=plm_DDK(lmcosi_data_raw,7);
                end

                if strcmp(Dataproduct{4},'ocean')
                    lmcosi_data_GAD_raw(:,1:4)=squeeze(potcoffs_gad_NoGIA(i,:,1:4));

                    lmcosi_pot_data(:,3:4)=lmcosi_pot_data(:,3:4)+lmcosi_data_GAD_raw(:,3:4);
                end

                % change to surface density kg/m^2, with GIA removal
                lmcosi_sd_data=plm2pot([lmcosi_pot_data(:,1:2) lmcosi_pot_data(:,3:4)*a],[],[],[],4);

                lmcosi_sd_data(:,3:4)=lmcosi_sd_data(:,3:4) -  squeeze(GIAt(i,1:size(potcoffs_NoGIA,2),3:4)) ;

                sdcoffs_data(i,:,:)=lmcosi_sd_data;

            end

            SHCcoffsdelta=sdcoffs_data(1:size(sdcoffs_data,1),:,:) - [repmat(mean(sdcoffs_data(1:size(sdcoffs_data,1),:,:),1),size(sdcoffs_data,1),1)];

            SHCcoffsdelta(:,:,1:2)=sdcoffs_data(:,:,1:2);

            [SHC_time,SHC_coefficient_SHC_values]=size(SHCcoffsdelta);
            sample_SHC=squeeze(SHCcoffsdelta(1,:,:));

            % reoder for first fit
            SHC_reoderfit=[];
            for i=1:SHC_time
                SHC_reoderfit(i,:)=reshape(squeeze(SHCcoffsdelta(i,:,3:4)),[],1);
            end

            %this process is for comparison with cases when GAD is also
            %smoothed
            %         %     for SHC (actually it is not the main purpose here)
            %         disp('Calculate the SHC, it may take a while. please wait.')
            %         [potcoffs,~,~]=grace2plmt_m(Dataproduct{1},Dataproduct{2},Dataproduct{3},'SD', ...
            %             0,Dataproduct{4},Dataproduct{5});
            %
            %         SHCpotcoffsdelta=potcoffs(1:size(potcoffs,1),:,:) - [repmat(mean(potcoffs(1:size(potcoffs,1),:,:),1),size(potcoffs,1),1)];
            %
            %         SHCpotcoffsdelta(:,:,1:2)=potcoffs(:,:,1:2);
            %
            %         [SHC_time,~]=size(SHCpotcoffsdelta);
            %         %         sample_SHC=squeeze(SHCpotcoffsdelta(1,:,:));
            %
            %         % reoder for first fit
            %         SHC_reoderfit2=[];
            %         for i=1:SHC_time
            %             SHC_reoderfit2(i,:)=reshape(squeeze(SHCpotcoffsdelta(i,:,3:4)),[],1);
            %         end

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

                for i=1:numel(missing_months)
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

            switch me-1
                case 0
                    mssaSort(ins).fill_SHCfilldelta=fill_SHCfilldelta;
                    %                     mssaSort(ins).fill_SHCresiddelta=fill_SHCresiddelta;
                    %                     mssaSort(ins).fill_SHCsignaldelta=fill_SHCsignaldelta;
                case 1
                    mssaSort(ins).fill_SHCfilldelta_Gau300=fill_SHCfilldelta;
                    %                     mssaSort(ins).fill_SHCresiddelta_Gau300=fill_SHCresiddelta;
                    %                     mssaSort(ins).fill_SHCsignaldelta_Gau300=fill_SHCsignaldelta;
                case 2
                    mssaSort(ins).fill_SHCfilldelta_Gau500=fill_SHCfilldelta;
                    %                     mssaSort(ins).fill_SHCresiddelta_Gau500=fill_SHCresiddelta;
                    %                     mssaSort(ins).fill_SHCsignaldelta_Gau500=fill_SHCsignaldelta;
                case 3
                    mssaSort(ins).fill_SHCfilldelta_DDK3=fill_SHCfilldelta;
                    %                     mssaSort(ins).fill_SHCresiddelta_DDK3=fill_SHCresiddelta;
                    %                     mssaSort(ins).fill_SHCsignaldelta_DDK3=fill_SHCsignaldelta;
                case 4
                    mssaSort(ins).fill_SHCfilldelta_DDK4=fill_SHCfilldelta;
                    %                     mssaSort(ins).fill_SHCresiddelta_DDK4=fill_SHCresiddelta;
                    %                     mssaSort(ins).fill_SHCsignaldelta_DDK4=fill_SHCsignaldelta;
                case 5
                    mssaSort(ins).fill_SHCfilldelta_DDK5=fill_SHCfilldelta;
                    %                     mssaSort(ins).fill_SHCresiddelta_DDK5=fill_SHCresiddelta;
                    %                     mssaSort(ins).fill_SHCsignaldelta_DDK5=fill_SHCsignaldelta;
                case 6
                    mssaSort(ins).fill_SHCfilldelta_DDK6=fill_SHCfilldelta;
                    %                     mssaSort(ins).fill_SHCresiddelta_DDK6=fill_SHCresiddelta;
                    %                     mssaSort(ins).fill_SHCsignaldelta_DDK6=fill_SHCsignaldelta;
                case 7
                    mssaSort(ins).fill_SHCfilldelta_DDK7=fill_SHCfilldelta;
                    %                     mssaSort(ins).fill_SHCresiddelta_DDK7=fill_SHCresiddelta;
                    %                     mssaSort(ins).fill_SHCsignaldelta_DDK7=fill_SHCsignaldelta;
            end

        end

        if isempty(SHC_degord)
            SHC_degord = SHCcoffsdelta(1,:,1:2);
        end

        if saveAddData
            mssaSort_ins=mssaSort(ins);
            method_order_cached=method_order;
            save(matFile1_ins,'mssaSort_ins','SHC_degord','method_order_cached')
        end

    end

    % search for the common months with real GRACE measurements among products
    common_use_months=1:fill_nmonths;
    for ins=1:intit_num
        common_use_months_index=ismember(mssaSort(ins).use_months,common_use_months);
        common_use_months=mssaSort(ins).use_months(common_use_months_index);
    end

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
    end

    if isempty(SHC_degord) && exist('SHCcoffsdelta','var')
        SHC_degord = SHCcoffsdelta(1,:,1:2);
    end
    if isempty(SHC_degord)
        error('SHC_degord was not loaded or computed for SMO1 preprocessing.')
    end

    if saveAddData
        save(matFile1, 'mssaSort','SHC_degord')
    end

end
%% spherical harmonics -> spatial maps and area-weighted time series (MSSA gap filling)
% mssaSort_Me=mssaSort;

[~,N_SHC_half,~]=size(SHC_degord);
N_SHC=N_SHC_half*2;

% time wasting
[~,lon,lat,~,~]=plm2xyz(CC{1},1,c11cmn,60);
[lonlon,latlat]=meshgrid(lon,lat);
[in, on]=check_polygon_in(XY,lonlon,latlat);

% area-weighted should be used to more precisely calculate the regional
% basin-averaged variables
for i=1:length(lat)
    for j=1:length(lon)
        c11cmn_area(i,j)=areaquad(lat(i)-0.5,lon(j)-0.5,lat(i)+0.5,lon(j)+0.5);
    end
end

mssaSort_Me=struct();
mssaSort_Me.name=[char(use_institu(intit_num+1)) '_SMO'];
mssaSort_Me.c11cmn=c11cmn;
mssaSort_Me.lon=lon;
mssaSort_Me.lat=lat;
mssaSort_Me.process=mssaSort(intit_num+1).process;

fill_dates_vec=datevec(mssaSort(intit_num+1).fill_dates);
mssaSort_Me.fill_dates=fill_dates_vec(:,1:3);
use_dates_vec=datevec(mssaSort(intit_num+1).use_dates);
mssaSort_Me.use_dates=use_dates_vec(:,1:3);
missing_dates_vec=datevec(mssaSort(intit_num+1).missing_dates);
mssaSort_Me.missing_dates=missing_dates_vec(:,1:3);

mssaSort_Me.EWH.unit='cm';
mssaSort_Me.MASS.unit='Gt';

matFile2=fullfile(ddir1,['MainMSSA_SMO2_' Attach.Attach_ALL '.mat']);

M=floor(numel(fill_dates) / 2);
N=8;
for me_seq=1:numel(method_order)

    me=method_order(me_seq);

    switch me-1
        case 0
            shcField = 'fill_SHCfilldelta';
            shcField_res = 'fill_SHCresiddelta';
            shcField_sig = 'fill_SHCsignaldelta';
        case 1
            shcField = 'fill_SHCfilldelta_Gau300';
            shcField_res = 'fill_SHCresiddelta_Gau300';
            shcField_sig = 'fill_SHCsignaldelta_Gau300';
        case 2
            shcField = 'fill_SHCfilldelta_Gau500';
            shcField_res = 'fill_SHCresiddelta_Gau500';
            shcField_sig = 'fill_SHCsignaldelta_Gau500';
        case 3
            shcField = 'fill_SHCfilldelta_DDK3';
            shcField_res = 'fill_SHCresiddelta_DDK3';
            shcField_sig = 'fill_SHCsignaldelta_DDK3';
        case 4
            shcField = 'fill_SHCfilldelta_DDK4';
            shcField_res = 'fill_SHCresiddelta_DDK4';
            shcField_sig = 'fill_SHCsignaldelta_DDK4';
        case 5
            shcField = 'fill_SHCfilldelta_DDK5';
            shcField_res = 'fill_SHCresiddelta_DDK5';
            shcField_sig = 'fill_SHCsignaldelta_DDK5';
        case 6
            shcField = 'fill_SHCfilldelta_DDK6';
            shcField_res = 'fill_SHCresiddelta_DDK6';
            shcField_sig = 'fill_SHCsignaldelta_DDK6';
        case 7
            shcField = 'fill_SHCfilldelta_DDK7';
            shcField_res = 'fill_SHCresiddelta_DDK7';
            shcField_sig = 'fill_SHCsignaldelta_DDK7';
    end

    matFile2_me=fullfile(ddir3,['MainMSSA_SMO2_' char(Default_Method(me)) ...
        '_' Attach.Attach_ALL '.mat']);

    disp(['Do gap-filling for ', char(Default_Method(me)) ' filter']);

    % for EWH
    mssaSort_Gap=struct();
    mssaSort_Gap.name=Default_Method(me);
    mssaSort_Gap.use_months=mssaSort(intit_num+1).use_months;
    mssaSort_Gap.missing_months=mssaSort(intit_num+1).missing_months;
    mssaSort_Gap.data_year_beg=mssaSort(intit_num+1).data_year_beg;

    % expand to one struct per institution
    mssaSort_Gap = repmat(mssaSort_Gap, intit_num, 1);

    for ins = 1:intit_num
        mssaSort_Gap(ins).name = mssaSort(ins).name;
    end

    if ~redo && exist(matFile2_me,'file')

        disp(['Loading: ' matFile2_me ])

        load(matFile2_me, 'EWH_EST_MSSA', ...
            'fill_EWH_MSSA_SHC', 'fill_MASS_MSSA_SHC');

    else

        fill_SHCfilldelta=0;fill_SHCresiddelta=0;fill_SHCsignaldelta=0;

        for ins = 1:intit_num
            fill_SHCfilldelta = fill_SHCfilldelta + mssaSort(ins).(shcField);
%             fill_SHCresiddelta = fill_SHCresiddelta + mssaSort(ins).(shcField_res);
%             fill_SHCsignaldelta = fill_SHCsignaldelta + mssaSort(ins).(shcField_sig);
        end

        fill_SHCfilldelta = fill_SHCfilldelta / intit_num;
        mssaSort(intit_num+1).(shcField) = fill_SHCfilldelta;

%         switch me-1
%             case 0
% 
%                 for ins=1:intit_num
%                     fill_SHCfilldelta=fill_SHCfilldelta+mssaSort(ins).fill_SHCfilldelta;
%                     %                 fill_SHCresiddelta=fill_SHCresiddelta+mssaSort(ins).fill_SHCresiddelta;
%                     %                 fill_SHCsignaldelta=fill_SHCsignaldelta+mssaSort(ins).fill_SHCsignaldelta;
%                 end
%                 fill_SHCfilldelta=fill_SHCfilldelta/intit_num;
%                 mssaSort(intit_num+1).fill_SHCfilldelta=fill_SHCfilldelta;
% 
%             case 1
% 
%                 for ins=1:intit_num
%                     fill_SHCfilldelta=fill_SHCfilldelta+mssaSort(ins).fill_SHCfilldelta_Gau300;
%                     %                 fill_SHCresiddelta=fill_SHCresiddelta+mssaSort(ins).fill_SHCresiddelta_Gau300;
%                     %                 fill_SHCsignaldelta=fill_SHCsignaldelta+mssaSort(ins).fill_SHCsignaldelta_Gau300;
%                 end
%                 fill_SHCfilldelta=fill_SHCfilldelta/intit_num;
%                 mssaSort(intit_num+1).fill_SHCfilldelta_Gau300=fill_SHCfilldelta;
% 
%             case 2
% 
%                 for ins=1:intit_num
%                     fill_SHCfilldelta=fill_SHCfilldelta+mssaSort(ins).fill_SHCfilldelta_Gau500;
%                     %                 fill_SHCresiddelta=fill_SHCresiddelta+mssaSort(ins).fill_SHCresiddelta_Gau500;
%                     %                 fill_SHCsignaldelta=fill_SHCsignaldelta+mssaSort(ins).fill_SHCsignaldelta_Gau500;
%                 end
%                 fill_SHCfilldelta=fill_SHCfilldelta/intit_num;
%                 mssaSort(intit_num+1).fill_SHCfilldelta_Gau500=fill_SHCfilldelta;
% 
%             case 3
% 
%                 for ins=1:intit_num
%                     fill_SHCfilldelta=fill_SHCfilldelta+mssaSort(ins).fill_SHCfilldelta_DDK3;
%                     %                 fill_SHCresiddelta=fill_SHCresiddelta+mssaSort(ins).fill_SHCresiddelta_DDK3;
%                     %                 fill_SHCsignaldelta=fill_SHCsignaldelta+mssaSort(ins).fill_SHCsignaldelta_DDK3;
%                 end
%                 fill_SHCfilldelta=fill_SHCfilldelta/intit_num;
%                 mssaSort(intit_num+1).fill_SHCfilldelta_DDK3=fill_SHCfilldelta;
% 
%             case 4
% 
%                 for ins=1:intit_num
%                     fill_SHCfilldelta=fill_SHCfilldelta+mssaSort(ins).fill_SHCfilldelta_DDK4;
%                     %                 fill_SHCresiddelta=fill_SHCresiddelta+mssaSort(ins).fill_SHCresiddelta_DDK4;
%                     %                 fill_SHCsignaldelta=fill_SHCsignaldelta+mssaSort(ins).fill_SHCsignaldelta_DDK4;
%                 end
%                 fill_SHCfilldelta=fill_SHCfilldelta/intit_num;
%                 mssaSort(intit_num+1).fill_SHCfilldelta_DDK4=fill_SHCfilldelta;
% 
%             case 5
% 
%                 for ins=1:intit_num
%                     fill_SHCfilldelta=fill_SHCfilldelta+mssaSort(ins).fill_SHCfilldelta_DDK5;
%                     %                 fill_SHCresiddelta=fill_SHCresiddelta+mssaSort(ins).fill_SHCresiddelta_DDK5;
%                     %                 fill_SHCsignaldelta=fill_SHCsignaldelta+mssaSort(ins).fill_SHCsignaldelta_DDK5;
%                 end
%                 fill_SHCfilldelta=fill_SHCfilldelta/intit_num;
%                 mssaSort(intit_num+1).fill_SHCfilldelta_DDK5=fill_SHCfilldelta;
% 
%             case 6
% 
%                 for ins=1:intit_num
%                     fill_SHCfilldelta=fill_SHCfilldelta+mssaSort(ins).fill_SHCfilldelta_DDK6;
%                     %                 fill_SHCresiddelta=fill_SHCresiddelta+mssaSort(ins).fill_SHCresiddelta_DDK6;
%                     %                 fill_SHCsignaldelta=fill_SHCsignaldelta+mssaSort(ins).fill_SHCsignaldelta_DDK6;
%                 end
%                 fill_SHCfilldelta=fill_SHCfilldelta/intit_num;
%                 mssaSort(intit_num+1).fill_SHCfilldelta_DDK6=fill_SHCfilldelta;
% 
%             case 7
% 
%                 for ins=1:intit_num
%                     fill_SHCfilldelta=fill_SHCfilldelta+mssaSort(ins).fill_SHCfilldelta_DDK7;
%                     %                 fill_SHCresiddelta=fill_SHCresiddelta+mssaSort(ins).fill_SHCresiddelta_DDK7;
%                     %                 fill_SHCsignaldelta=fill_SHCsignaldelta+mssaSort(ins).fill_SHCsignaldelta_DDK7;
%                 end
%                 fill_SHCfilldelta=fill_SHCfilldelta/intit_num;
%                 mssaSort(intit_num+1).fill_SHCfilldelta_DDK7=fill_SHCfilldelta;
% 
%         end

        for ins=1:intit_num

            fill_SHCfilldelta_ins = mssaSort(ins).(shcField);

            % reoder back
            fill_SHCdeltacoffs=[];fill_SHCsignaldelta=[];fill_SHCresiddelta=[];
            for i=1:numel(fill_dates)
                fill_SHCdeltacoffs(i,:,1:2)=SHC_degord;
                fill_SHCdeltacoffs(i,:,3:4)=reshape(fill_SHCfilldelta_ins(i,1:N_SHC),[],2);

                %         fill_SHCsignaldelta(i,:,1:2)=SHC_degord;
                %         fill_SHCsignaldelta(i,:,3:4)=reshape(use_SHCsignaldelta(i,1:N_SHC),[],2);
                %
                %         fill_SHCresiddelta(i,:,1:2)=SHC_degord;
                %         fill_SHCresiddelta(i,:,3:4)=reshape(use_SHCresiddelta(i,1:N_SHC),[],2);

            end

            %         fill_EWH_SHC=zeros(1,numel(fill_dates));
            %         fill_Grid_EWH_SHC=zeros(numel(fill_dates),numel(lat),numel(lon));

            for i=1:numel(fill_dates)
                TWSlmcosi=squeeze(fill_SHCdeltacoffs(i,:,:));
                [fill_Grid_SHC_mon,~,~,~,~]=plm2xyz(TWSlmcosi,1,c11cmn,Lwindow);
                r_ewh_each=fill_Grid_SHC_mon/1000*100; % change from kg/m^2 to cm (/ 1000 kg/m*3 * 100)
                mssaSort(ins).fill_EWH_SHC(i)=sum(r_ewh_each(in).*c11cmn_area(in))/sum(c11cmn_area(in)); % cm

                mssaSort(ins).fill_Grid_EWH_SHC(i,:,:)=r_ewh_each;

            end

            if strcmp(mssaSort_Me.process{3},'ocean')

%                 fill_MSL_SHC=zeros(1,numel(fill_dates));
%                 fill_Grid_MSL_SHC=zeros(numel(fill_dates),numel(lat),numel(lon));

                for i=1:numel(fill_dates)
                    mssaSort(ins).fill_Grid_MSL_SHC(i,:,:)=mssaSort(ins).fill_Grid_EWH_SHC(i,:,:)*1000/1028-...
                        mssaSort(ins).fill_IBdeltacoffs(i)'*1000/1028/10; % change unit from mmH2O to cm(salty H2O)

                end

                mssaSort(ins).fill_MSL_SHC = mssaSort(ins).fill_EWH_SHC*1000/1028-...
                    mssaSort(ins).fill_IBdeltacoffs'*1000/1028/10;%from mmH20 to cm

            end


        end

        %% gap-filling

%         % for EWH
%         mssaSort_Gap=struct();
%         mssaSort_Gap.name=Default_Method(me);
%         mssaSort_Gap.use_months=mssaSort(intit_num+1).use_months;
%         mssaSort_Gap.missing_months=mssaSort(intit_num+1).missing_months;
%         mssaSort_Gap.data_year_beg=mssaSort(intit_num+1).data_year_beg;
% 
%         % expand to one struct per institution
%         mssaSort_Gap = repmat(mssaSort_Gap, intit_num, 1);
% 
%         for ins = 1:intit_num
%             mssaSort_Gap(ins).name = mssaSort(ins).name;
%         end

        [lat_m, lon_n] = size(latlat);

        time_num = numel(fill_dates);
        grid_num = lat_m * lon_n;

        EST_MSSA = nan(time_num, intit_num, lat_m, lon_n);
        EST_MSSA_iter_num = nan(intit_num, lat_m, lon_n);

        % time x institution x grid
        TS_all = zeros(time_num, intit_num, grid_num);
        valid_grid = true(grid_num, 1);

        for ins = 1:intit_num
            A = mssaSort(ins).fill_Grid_EWH_SHC;          % time x lat x lon
            A2 = reshape(A, time_num, grid_num);          % time x grid

            TS_all(:,ins,:) = reshape(A2, time_num, 1, grid_num);
            valid_grid = valid_grid & all(~isnan(A2), 1).';
        end

        opts = struct();
        opts.maxIter = 100;
        opts.chiThreshold = 0.1;
        opts.doTest = false;
        opts.verbose = false;
        opts.useSvds = true;
        opts.fig = false;        % 批量跑时保持 false

        EST_MSSA_cell = cell(grid_num, 1);
        EST_iter_cell = cell(grid_num, 1);

        % 如果没有 Parallel Computing Toolbox，把 parfor 改成 for
        parfor g = 1:grid_num

            if ~valid_grid(g)
                continue
            end

            EWHA1_CJGI = TS_all(:,:,g);  % time x institution

            [MSSAn_CGJI_fill_reconst, ~, iter_num] = ...
                MSSA_gap_fitting_only_correct_fast(mssaSort_Gap, EWHA1_CJGI, M, N, opts);

            replace_CJGI = MSSAn_CGJI_fill_reconst{max(iter_num) + 1};

            EWHAf_CJGI = EWHA1_CJGI;

            for ins = 1:intit_num
                missing_months = mssaSort_Gap(ins).missing_months;
                EWHAf_CJGI(missing_months, ins) = replace_CJGI(missing_months, ins);
            end

            EST_MSSA_cell{g} = EWHAf_CJGI;
            EST_iter_cell{g} = iter_num;
        end

        % collect parfor outputs
        for g = 1:grid_num
            if isempty(EST_MSSA_cell{g})
                continue
            end

            [mm, nn] = ind2sub([lat_m, lon_n], g);

            EST_MSSA(:,:,mm,nn) = EST_MSSA_cell{g};
            EST_MSSA_iter_num(:,mm,nn) = EST_iter_cell{g}(:);
        end

        fprintf('EWH valid grids: %d / %d\n', sum(valid_grid), grid_num);
        fprintf('EWH skipped NaN grids: %d\n', grid_num - sum(valid_grid));

        % small bug here
        %         for i=1:intit_num
        %             EWH_EST_MSSA=squeeze(EST_MSSA(:,i,:,:));
        %             EWH_EST_MSSA_IterationN=squeeze(EST_MSSA_iter_num(i,:,:));
        %         end

        % Keep individual channel results for checking
        EWH_EST_MSSA_All = EST_MSSA;
        EWH_EST_MSSA_IterationN_All = EST_MSSA_iter_num;

        % Final product: average over real institutions/channels
        % EST_MSSA: time_num x intit_num x lat_m x lon_n
        EWH_EST_MSSA = squeeze(mean(EST_MSSA, 2, 'omitnan'));

        % EST_MSSA_iter_num: intit_num x lat_m x lon_n
        EWH_EST_MSSA_IterationN = squeeze(mean(EST_MSSA_iter_num, 1, 'omitnan'));

        fill_MASS_MSSA_SHC=zeros(1,numel(fill_dates));
        fill_EWH_MSSA_SHC=zeros(1,numel(fill_dates));
        for i=1:numel(fill_dates)
            %         TWSlmcosi=squeeze(fill_SHCdeltacoffs(i,:,:));
            %         [fill_Grid_SHC_mon,~,~,~,~]=plm2xyz(TWSlmcosi,1,c11cmn_combine,Lwindow);
            r_ewh_each=squeeze(EWH_EST_MSSA(i,:,:)); % cm
            fill_MASS_MSSA_SHC(i)=sum(r_ewh_each(in).*c11cmn_area(in))/sum(c11cmn_area(in))*10/10^3*BasinArea*10^3/10^3/10^9; % change from cm to Gt
            fill_EWH_MSSA_SHC(i)=sum(r_ewh_each(in).*c11cmn_area(in))/sum(c11cmn_area(in)); % cm

        end

        if saveAddData
            save(matFile2_me, ...
                'EWH_EST_MSSA', ...
                'EWH_EST_MSSA_All', ...
                'EWH_EST_MSSA_IterationN', ...
                'EWH_EST_MSSA_IterationN_All', ...
                'fill_EWH_MSSA_SHC', ...
                'fill_MASS_MSSA_SHC')
        end

    end

    switch me-1
        case 0

            mssaSort_Me.EWH.fillM_Grid_None=EWH_EST_MSSA;
            mssaSort_Me.EWH.fillM_None=fill_EWH_MSSA_SHC;
            mssaSort_Me.MASS.fillM_None=fill_MASS_MSSA_SHC;
        case 1

            mssaSort_Me.EWH.fillM_Grid_Gau300=EWH_EST_MSSA;
            mssaSort_Me.EWH.fillM_Gau300=fill_EWH_MSSA_SHC;
            mssaSort_Me.MASS.fillM_Gau300=fill_MASS_MSSA_SHC;
        case 2

            mssaSort_Me.EWH.fillM_Grid_Gau500=EWH_EST_MSSA;
            mssaSort_Me.EWH.fillM_Gau500=fill_EWH_MSSA_SHC;
            mssaSort_Me.MASS.fillM_Gau500=fill_MASS_MSSA_SHC;
        case 3

            mssaSort_Me.EWH.fillM_Grid_DDK3=EWH_EST_MSSA;
            mssaSort_Me.EWH.fillM_DDK3=fill_EWH_MSSA_SHC;
            mssaSort_Me.MASS.fillM_DDK3=fill_MASS_MSSA_SHC;
        case 4

            mssaSort_Me.EWH.fillM_Grid_DDK4=EWH_EST_MSSA;
            mssaSort_Me.EWH.fillM_DDK4=fill_EWH_MSSA_SHC;
            mssaSort_Me.MASS.fillM_DDK4=fill_MASS_MSSA_SHC;
        case 5

            mssaSort_Me.EWH.fillM_Grid_DDK5=EWH_EST_MSSA;
            mssaSort_Me.EWH.fillM_DDK5=fill_EWH_MSSA_SHC;
            mssaSort_Me.MASS.fillM_DDK5=fill_MASS_MSSA_SHC;
        case 6

            mssaSort_Me.EWH.fillM_Grid_DDK6=EWH_EST_MSSA;
            mssaSort_Me.EWH.fillM_DDK6=fill_EWH_MSSA_SHC;
            mssaSort_Me.MASS.fillM_DDK6=fill_MASS_MSSA_SHC;
        case 7

            mssaSort_Me.EWH.fillM_Grid_DDK7=EWH_EST_MSSA;
            mssaSort_Me.EWH.fillM_DDK7=fill_EWH_MSSA_SHC;
            mssaSort_Me.MASS.fillM_DDK7=fill_MASS_MSSA_SHC;
    end

    if strcmp(mssaSort_Me.process{3},'ocean')

        mssaSort_Me.MSL.unit='cm';

        if exist(matFile2_me,'file')

            matObj = matfile(matFile2_me);

            fileInfo = whos(matObj);
            varNames = {fileInfo.name};

            varExists1 = ismember('MSL_EST_MSSA', varNames);
            varExists2 = ismember('fill_MSL_MSSA_SHC', varNames);

            flag=varExists1 & varExists2;

        else
            flag=0;
        end

        if ~redo && flag

            disp(['Loading: ' matFile2_me ])

            load(matFile2_me, 'MSL_EST_MSSA', ...
                'fill_MSL_MSSA_SHC');

        else

            [lat_m, lon_n] = size(latlat);

            time_num = numel(fill_dates);
            grid_num = lat_m * lon_n;

            EST_MSSA = nan(time_num, intit_num, lat_m, lon_n);
            EST_MSSA_iter_num = nan(intit_num, lat_m, lon_n);

            % time x institution x grid
            TS_all = zeros(time_num, intit_num, grid_num);
            valid_grid = true(grid_num, 1);

            for ins = 1:intit_num
                A = mssaSort(ins).fill_Grid_MSL_SHC;          % time x lat x lon
                A2 = reshape(A, time_num, grid_num);          % time x grid

                TS_all(:,ins,:) = reshape(A2, time_num, 1, grid_num);
                valid_grid = valid_grid & all(~isnan(A2), 1).';
            end

            opts = struct();
            opts.maxIter = 100;
            opts.chiThreshold = 0.1;
            opts.doTest = false;
            opts.verbose = false;
            opts.useSvds = true;
            opts.fig = false;

            EST_MSSA_cell = cell(grid_num, 1);
            EST_iter_cell = cell(grid_num, 1);

            %
            parfor g = 1:grid_num

                if ~valid_grid(g)
                    continue
                end

                MSLA1_CJGI = TS_all(:,:,g);  % time x institution

                [MSSAn_CGJI_fill_reconst, ~, iter_num] = ...
                    MSSA_gap_fitting_only_correct_fast(mssaSort_Gap, MSLA1_CJGI, M, N, opts);

                replace_CJGI = MSSAn_CGJI_fill_reconst{max(iter_num) + 1};

                MSLAf_CJGI = MSLA1_CJGI;

                for ins = 1:intit_num
                    missing_months = mssaSort_Gap(ins).missing_months;
                    MSLAf_CJGI(missing_months, ins) = replace_CJGI(missing_months, ins);
                end

                EST_MSSA_cell{g} = MSLAf_CJGI;
                EST_iter_cell{g} = iter_num;
            end

            % collect parfor results
            for g = 1:grid_num
                if isempty(EST_MSSA_cell{g})
                    continue
                end

                [mm, nn] = ind2sub([lat_m, lon_n], g);

                EST_MSSA(:,:,mm,nn) = EST_MSSA_cell{g};
                EST_MSSA_iter_num(:,mm,nn) = EST_iter_cell{g}(:);
            end

            % optional summary
            fprintf('Valid grids: %d / %d\n', sum(valid_grid), grid_num);
            fprintf('Skipped NaN grids: %d\n', grid_num - sum(valid_grid));


            %             for i=1:intit_num
            %                 MSL_EST_MSSA=squeeze(EST_MSSA(:,i,:,:));
            %                 MSL_EST_MSSA_IterationN=squeeze(EST_MSSA_iter_num(i,:,:));
            %             end

            % Keep individual channel results for checking
            MSL_EST_MSSA_All = EST_MSSA;
            MSL_EST_MSSA_IterationN_All = EST_MSSA_iter_num;

            % Final product: average over real institutions/channels
            % EST_MSSA: time_num x intit_num x lat_m x lon_n
            MSL_EST_MSSA = squeeze(mean(EST_MSSA, 2, 'omitnan'));

            % EST_MSSA_iter_num: intit_num x lat_m x lon_n
            MSL_EST_MSSA_IterationN = squeeze(mean(EST_MSSA_iter_num, 1, 'omitnan'));

            for i=1:numel(fill_dates)
                %         TWSlmcosi=squeeze(fill_SHCdeltacoffs(i,:,:));
                %         [fill_Grid_SHC_mon,~,~,~,~]=plm2xyz(TWSlmcosi,1,c11cmn_combine,Lwindow);
                r_ewh_each=squeeze(MSL_EST_MSSA(i,:,:)); % cm
                %             fill_MASS_SHC(i)=sum(r_ewh_each(in).*c11cmn_area(in))/sum(c11cmn_area(in))*10/10^3*BasinArea*10^3/10^3/10^9; % change from cm to Gt
                fill_MSL_MSSA_SHC(i)=sum(r_ewh_each(in).*c11cmn_area(in))/sum(c11cmn_area(in)); % cm

            end

            if saveAddData
                save(matFile2_me, ...
                    'MSL_EST_MSSA', ...
                    'MSL_EST_MSSA_All', ...
                    'MSL_EST_MSSA_IterationN', ...
                    'MSL_EST_MSSA_IterationN_All', ...
                    'fill_MSL_MSSA_SHC', ...
                    '-append')
            end

        end

        switch me-1

            case 0

                mssaSort_Me.MSL.fillM_Grid_None=MSL_EST_MSSA;
                mssaSort_Me.MSL.fillM_None=fill_MSL_MSSA_SHC;

            case 1

                mssaSort_Me.MSL.fillM_Grid_Gau300=MSL_EST_MSSA;
                mssaSort_Me.MSL.fillM_Gau300=fill_MSL_MSSA_SHC;

            case 2

                mssaSort_Me.MSL.fillM_Grid_Gau500=MSL_EST_MSSA;
                mssaSort_Me.MSL.fillM_Gau500=fill_MSL_MSSA_SHC;

            case 3

                mssaSort_Me.MSL.fillM_Grid_DDK3=MSL_EST_MSSA;
                mssaSort_Me.MSL.fillM_DDK3=fill_MSL_MSSA_SHC;

            case 4

                mssaSort_Me.MSL.fillM_Grid_DDK4=MSL_EST_MSSA;
                mssaSort_Me.MSL.fillM_DDK4=fill_MSL_MSSA_SHC;

            case 5

                mssaSort_Me.MSL.fillM_Grid_DDK5=MSL_EST_MSSA;
                mssaSort_Me.MSL.fillM_DDK5=fill_MSL_MSSA_SHC;

            case 6

                mssaSort_Me.MSL.fillM_Grid_DDK6=MSL_EST_MSSA;
                mssaSort_Me.MSL.fillM_DDK6=fill_MSL_MSSA_SHC;

            case 7

                mssaSort_Me.MSL.fillM_Grid_DDK7=MSL_EST_MSSA;
                mssaSort_Me.MSL.fillM_DDK7=fill_MSL_MSSA_SHC;

        end

    end

end

save(matFile2, 'mssaSort_Me')

end

function tf = cache_has_smooth_fields(mssaSort_ins, method_order)

tf=true;
for me=method_order
    shcField=smooth_method_shc_field(me);
    if ~isfield(mssaSort_ins,shcField)
        tf=false;
        return
    end
end

end

function shcField = smooth_method_shc_field(me)

switch me-1
    case 0
        shcField = 'fill_SHCfilldelta';
    case 1
        shcField = 'fill_SHCfilldelta_Gau300';
    case 2
        shcField = 'fill_SHCfilldelta_Gau500';
    case 3
        shcField = 'fill_SHCfilldelta_DDK3';
    case 4
        shcField = 'fill_SHCfilldelta_DDK4';
    case 5
        shcField = 'fill_SHCfilldelta_DDK5';
    case 6
        shcField = 'fill_SHCfilldelta_DDK6';
    case 7
        shcField = 'fill_SHCfilldelta_DDK7';
    otherwise
        error('Unsupported smoothing method index: %d', me)
end

end
