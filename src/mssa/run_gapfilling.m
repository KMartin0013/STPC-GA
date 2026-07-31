function [slepcoff_gapfilling, iternum_gapfilling] = run_gapfilling(mssaSort, S, S_ol, ...
    mssaInfo, Attach, doPlot);

intit_num   = mssaInfo.intit_num;

fill_dates  = mssaSort(1).fill_dates;
fill_nomths = numel(fill_dates);
ddr1        = Attach.fig_path_ALL;
% doPlot      = true;
%% Should we fill the gap?

%%%%%% The determination of (M, N) for gap-filling procedure is by other simulation

if any(mssaSort(1).missing_months)

    % consistent M
    M_gap=mssaInfo.M_gap; % 96 months

    if M_gap > fill_nomths / 2
        error('M should not be larger than half of your data length')
    end

    % cutoff number N
    N_gap=mssaInfo.N_gap;

    disp(['The gap-filling procedure parameter: M = ' num2str(M_gap) ', N = ' num2str(N_gap)])

else

    disp('No gaps in the data.');
end

%% Gap filling (M-SSA) following (Gauer et al., 2022) and calculate the SSF

%
slepcoff_gapfilling=zeros(fill_nomths,intit_num,S_ol);
iternum_gapfilling=zeros(intit_num,S_ol);
for ss=1:S_ol

    ss

    NaNflag=0;

    mssaSort_gap=mssaSort(1:intit_num);

    % substitute the MSSA_TS as the fill_deltacoffs for
    % each grid
    for ins=1:intit_num

        if ss<S+1
            XXXX_mn=mssaSort_gap(ins).fill_deltacoffs;
            mssaSort_gap(ins).MSSA_TS=XXXX_mn(:,ss)';
            mssaSort_gap(ins).MSSA_TS_order=ss;
        else
            disp('we need to consider IB for ocean study.')
            XXXX_mn=mssaSort_gap(ins).fill_IBdeltacoffs;
            mssaSort_gap(ins).MSSA_TS=XXXX_mn';
            mssaSort_gap(ins).MSSA_TS_order='IB';
        end

        if sum(isnan(mssaSort_gap(ins).MSSA_TS))
            NaNflag=1;
            break
        end

    end

    % don't do MSSA for NaN
    if NaNflag
        warning(['NaN exists for m: ' num2str(ss)])
        continue
    end

    % if no missing months exist, just expand the coefficient
    MSLAf_Centers=zeros(numel(fill_dates),intit_num);
    for ins=1:intit_num
        MSLAf_Centers(:,ins)=mssaSort_gap(ins).MSSA_TS';
    end

    % if missing months exist, replace the gap filling
    if any(mssaSort(1).missing_months)

        % new version 1.2, correct for calculating std only for missing values
        % 3 sigma
        if doPlot

%             if mod(ss-1,5)==0
                
                [MSSAn_CGJI_fill_reconst,~,iter_num,chi,~]=MSSA_gap_fitting_correct(mssaSort_gap,M_gap,N_gap,...
                    fullfile(ddr1,['MSSA_N' num2str(ss) '_Final']));
%             else

%                 [MSSAn_CGJI_fill_reconst,~,iter_num,chi,~]=MSSA_gap_fitting_correct(mssaSort_gap,M_gap,N_gap);

%             end

        else

            [MSSAn_CGJI_fill_reconst,~,iter_num,chi,~]=MSSA_gap_fitting_correct(mssaSort_gap,M_gap,N_gap);
        end
        replace_Centers=MSSAn_CGJI_fill_reconst{max(iter_num)+1};

        for ins=1:intit_num

            MSLAf_Centers(mssaSort(ins).missing_months,ins)=replace_Centers(mssaSort(ins).missing_months,ins);

        end

    else

        iter_num=zeros(intit_num,1);

    end

    slepcoff_gapfilling(:,:,ss)=MSLAf_Centers;
    iternum_gapfilling(:,ss)=iter_num;

end

end
