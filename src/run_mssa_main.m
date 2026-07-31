function mssa_result = run_mssa_main(basicInfo, slepianInfo, mssaInfo, ...
    slepian_results)
% RUN_MSSA_MAIN  Do MSSA-based gap filling and reconstruction on Slepian results
%
% Usage:
%   mssaRes = run_mssa_main(basicInfo, slepianInfo, mssaInfo, instResults)
%
% Inputs:
%   basicInfo   - struct, same definition as in run_slepian_main:
%       .Code_Version
%       .Area
%       .Lwindow
%       .c11cmn
%       .ddir1, .ddir2
%       .ifilesRoot
%
%   slepianInfo - struct, same as in run_slepian_main:
%       .Dataproduct
%       .buffer_str
%       .Max_S
%       .Radius
%       .use_institu     % (recommended) list of institution names
%
%   mssaInfo    - struct, corresponding to original MSSA_Information.mat:
%       .S_bou          % lower/upper lambda thresholds for S
%       .N_bou          % lower/upper cumulative EV thresholds for N
%       .Turning_number % number of change points to detect
%       .M_gap_default  % window size M for MSSA gap filling
%       .N_gap_default  % number of RCs N for MSSA gap filling
%       .redo           % optional flag to force recomputation (not used here)
%
%   instResults - struct array, length = intit_num or intit_num+1
%       Each element is the output of run_slepian_main for one institution:
%       .Dataproduct
%       .TH, .XY, .BasinArea
%       .S, .S_sig, .S_shannon, .CC, .V, .Radius
%       .time, .fill, .grid, .date, .ocean
%
% Output:
%   mssaRes - struct, containing:
%       .A5_Centers_Sort, .A5_Centers_Fit
%       .common_use_months
%       .turn_V_four, .turn_MSSA_four
%       .M_rec, .N_gap, .M_gap
%       .MSSAf_evalues_sumup, .RC_RMS_Seq, .V, .Max_S, .S_ol, .BasinArea
%       .lon, .lat, .r_record
%
% References: see run_slepian_main.

arguments
    basicInfo       struct
    slepianInfo     struct
    mssaInfo        struct
    slepian_results     struct
end

%% 0. Unpack basic information
Code_Version = basicInfo.Code_Version;
Area         = basicInfo.Area;
Lwindow      = basicInfo.Lwindow;
c11cmn       = basicInfo.c11cmn;
ddir1        = basicInfo.ddir1;      % output results directory (Results_*)
ddir2        = basicInfo.ddir2;      % output figure directory (Figures_*)
ifilesRoot   = basicInfo.ifilesRoot; % root directory of original IFILES
redo         = basicInfo.redo;
plotProcess  = basicInfo.plotProcess;
saveAddData  = basicInfo.saveAddData;
% Dataproduct  = slepianInfo.Dataproduct;
% buffer_str   = slepianInfo.buffer_str;
% Max_S        = slepianInfo.Max_S;
% Radius       = slepianInfo.Radius;

buffer_str   = slepianInfo.buffer_str;

use_institu  = mssaInfo.use_institu;
intit_num    = mssaInfo.intit_num;
M_gap        = mssaInfo.M_gap;
N_gap        = mssaInfo.N_gap;
M            = mssaInfo.M;

%% 1. Build paths and filename tags (similar to FIG_Attach in run_slepian_main)
[Attach] = build_mssa_paths_and_tags(basicInfo, buffer_str, ...
    mssaInfo);

% Ensure all directories exist
ensure_dirs_exist(Attach);

%% 2. Collect institution results into A5_Centers_Sort and compute mean product

disp(['==== Averaging multiple data centers (e.g., ' ...
    char(use_institu(1)) ', ' char(use_institu(2)) ') to ' ...
    char(use_institu(end)) '. ===='])

% redo=true;
if redo || ~exist(fullfile(ddir1,['MainMSSA_' Attach.Attach_ALL '.mat']),'file')

    [mssa_Sort, S, S_ol, SHC_degord] = collect_and_average_institutions(basicInfo, ...
        slepian_results, Attach);

    %% 3. MSSA gap-filling (MSSA_gap_filling_correct)
    [fillM_deltacoffs, fillM_iternum] = run_gapfilling(mssa_Sort, S, S_ol, ...
        mssaInfo, Attach, plotProcess);

    %% 4. Search for optimal M for reconstruction (separation index)
    if (ischar(M) || (isstring(M) && isscalar(M))) && ...
            strcmpi(char(string(M)), 'Auto')

        M_rec = det_opt_M(mssa_Sort, mssaInfo, fillM_deltacoffs, ...
            S_ol, S, Attach, plotProcess);

    elseif isnumeric(M) && isscalar(M) && isfinite(M) && ...
            M > 0 && M == fix(M)

        if M > size(fillM_deltacoffs, 1)
            error('V3:FixedMOutOfRange', ...
                ['Requested M=%d exceeds the filled time-series length ', ...
                 'of %d months.'], M, size(fillM_deltacoffs, 1));
        end
        M_rec = double(M);

    else

        error('V3:InvalidM', ...
            'M must be "STPC" (normalized to "Auto") or a positive integer.')
    end

    %%

    mssa_result.use_institu    = mssaInfo.use_institu;
    mssa_result.intit_num  = mssaInfo.intit_num;
    mssa_result.M_gap      = mssaInfo.M_gap;
    mssa_result.N_gap      = mssaInfo.N_gap;

    mssa_result.Attach     = Attach;
    mssa_result.S          = S;
    mssa_result.S_ol       = S_ol;
    mssa_result.SHC_degord = SHC_degord;
    mssa_result.M_rec      = M_rec;
    mssa_result.fillM_deltacoffs    = fillM_deltacoffs;
    mssa_result.fillM_iternum       = fillM_iternum;
    mssa_result.mssa_Sort  = mssa_Sort;
    
    matFile = fullfile(ddir1, ['MainMSSA_' Attach.Attach_ALL '.mat']);
    save(matFile, '-struct', 'mssa_result');

else
    
    % directly load this variables
    disp(['MSSA: Loading ',fullfile(ddir1, ['MainMSSA_' Attach.Attach_ALL '.mat'])]);
    mssa_result = load(fullfile(ddir1,['MainMSSA_' Attach.Attach_ALL '.mat']));

end

disp(['==== End for preparing ' char(use_institu(end)) '. ===='])

end
