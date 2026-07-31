function center_stpc = run_stpc_center_specific_selected(basicInfo, STPCInfo, ...
    slepian_results, mssa_results, p_S, p_N, used_sta_SN, varargin)
%RUN_STPC_CENTER_SPECIFIC_SELECTED Export center-specific STPC results.
%
% This lightweight workflow reuses the selected S/N choices from the main
% STPC run and saves the center-specific outputs that are normally averaged
% into the final CJGI result. It does not repeat Slepian construction or
% search all S/N combinations.

if nargin < 7 || isempty(used_sta_SN)
    used_sta_SN = 4;
end
opts = local_parse_options(varargin{:});

ddir1 = basicInfo.ddir1;
if isfield(basicInfo, 'ddir3')
    ddir3 = basicInfo.ddir3;
else
    ddir3 = fullfile(ddir1, 'AddData');
end
if exist(ddir3, 'dir') ~= 7
    mkdir(ddir3);
end

Attach = mssa_results.Attach;
default_output = fullfile(ddir3, ['MainSTPC_' Attach.Attach_ALL ...
    '_center_specific_selected.mat']);
if isempty(opts.output_file)
    opts.output_file = default_output;
end
if ~opts.redo && exist(opts.output_file, 'file') == 2
    load(opts.output_file, 'center_stpc');
    disp(['STPC center-specific: Loading ' opts.output_file]);
    return
end

tp_file = fullfile(ddir1, ['MainSTPC_' Attach.Attach_ALL '_TP.mat']);
if exist(tp_file, 'file') ~= 2
    error('run_stpc_center_specific_selected:MissingTP', ...
        'Missing selected turning-point file: %s', tp_file);
end
load(tp_file, 'turn_V_four', 'turn_MSSA_four');

p_use = STPCInfo.p_use;
if isempty(opts.p_indices)
    opts.p_indices = 1:numel(p_use);
end

lon = slepian_results(1).grid.lon;
lat = slepian_results(1).grid.lat;
XY = slepian_results(1).XY;
BasinArea = slepian_results.BasinArea;
r_record = slepian_results(1).grid.r_record;

intit_num = mssa_results.intit_num;
mssa_Sort = mssa_results.mssa_Sort;
fillM_deltacoffs = mssa_results.fillM_deltacoffs;
M_rec = mssa_results.M_rec;
fill_dates = mssa_Sort(1).fill_dates;
fill_nmonths = numel(fill_dates);
center_names = string(mssa_results.use_institu(1:intit_num));

center_stpc = struct();
center_stpc.method = 'center-specific STPC selected-S/N export';
center_stpc.center_names = center_names;
center_stpc.n_centers = intit_num;
center_stpc.p_use = p_use(opts.p_indices);
center_stpc.p_indices = opts.p_indices;
center_stpc.p_S = p_S(opts.p_indices);
center_stpc.p_N = p_N(opts.p_indices);
center_stpc.used_sta_SN = used_sta_SN;
center_stpc.fill_dates = fill_dates;
center_stpc.lon = lon;
center_stpc.lat = lat;
center_stpc.EWH.unit = 'cm';

for out_i = 1:numel(opts.p_indices)
    p_i = opts.p_indices(out_i);
    Noise_SigLev = p_use(p_i);
    TP_ss = p_S(p_i);
    TP_nn = p_N(p_i);

    S_rec = turn_V_four(TP_ss, used_sta_SN);
    N_rec = zeros(1, S_rec);
    for k = 1:S_rec
        N_rec(k) = turn_MSSA_four{k}(TP_nn, used_sta_SN);
    end
    if strcmp(mssa_Sort(1).process{3}, 'ocean')
        S_ol = S_rec + 1;
        N_rec(S_rec + 1) = turn_MSSA_four{S_rec + 1}(TP_nn, used_sta_SN);
    else
        S_ol = S_rec;
    end

    fprintf('STPC center-specific selected export: p=%g, S index=%d, N index=%d, S_rec=%d.\n', ...
        Noise_SigLev, TP_ss, TP_nn, S_rec);

    [MSSAcoffs_rec, MSSAcoffs_evalues, MSSAcoffs_RC, MSSAcoffs_RCTest] = ...
        run_MSSA_decompose(mssa_Sort, fillM_deltacoffs, M_rec, N_rec, ...
        Noise_SigLev, fill_nmonths, intit_num, S_ol, S_rec);

    mssa_Sort_coff_sn = coffs_reconstruction(MSSAcoffs_RCTest, ...
        MSSAcoffs_RC, fillM_deltacoffs, MSSAcoffs_rec, ...
        MSSAcoffs_evalues, S_ol, fill_nmonths, intit_num);

    [mssa_Sort_EWHorMSL_sn, ~, ~, ~] = grids_reconstruction_EWHorMSL( ...
        mssa_Sort_coff_sn, intit_num, fill_nmonths, S_rec, S_ol, ...
        r_record, lat, lon, XY, BasinArea);

    for center_i = 1:intit_num
        center_stpc.EWH.fillM_Slepian_center(:, center_i, out_i) = ...
            mssa_Sort_EWHorMSL_sn(center_i).fillM_EWH(:);
        center_stpc.EWH.fillM_MSSA_center(:, center_i, out_i) = ...
            mssa_Sort_EWHorMSL_sn(center_i).fillM_EWH_MSSA(:);
        center_stpc.EWH.fillM_STPC_center(:, center_i, out_i) = ...
            mssa_Sort_EWHorMSL_sn(center_i).fillM_EWH_MSSA_both(:);

        center_stpc.EWH.fillM_Grid_STPC_center(:,:,:,center_i,out_i) = ...
            mssa_Sort_EWHorMSL_sn(center_i).fillM_Grid_EWH_MSSA_both;
    end
end

save(opts.output_file, 'center_stpc', '-v7.3');
disp(['STPC center-specific: Saved ' opts.output_file]);
end

function opts = local_parse_options(varargin)
opts = struct();
opts.output_file = '';
opts.redo = false;
opts.p_indices = [];

if mod(numel(varargin), 2) ~= 0
    error('Options must be name/value pairs.');
end
for i = 1:2:numel(varargin)
    name = lower(string(varargin{i}));
    value = varargin{i + 1};
    switch name
        case "output_file"
            opts.output_file = char(value);
        case "redo"
            opts.redo = logical(value);
        case "p_indices"
            opts.p_indices = value;
        otherwise
            error('Unknown option: %s', name);
    end
end
end
