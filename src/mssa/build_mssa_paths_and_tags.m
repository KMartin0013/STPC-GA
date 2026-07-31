function [Attach] = build_mssa_paths_and_tags(basicInfo, buffer_str, ...
    mssaInfo)

Code_Version = basicInfo.Code_Version;
Area         = basicInfo.Area;
Lwindow      = basicInfo.Lwindow;
use_institu  = basicInfo.use_institu;
ifilesRoot   = basicInfo.ifilesRoot;
Institu_ver  = basicInfo.Institu_ver;
resultDir    = basicInfo.ddir1;
figureDir    = basicInfo.ddir2;

intit_num    = mssaInfo.intit_num;

% S_bou        = STPCInfo.S_bou;
% N_bou        = STPCInfo.N_bou;

% suffix_area: e.g. 'RL_XX'
suffix_area = [Institu_ver(3:4) '_' Area];

suffix_XY_use = [num2str(Lwindow) '_' buffer_str];
note_suffix   = '';
if isfield(basicInfo, 'note')
    note_suffix = stpc_note_suffix(basicInfo.note);
end
suffix_str    = [suffix_area '_' suffix_XY_use note_suffix];

% if S_bou==0 && N_bou==0
%     suffix_bou = '';
% elseif S_bou==0
%     suffix_bou = ['_Hn' num2str(N_bou)];
% elseif N_bou==0
%     suffix_bou = ['_Hs' num2str(S_bou)];
% else
%     suffix_bou = ['_Hs' num2str(S_bou) 'n' num2str(N_bou)];
% end

Attach_ALL = [char(use_institu(intit_num+1)) '_' suffix_str];

% If you still follow the original IFILES/Text_* and Figure_* structure:
txt_path_ALL = fullfile(ifilesRoot, ['Text_'   Code_Version], Attach_ALL);
% fig_path_ALL = fullfile(ifilesRoot, ['Figure_' Code_Version], Attach_ALL);
fig_path_ALL = fullfile(figureDir, Attach_ALL);

Attach_each   = strings(intit_num+1,1);
txt_path_each = strings(intit_num+1,1);
fig_path_each = strings(intit_num+1,1);

for ins = 1:intit_num
    Attach_each(ins)   = string([char(use_institu(ins)) '_' suffix_str]);
    txt_path_each(ins) = string(fullfile(ifilesRoot, ...
        ['Text_'   Code_Version], [char(Attach_each(ins))]));
%     fig_path_each(ins) = string(fullfile(ifilesRoot, ...
%         ['Figure_' Code_Version], [char(Attach_each(ins))]));
    fig_path_each(ins) = string(fullfile(figureDir, ...
        [char(Attach_each(ins))]));

end

Attach_each(intit_num+1)   = string(Attach_ALL);
txt_path_each(intit_num+1) = string(txt_path_ALL);
fig_path_each(intit_num+1) = string(fig_path_ALL);

Attach=struct();
Attach.Attach_ALL   = Attach_ALL;
Attach.txt_path_ALL = txt_path_ALL;
Attach.fig_path_ALL = fig_path_ALL;
Attach.Attach_each  = Attach_each;
Attach.txt_path_each= txt_path_each;
Attach.fig_path_each= fig_path_each;

Attach.result_path_ALL  = resultDir;

end
