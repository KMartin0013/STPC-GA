function basicInfo = prepare_basic_info(codeVersion, config)
% PREPARE_BASIC_INFO  Build struct to be saved as Basic_Information.mat
%
% Usage:
%   basicInfo = prepare_basic_info(codeVersion, config)

    config = stpc_configure_center_list(config);
    if ~isfield(config, 'note')
        config.note = '';
    end
    config.note = stpc_note_suffix(config.note);

    basicInfo = struct();
    basicInfo.Code_Version  = codeVersion;
    basicInfo.TH_ori        = config.TH_ori;
    basicInfo.Area          = config.areaName;
    basicInfo.Lwindow       = config.Lwindow;
    basicInfo.land_or_ocean = config.landOrOcean;
    basicInfo.c11cmn        = config.c11cmn;
    basicInfo.Max_S         = config.Max_S;
    basicInfo.ddir1         = config.resultDir;
    basicInfo.ddir2         = config.figureDir;
    basicInfo.ddir3         = config.AddresultDir;
    basicInfo.resultRoot    = config.resultRoot;
    basicInfo.sharedResultDir = config.sharedResultDir;
    basicInfo.figureRoot    = config.figureRoot;
    basicInfo.caseName      = config.caseName;
    basicInfo.ifilesRoot    = config.ifilesRoot;
    basicInfo.redo          = config.redo;
    basicInfo.plotProcess   = config.plotProcess;
    basicInfo.saveAddData   = config.saveAddData;
    basicInfo.note          = config.note;
    basicInfo.Smooth        = config.Smooth;

    basicInfo.Lwindow       = config.Lwindow;
    basicInfo.GIA           = config.GIA;
    basicInfo.Institu_ver   = config.Institu_ver;
    basicInfo.use_institu   = config.use_institu;
    basicInfo.center_list   = config.center_list;
    basicInfo.ensInstitu    = config.ensInstitu;
    basicInfo.v3Selection   = config.v3Selection;

end
