function configuration = v3_configuration_snapshot(config)
%V3_CONFIGURATION_SNAPSHOT Return science-affecting V3 configuration only.

configuration = struct();
configuration.schema_version = 'V3-configuration-1';
configuration.selection = config.v3Selection;

fields = { ...
    'codeVersionBase', 'TH_ori', 'areaName', 'Lwindow', ...
    'landOrOcean', 'c11cmn', ...
    'Max_S', 'groupBuffer', 'M_gap', 'N_gap', 'turningNumber', ...
    'S_bound', 'N_bound', 'p_use', 'GIA', 'Institu_ver', ...
    'givInstitu', 'Smooth', 'artificialMonths', ...
    'phi', 'theta', 'omega'};

for i = 1:numel(fields)
    name = fields{i};
    if isfield(config, name)
        configuration.(name) = config.(name);
    end
end
end
