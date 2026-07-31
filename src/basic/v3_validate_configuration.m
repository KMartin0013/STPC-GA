function configuration = v3_validate_configuration(config)
%V3_VALIDATE_CONFIGURATION Guard a configuration directory against cache mixing.

configuration = v3_configuration_snapshot(config);
configurationFile = fullfile(config.resultDir, 'Configuration.mat');

if exist(configurationFile, 'file')
    saved = load(configurationFile, 'configuration');
    if ~isfield(saved, 'configuration')
        error('V3:InvalidConfigurationFile', ...
            'Existing %s does not contain variable "configuration".', ...
            configurationFile);
    end
    if ~isequaln(saved.configuration, configuration)
        if config.redo
            warning('V3:ConfigurationChangedWithRedo', ...
                ['Configuration differs from the existing case "%s". ', ...
                 'redo=true allows recomputation and replaces Configuration.mat.'], ...
                config.caseName);
        else
            error('V3:ConfigurationMismatch', ...
                ['Configuration differs from existing case "%s". ', ...
                 'Choose a new config.caseName or set redo=true to recompute.'], ...
                config.caseName);
        end
    else
        return
    end
end

save(configurationFile, 'configuration');
end
