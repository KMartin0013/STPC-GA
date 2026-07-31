function slepianInfo = prepare_slepian_info(config)
% PREPARE_SLEPIAN_INFO  Build struct to be saved as Slepian_Information.mat
%
% Usage:
%   slepianInfo = prepare_slepian_info(config)

    slepianInfo = struct();
    slepianInfo.buffer_deg = config.buffer_deg;
    if strcmp(config.buffer_deg, 'Auto')
        if isfield(config, 'groupBuffer')
            if numel(config.groupBuffer) < 2
                error('V3:InvalidGroupBuffer', ...
                    'Provide at least two candidate buffer zones.');
            elseif ~(all(config.groupBuffer >= 0) || all(config.groupBuffer <= 0))
                error('V3:InvalidGroupBufferSigns', ...
                    'Candidate buffer zones must use consistent signs.');
            end
            slepianInfo.group_buffer = config.groupBuffer;
        elseif strcmp(config.landOrOcean, 'ocean')
            slepianInfo.group_buffer = 0:-0.5:-1.5;
        else
            slepianInfo.group_buffer = 0:0.5:1.5;
        end
    end

    slepianInfo.S_selection        = config.v3Selection.S;
    slepianInfo.Radius             = config.Radius;
    slepianInfo.phi                = config.phi;
    slepianInfo.theta              = config.theta;
    slepianInfo.omega              = config.omega;
    slepianInfo.artificial_months  = config.artificialMonths;
end
