function mssa_Sort_RMS_sn = grids_reconstruction_EWHorMSL_RMS( ...
    mssa_Sort_coff_sn,intit_num,fill_nmonths,S_rec,...
    r_record,lonlon)

% -------------------------------------------------------------------------
% compute_STD_RMS
%
% INPUT:
%   mssa_Sort_coff_sn   : struct array with MSSA coefficients:
%                       fields must include fillM_deltacoffs, fillM_reconcoffs,
%                       fillM_reconcoffs_both, etc.
%   intit_num      : number of MSSA iterations (scalar), loop 1:intit_num+1
%   fill_nmonths   : number of months in the time series
%   S_rec          : number of Slepian functions used
%   r_record       : precomputed Slepian basis on grid, size [S_rec, nlat, nlon]
%   lonlon, latlat : 2-D longitude/latitude meshgrid used to define grid size
%
% OUTPUT:
%   mssa_Sort_RMS_sn : struct array containing RMS fields:
%                      fillM_Grid_EWH_RMS,
%                      fillM_Grid_EWH_MSSA_RMS,
%                      fillM_Grid_EWH_MSSA_both_RMS,
%                      fillM_Grid_EWH_MSSA_bothfail_RMS,
%                      fillM_Grid_EWH_MSSA_resid_RMS, etc.
%
% NOTE:
%   - Direct translation of the original STD and RMS block (active lines only).
%   - All parameter names and computational steps are preserved.
% -------------------------------------------------------------------------

% calculate the STD for the reject and residuals (maybe?)
% notice that here we do not consider IB, and only consider the EWH (e.g., use 1000 kg/m^3)

mssa_Sort_RMS_sn = struct();

for ins = 1:intit_num+1

    % the original coefficients without MSSA has already been calculated
    % for the spatial and temporal results.
    fillM_deltacoffs_ins      = mssa_Sort_coff_sn(ins).fillM_deltacoffs;
%   fill_signaldelta_ins      = mssa_Sort_coff_sn(ins).fill_signaldelta;

    % (spatial)
    % the gap-filling part with MSSA
    fillM_reconcoffs_ins      = mssa_Sort_coff_sn(ins).fillM_reconcoffs;
    fillM_reconcoffs_both_ins = mssa_Sort_coff_sn(ins).fillM_reconcoffs_both;

    % These fields are intentionally per-mode (their names do not contain
    % "uptoS"), so the work arrays must be reset for every j.
    for j = 1:S_rec
        sp_ewh          = zeros(fill_nmonths,size(lonlon,1),size(lonlon,2));
        sp_ewh_mssa     = zeros(fill_nmonths,size(lonlon,1),size(lonlon,2));
        sp_ewh_mssa_both= zeros(fill_nmonths,size(lonlon,1),size(lonlon,2));

        r = squeeze(r_record(j,:,:));

        for i = 1:fill_nmonths
            % original EWH (cm)
            sp_ewh(i,:,:) = r * fillM_deltacoffs_ins(i,j)' / 1000 * 100;

            % MSSA reconstruction (cm)
            sp_ewh_mssa(i,:,:)      = r * fillM_reconcoffs_ins(i,j)'      / 1000 * 100;
            sp_ewh_mssa_both(i,:,:) = r * fillM_reconcoffs_both_ins(i,j)' / 1000 * 100;
        end

        % compute RMS at each grid cell
        for p = 1:size(lonlon,1)
            for q = 1:size(lonlon,2)

                % RMS of original EWH
                mssa_Sort_RMS_sn(ins).fillM_Grid_EWH_RMS(j,p,q) = ...
                    rms(squeeze(sp_ewh(:,p,q)));
%               
                % RMS of MSSA reconstruction
                mssa_Sort_RMS_sn(ins).fillM_Grid_EWH_MSSA_RMS(j,p,q) = ...
                    rms(squeeze(sp_ewh_mssa(:,p,q)));
                mssa_Sort_RMS_sn(ins).fillM_Grid_EWH_MSSA_both_RMS(j,p,q) = ...
                    rms(squeeze(sp_ewh_mssa_both(:,p,q)));

                % RMS of MSSA residuals and both-fail terms
                mssa_Sort_RMS_sn(ins).fillM_Grid_EWH_MSSA_bothfail_RMS(j,p,q) = ...
                    rms(squeeze(sp_ewh_mssa(:,p,q) - sp_ewh_mssa_both(:,p,q)));
                mssa_Sort_RMS_sn(ins).fillM_Grid_EWH_MSSA_resid_RMS(j,p,q) = ...
                    rms(squeeze(sp_ewh(:,p,q) - sp_ewh_mssa(:,p,q)));
            end
        end
    end
end

end
