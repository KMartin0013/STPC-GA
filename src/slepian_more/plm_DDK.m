function lmcosi_DDK=plm_DDK(lmcosi,number)
% plm_Gauss=plm_Gaussian_smoothing(lmcosi,order,radius)
%
% Converts the real spherical harmonic coefficients expressing a geoPOTENTIAL
% field into a smoothed version with DDK filter of the smoothing number N
%
% INPUT:
%
% lmcosi     [l m Ccos Csin] degrees, order, coefficients, as in PLM2XYZ
%            Units are properly those of potential, J/kg or m^2/s^2
% radius     The averaging radius, in km [default: 500]
% a          Reference radius [default: a_EGM96]
%
% OUPUT:
%
% lmcosig    [l m Ccos Csin] degrees, order, smoothed coefficients of 
%            output in the standard units of what's requested [J/kg, m/s^2, m]
%
% DDK1d11: filtered with inverse signal degree power law 1e11*deg^4 (DDK5) weakest smoothing
% DDK5d11: filtered with inverse signal degree power law 5e11*deg^4 (DDK4)         |
% DDK1d12: filtered with inverse signal degree power law 1e12*deg^4 (DDK3)         |
% DDK1d13: filtered with inverse signal degree power law 1e13*deg^4 (DDK2)         |
% DDK1d14: filtered with inverse signal degree power law 1e14*deg^4 (DDK1) strongest smoothing
%
% Last modified by zhongtian.ma@connect.polyu.hjk, 15/07/2025

% This function is originally created based on the code from Roelof Rietbroek.
% Copyright Roelof Rietbroek 2016
% This software is licensed under the MIT license see https://github.com/strawpants/GRACE-filter/blob/master/LICENSE
% URL: https://github.com/strawpants/GRACE-filter

% This function also partly references to GRAMAT toolbox by Feng Wei.
% https://github.com/fengweiigg/GRACE_Matlab_Toolbox

% get the DDK filtering document
defval('number',3)

dependencyRoot = stpc_dependency_root();
defval('ddir_Code',fullfile(dependencyRoot,...
    'GRACE-filter-master','src','matlab'));
defval('ddir_DDK',fullfile(dependencyRoot,...
    'GRACE-filter-master','data','DDK'));

if exist(ddir_Code, 'dir') ~= 7 || exist(ddir_DDK, 'dir') ~= 7
    error('STPCGA:MissingDDKDependency', ...
        ['GRACE-filter-master is incomplete under "%s". ', ...
         'Install its src/matlab and data/DDK directories.'], ...
        dependencyRoot);
end
addpath(ddir_Code)
switch number
    case 1 %strongest DDK1
        file='Wbd_2-120.a_1d14p_4';
    case 2 % DDK2
        file='Wbd_2-120.a_1d13p_4';
    case 3 % DDK3
        file='Wbd_2-120.a_1d12p_4';
    case 4 % DDK4
        file='Wbd_2-120.a_5d11p_4';
    case 5 % DDK5
        file='Wbd_2-120.a_1d11p_4';
	case 6
		file='Wbd_2-120.a_5d10p_4';
	case 7
		file='Wbd_2-120.a_1d10p_4';
	case 8
		file='Wbd_2-120.a_5d9p_4';		
end

dat=read_BIN(fullfile(ddir_DDK,[file]));

degree_max=lmcosi(end,1);
cnm=zeros(degree_max);
snm=zeros(degree_max);

for i=1:size(lmcosi,1)
    cnm(lmcosi(i,1)+1,lmcosi(i,2)+1)=lmcosi(i,3);
    snm(lmcosi(i,1)+1,lmcosi(i,2)+1)=lmcosi(i,4);
end

% DDK filter
[cnmfilt,snmfilt]=filterSH_m(dat,cnm,snm);

lmcosi_DDK=zeros(size(lmcosi));
lmcosi_DDK(:,1:2)=lmcosi(:,1:2);
for i=1:size(lmcosi,1)

    lmcosi_DDK(i,3)=cnmfilt(lmcosi_DDK(i,1)+1,lmcosi_DDK(i,2)+1);
    lmcosi_DDK(i,4)=snmfilt(lmcosi_DDK(i,1)+1,lmcosi_DDK(i,2)+1);

end

end


