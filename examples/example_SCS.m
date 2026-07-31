% EXAMPLE_SCS  STPC-GA (Version 1) South China Sea example.
%
% Each caseName is an independent cache boundary under resultRoot and
% figureRoot.
%
% Shared simulation inputs remain under <resultRoot>/Simulate and are
% reused only after their embedded metadata passes validation.

clear; clc;
tic

config = struct();

%% Project, dependencies, and internal result tag
exampleFile = mfilename('fullpath');
config.ifilesRoot = fileparts(fileparts(exampleFile));
addpath(config.ifilesRoot);
startup_STPC_GA(config.ifilesRoot);
% This tag is written to output metadata. The public software release is
% STPC-GA (Version 1).
config.codeVersionBase = 'STPCGA1_';

%% Region
config.TH_ori        = 'SCSpTH';
config.areaName      = config.TH_ori;
config.landOrOcean   = 'ocean';
config.c11cmn        = [95.5 29.5 124.5 -4.5];

%% GRACE / GRACE-FO
config.Lwindow       = 60;
config.GIA           = 'Pelt17';
config.Institu_ver   = 'RL06';
config.givInstitu    = ...
    ["CSR0324","JPL0324","GFZ0324","ITSG0324"];

%% Public R/S/M/N interface
% R: "STPC" searches groupBuffer; a number fixes the buffer in degrees.
config.R = "STPC";
% config.R = 1;

% S: choose exactly one form.
config.S = "STPC";                         % default STPC search
% config.S = "Shannon";                    % Shannon number
% config.S = ["Fixed","0.1"];              % retain eigenvalues > 0.1
% config.S = ["Sel_Gau_default","500"];    % selective 500-km Gaussian
% config.S = 30;                           % directly retain 30 SSFs

% M: "STPC" searches the reconstruction window; a number fixes it.
config.M = "STPC";
% config.M = 120;

% N: "STPC" uses cumulative-eigenvalue turning points; "Wcorr" generates
% candidates from joint-center w-correlation and lets STPC choose among
% them; a number applies the same retained count to every coefficient.
config.N = "STPC";
% config.N = "Wcorr";
% config.N = 8;
%
% Optional advanced Wcorr settings (the defaults shown here are normally
% sufficient and are recorded in Configuration.mat):
% config.Wcorr = struct('blockWidth',6,'minSpacing',2,'maxModes',[]);

%% Search and method parameters
% Advanced setting. Default is 5. Some legacy process figures assume the
% default layout, so change it only when those plotting routines are valid.
config.turningNumber = 5;
config.groupBuffer   = [0,-0.5,-1,-1.5];
config.Max_S         = 50;
config.S_bound       = 0.05;
config.N_bound       = 0;
config.p_use         = [0.05, 0.1, 0.3];

% Selective Gaussian radius is supplied by config.S. Radius remains here
% for legacy helpers and is ignored by non-Gaussian S modes.
config.Radius        = 500;
config.phi           = 0;
config.theta         = 0;
config.omega         = 0;
config.artificialMonths = 5*12+1:6*12;

%% Gap filling (independent of the main reconstruction M/N)
config.M_gap         = 13;
config.N_gap         = 8;

%% Institution display name
ensInstitu = [];
for i = 1:numel(config.givInstitu)
    product = char(config.givInstitu(i));
    center = product(1:end-4);
    productTime = product(end-3:end);
    ensInstitu = [ensInstitu center(1)]; %#ok<AGROW>
end
ensInstitu = [ensInstitu productTime];
config.ensInstitu  = ensInstitu;
config.use_institu = [config.givInstitu config.ensInstitu];

%% Shared roots and independent configuration case
config.resultRoot = fullfile(config.ifilesRoot, 'Results_STPCGA1_SCS');
config.figureRoot = fullfile(config.ifilesRoot, 'Figure_STPCGA1_SCS');
config.caseName   = 'paper_default';

%% Execution options
config.Smooth       = ["None","Gaussian300km","DDK3","DDK5"];
config.redo         = true;
config.plotProcess  = true;
config.saveAddData  = true;
config.note         = '';

results = run_case_fast(config);
disp([config.areaName ' STPC-GA case finished.']);
toc
