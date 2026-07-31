function varargout=TWSt2slept_model(Dataproduct,TH,XY_buffer,Lwindow,forcenew)
% [slepcoffs,thedates,Coeff,TH,G,CC,V,N]=TWSt2slept(Dataproduct,TH,XY_buffer,Lwindow,masked,forcenew)
%
% This program reads in the Terrestrial Water Storage we have modeled and 
% converts to a Slepian basis for analysis.
%
% INPUT:
% 
% Dataproduct   What version of model you want (ii,oo,rr). Such as:
%                 'GLDAS1-NOAH-M-1-masked'  [default]
%                 'GLDAS1-NOAH-M-025-masked'
%                 etc.
% TH         The concentration region:
%              'england', 'eurasia',  'namerica', 'australia', 'greenland'
%              'africa', 'samerica', 'amazon', 'orinoco' OR
%              [lon lat] an ordered list defining a closed curve [degrees]
% XY_buffer  Distance in degrees that the region outline will be enlarged
%              by BUFFERM [default: 0]
% Lwindow    Bandwidth of the window [default: bandwidth of the data], 
%              or bandpass (give two degrees)
% forcenew   Wether or not you want to force new generation of a save file
%               (1) or just use the one we already have (0) [default].
%
% OUTPUT:
%
% Returns these variables and saves the first two in a .mat 
% file (kernels are not saved as they already should be elsewhere):
%
% slepcoffs    The expansion coefficients of the geopotential (or surface 
%                   density) into the Slepian basispotential Slepian coefficients 
%                   [nmonths x addmoff(Ldata)]
% thedates     Time stamps in Matlab time
% Coeff        The settings of inner and outer model.
% TH           The region back to you.  If there was buffering, this will
%                   be a XY array of coordinates, which you can use with 
%                   SPHAREA to get the Shannon number.
% G            The unitary matrix of localization coefficients
% CC           A cell array with cosine/sine coefficients eigenfunctions
% V            The eigenvalues in this ordering
% N            The Shannon number
%
%
% EXAMPLE: 
%
% See also: GRACE2SLEPT, READGLDAS, GETGLACIERS, MASKFROMGMT, COARSENMASK,
%           GLDAS2TWST
%
% Last modified by charig-at-princeton.edu, 03/16/2016

% Set some defaults
defval('Dataproduct',{'1','1'});
defval('TH','greenland');
defval('Lwindow',60);
defval('XY_buffer',0.5);
defval('forcenew',0);

% Figure out if it's lowpass or bandpass
lp=length(Lwindow)==1;
bp=length(Lwindow)==2;
maxL=max(Lwindow);
% The spherical harmonic dimension
ldim=(Lwindow(2-lp)+1)^2-bp*Lwindow(1)^2;
defval('J',ldim)

% Generated model data and Slepian expansions belong to the writable
% project data root, not to the source directory.
projectRoot = getenv('IFILES');
if isempty(projectRoot)
    error('STPCGA:ProjectRootNotConfigured', ...
        'Run startup_STPC_GA before creating simulation models.');
end
defval('ddir1',fullfile(projectRoot,'MOD'));
defval('ddir2',fullfile(projectRoot,'MOD','SlepianExpansions'));

if exist(ddir1, 'dir') ~= 7
    mkdir(ddir1);
end
if exist(ddir2, 'dir') ~= 7
    mkdir(ddir2);
end

% GEOGRAPHICAL REGIONS and XY REGIONS
% Check if you want the Shannon number of eigenfunctions
if strcmp(J,'N')
   J = round((Lwindow+1)^2*spharea(TH));
end
if isstr(TH) % Geographic
    % Here, TH gets passed to glmalpha, and glmalpha will interpret
    % either the cell of the region
    TH_origin=TH;
    if XY_buffer ~= 0
        TH = {TH XY_buffer};
        h = [TH{1} num2str(XY_buffer)];
    else
        h=TH;
    end
else % Closed coordinates (make a hash)
    h=hash(TH,'sha1');
end
% The name of the save file
if length(Dataproduct)==2
    if lp
        fnpl=sprintf('%s/TWSt2slept-I%sO%s-%s-%i-%i.mat',...
            ddir2,Dataproduct{1},Dataproduct{2},h,Lwindow,J);

        % notice that this has been changed due to the usage of land
%         fnmod=sprintf('%s/MODEL-I%sO%s-%s-%i.mat',...
%             ddir1,Dataproduct{1},Dataproduct{2},TH_origin,Lwindow);
        fnmod=sprintf('%s/MODEL-I%sO%s-%s-%i.mat',...
            ddir1,Dataproduct{1},Dataproduct{2},h,Lwindow);
    elseif bp
        error('Not feasible by now. [MA]')
        fnpl=sprintf('%s/%s-%s-%s-%i-%i-%i.mat',...
            ddir2,Dataproduct,h,Lwindow(1),Lwindow(2),J);
    else
        error('The degree range is either one or two numbers')
    end

elseif length(Dataproduct)==3

        if lp
        fnpl=sprintf('%s/TWSt2slept-I%sO%sR%s-%s-%i-%i.mat',...
            ddir2,Dataproduct{1},Dataproduct{2},Dataproduct{3},h,Lwindow,J);

        % notice that this has been changed due to the usage of land
%         fnmod=sprintf('%s/MODEL-I%sO%sR%s-%s-%i.mat',...
%             ddir1,Dataproduct{1},Dataproduct{2},Dataproduct{3},TH_origin,Lwindow);
        fnmod=sprintf('%s/MODEL-I%sO%sR%s-%s-%i.mat',...
            ddir1,Dataproduct{1},Dataproduct{2},Dataproduct{3},h,Lwindow);
    elseif bp
        error('Not feasible by now. [MA]')
        fnpl=sprintf('%s/%s-%s-%s-%s-%i-%i-%i.mat',...
            ddir2,Dataproduct,h,Lwindow(1),Lwindow(2),J);
    else
        error('The degree range is either one or two numbers')
        end

end

% GET THE SLEPIAN BASIS WE WANT
[~,~,~,lmcosiW,~,~,~,~,~,ronmW]=addmon(maxL);
% NOTE: a) Could have used PLM2SLEP but since we need the same operation for many
%         months, slightly better to load the Slepian basis once, and repeatedly
%         multiply by G.  Mostly copied from PLM2SLEP.
%       b) The kernel and eigenfunctions are large to save for large Lwindow.  
%          So instead just load them from either GLMALPHA or GLMALPHAPTO.
%       c) GLMALPHA now handles bandpass

% If it is the standard North-Polar cap or a geographic region, it's easy
[G,V,EL,EM,N,GM2AL,MTAP,IMTAP]=glmalpha(TH,Lwindow,[],0,[],[],J);
% Sort by decreasing eigenvalue
[V,vi]=sort(V,'descend');
G=G(:,vi); if ~isnan(MTAP); MTAP=MTAP(vi); end
% If you don't do this, the eigenfunctions are ordered in the way
%   that they correspond to single-orders back when, unrotated, they
%   belonged to a polar cap, and the eigenvalues are sorted within
%   these blocks. This is useful for, e.g. SPIE2009_1 a la SDSNEEUW. 
% Collect the eigenvector output into a format that PLM2XYZ knows how to interpret
for j=1:size(G,2)
   % Create the blanks
   cosi=lmcosiW(:,3:4);
   % Stick in the coefficients of the 1st eigentaper
   cosi(ronmW)=G(:,j);
   % Construct the full matrix
   CC{j} = [lmcosiW(:,1:2) cosi]; 
end

%%%
% Setup complete
%%%

% If this expansion already exists, load it.  Otherwise, or if we force 
% it, make a new one (e.g. if you added extra months to the database).
if exist(fnpl,'file')==2 && forcenew==0
     load(fnpl)
     disp(sprintf('%s loaded by TWST2SLEPT',fnpl))
else

    if exist(fnmod,'file')==2 && forcenew==0
        load(fnmod);
    else
        % notice that this has been changed due to the usage of land
        if length(Dataproduct)==2
%             [thedates,TWSlmcosi,Coeff]=MODEL_Make_m(Dataproduct,TH_origin,Lwindow);
            [thedates,TWSlmcosi,Coeff]=MODEL_Make_m(Dataproduct,TH_origin,Lwindow,XY_buffer);
        elseif length(Dataproduct)==3
%             [thedates,TWSlmcosi,Coeff,Check_board]=MODEL_Checkboard_Make_m(Dataproduct,TH_origin,Lwindow);
            [thedates,TWSlmcosi,Coeff,Check_board]=MODEL_Checkerboard_Make_new_m(Dataproduct,TH_origin,Lwindow,XY_buffer);
        else
            error('Wrong input format of Dataproduct')
        end
    end

%     [version,forcing,timeres,degres] = strread(Dataproduct,'%s%s%s%n','delimiter','-');
%     % Use GLDAS2TWST to get the Terrestrial Water Storage
%     [thedates,TWSlmcosi] = gldas2TWSt(version{1},degres,timeres{1},forcing{1},masked,forcenew);
    
    % Initialize new coefficients
    nmonths=length(thedates);
    slepcoffs=nan(nmonths,J);
    
    % Limit everything to the window bandwidth
    TWSlmcosiW = cellfun(@(x) x(1:size(lmcosiW,1),1:4),TWSlmcosi,'UniformOutput',false);
    
    % Loop over the months
    for index = 1:nmonths
        % Expand this month's TWS into the Slepian basis
        %potcoffs_month=squeeze(potcoffsW(index,:,:));
        slepcoffs(index,:) = ...
            TWSlmcosiW{index}(2*size(TWSlmcosiW{index},1)+ronmW(1:(maxL+1)^2))'*G;
    end
        
    % SAVE
    % Don't save the kernel and eigenfunctions because we already have 
    % this info saved and can load from GLMALPHA
    % save(fnpl,'slepcoffs','thedates','G','CC','V');
    if length(Dataproduct)==2
        save(fnpl,'slepcoffs','thedates','Coeff');
    elseif length(Dataproduct)==3
        save(fnpl,'slepcoffs','thedates','Coeff','Check_board');
    end

end % End if we have a save file already
    


% Collect output
if length(Dataproduct)==2
    varns={slepcoffs,thedates,Coeff,TH,G,CC,V,N};
elseif length(Dataproduct)==3
    varns={slepcoffs,thedates,Coeff,TH,G,CC,V,N,Check_board};
end
varargout=varns(1:nargout);


