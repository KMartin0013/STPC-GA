function varargout=grace2plmt_m(Pcenter,Rlevel,Ldata,units,forcenew,land_or_ocean,GIA)
% [potcoffs,cal_errors,thedates]=GRACE2PLMT(Pcenter,Rlevel,units,forcenew)
%
% This program reads in the Level-2 GRACE geoid products from either the CSR or
% GFZ data centers, does some processing, and saves them as a plmt matrix
% in a .mat file.  In particular, the coefficients are reordered to our
% prefered lmcosi format, they are referenced to the WGS84 ellipsoid,
% the C2,0 coefficients are replaced with more accurate measurements from
% satellite laser ranging from Loomis et al, (2020), and the degree one coefficients are
% substituted with those from Sun et al., (2016).  You have the option
% of leaving them as geopotential
% or converting them to surface mass density using the method of
% Wahr et al. 1998, based on Love numbers (see PLM2POT).
%
% INPUT:
%
% Pcenter     'CSR' data center at the Center for Space Research
%             'GFZ' data center at the GeoForschungsZentrum Potsdam
% Rlevel      The release level of the solution you want.
%              Either 'RL04','RL05', or 'RL06'
% Ldata       The bandwidth of the dataproduct that you want [default: 60].
%              In the case where there are more than one product from a
%              datacenter (such as BA 60 or BB 96 standard L2 products)
%              this allows you to choose between them.
% units       'POT' or 'SD' for whether you want geopotential or surface
%               mass density
% forcenew    Whether or not you want to force new generation of a save file
%              (1) or just use the one we already have (0) [default].
% land_or_ocean
%             Whether the data was used on land or ocean. This will
%             change the process on degree 1
% GIA
%             Whether you want to remove GIA (e.g., 'Pelt17' for RL06)
%             or not [default].
%
% OUTPUT:
%
% Returns these variables and saves them in a .mat file:
%    potcoffs       potential coefficients [nmonths x addmup(Ldata) x 4]
%                    these could also be in surface mass density
%    thedates       time stamps in Matlab time
%
% NOTE:
%
%   10/15/2023 Added coefficient 'land_or_ocean' for differenrt applications
%   with different process on Degree 1 (Refer to GRACE Technical Note 13)
%
%   5/18/2022 Added L to the inputs so that we can utilize more than one
%    product from a data center. A corresponding change has been made in
%    GRACE2SLEPT
%
%   2/19/2021 Formal or calibrated uncertainties have not been reported
%    since RL04 so we will discontinue the output of these. The output
%    variables will be altered so existing scripts will need to be adjusted
%    to take this into account. A corresponding change has been made in
%    GRACE2SLEPT
%
%	GRACE data available from NASA PODAAC:
%	https://podaac.jpl.nasa.gov/
%
%   See gracedeg1.m and gracedeg2.m to see how we update the spherical
%   harmonic degree 1, degree 2, and degree 3 (when available)
%   coefficients. These are usually read from the Technical Notes documents
%   distributed alongside GRACE data.
%
%
%
% EXAMPLE: to make a new save file when you have added more months
% [potcoffs,thedates]=grace2plmt('CSR','RL06',60,'SD',1);
%
% Last modified by zhongtian.ma@connect.polyu.hk, 10/15/2023
% Last modified by charig-at-email.arizona.edu, 5/18/2022
% Last modified by lashokkumar-at-arizona.edu, 11/09/2020
% Last modified by mlubeck-at-email.arizona.edu, 03/18/2019
% Last modified by fjsimons-at-alum.mit.edu, 05/17/2011





% Determine parameters and set defaults
defval('Pcenter','CSR')
defval('Rlevel','RL06')
defval('Ldata','60')
defval('units','SD')
defval('forcenew',1)
defval('land_or_ocean','land')
defval('GIA','No')

% Where the original data files are kept
defval('ddir1',fullfile(getenv('IFILES'),'GRACE','Originals',Rlevel,Pcenter));

% Where you would like to save the new .mat file
defval('ddir2',fullfile(getenv('IFILES'),'GRACE'));
% And the name of that save file
if strcmp(units,'SD')
    fnpl=sprintf('%s/%s_%s_alldata_%s_%s_%s_%s.mat',ddir2,Pcenter,Rlevel,num2str(Ldata),units,land_or_ocean,GIA(1:2));
elseif strcmp(units,'POT')
    fnpl=sprintf('%s/%s_%s_alldata_%s_%s_%s.mat',ddir2,Pcenter,Rlevel,num2str(Ldata),land_or_ocean,GIA(1:2));
else
    error ('The unit input is not valid.')
end

%m in case different time span of data is used
if length(Pcenter)>=3
    Pcenters=Pcenter(1:end-4);
end

% Calculate the IB correction for specific products
% It is officially provided by GRACE as TN-12.
% However, no TN-12 is provided in GRACE-FO, so it is necessary now.
% GAA was used for GRGS, while GAD was used for other institutions

% Where you have the land-ocean mask file
defval('ddir3',fullfile(getenv('IFILES'),'Data'));
fnpl_IB=sprintf('%s/%s_%s_IB_%s_%s_%s.mat',ddir2,Pcenter,Rlevel,num2str(Ldata),units,'GSHHC');

c11cmn2=[-179.5 89.5 179.5 -89.5];
lon=c11cmn2(1):c11cmn2(3);lat=c11cmn2(2):-1:c11cmn2(4);
[lonlon,latlat] = meshgrid(lon,lat);

load(fullfile(ddir3, 'coast_GSHHS_c_globe'))
[land_in,land_on] = inpolygon(lonlon,latlat,long,lat);

leakage_sp=[];

% If this file already exists, load it.  Otherwise, or if we force it, make
% a new one (e.g. you added extra months to the database).
if exist(fnpl,'file')==2 && (strcmp(land_or_ocean,'land') || exist(fnpl_IB,'file')==2) && forcenew==0
    load(fnpl)
    disp(sprintf('%s loaded by GRACE2PLMT',fnpl))
else
    if ~exist(ddir1,'dir')
        error ('The data you asked for are not currently stored.')
        % If you received this error it means the directory datanames does
        % not exist which means you have not downloaded the data, or the data
        % does not exist where the function expects.
    end

    % DATA CENTER
    datanames_D = {};
    Ldata_D = Ldata;
    ddir1_D = ddir1;


    if strcmp(Pcenters,'GFZ')
        if strcmp(Rlevel,'RL04')
            % Find the coefficient files
            datanames=ls2cell(fullfile(ddir1,'GSM*G---_0004'));
            % Find the error files
            errornames=ls2cell(fullfile(ddir1,'GSM*G---_0004.txt'));
            % Know a priori what the bandwidth of the coefficients is
            Ldata=120;
        elseif strcmp(Rlevel,'RL05')
            % Find the coefficient files
            datanames=ls2cell(fullfile(ddir1,'GSM*G---_005a'));
            % Know a priori what the bandwidth of the coefficients is
            Ldata=90;
            datanames_D=ls2cell(fullfile(ddir1,'GAD*0005'));
            Ldata_D=100; % upperlimit of GAD
        elseif strcmp(Rlevel,'RL06')
            % 10/17/2023 GR-FO only now version 6.1
            if Ldata == 60
                datanames=[ls2cell(fullfile(ddir1,'GSM*BA01_060*'))];
            elseif Ldata == 96
                datanames=[ls2cell(fullfile(ddir1,'GSM*BB01_060*'))];
            else
                error(['Solutions with requested L=' num2str(Ldata) ' not currently stored']);
            end
            datanames_D=ls2cell(fullfile(ddir1,'GAD*060*'));
            Ldata_D=180;
        end
    elseif  strcmp(Pcenters,'CSR')
        if strcmp(Rlevel,'RL04')
            datanames=ls2cell(fullfile(ddir1,'GSM*0060_0004'));
            errornames=ls2cell(fullfile(ddir1,'GSM*0060_0004.txt'));
        elseif strcmp(Rlevel,'RL05')
            datanames=ls2cell(fullfile(ddir1,'GSM*0060_0005'));
            datanames_D=ls2cell(fullfile(ddir1,'GAD*0000_0005'));
            Ldata_D=100;
        elseif strcmp(Rlevel,'RL06')
            % 5/12/2022 GR-FO only now version 6.1
            % 10/26/2024 GR-FO only now version 6.3 %m
            if Ldata == 60
                datanames=[ls2cell(fullfile(ddir1,'GSM*BA01_06*'))];
            elseif Ldata == 96
                datanames=[ls2cell(fullfile(ddir1,'GSM*BB01_06*'))];
            else
                error(['Solutions with requested L=' num2str(Ldata) ' not currently stored']);
            end
            % Naming convention was changed for RL06 where BA stands
            % for degree 60 gravity solution and 01 represents unconstrained
            % spherical harmonic solution with a boxcar windowing function
            % (see L-2 UserHandbook_v4.0).
            % In addition, the BB stands for degree 96 gravity solution with a
            % boxcar windowing function.
            % The other change was made to the last naming entry, which is now
            % in the form rrvv. In this case rr represents the release number
            % maximum 2 digits and vv represents the maximum 2 digit version
            % number. So for RL06 the nomenclature is 0600 instead of 0006 for
            % Rl05 previosly. The PID naming convention stays the same.
            datanames_D=ls2cell(fullfile(ddir1,'GAD*_06*'));
            Ldata_D=180;
        end
        % Know a priori what the bandwidth of the coefficients is
        % Ldata=60;
    elseif  strcmp(Pcenters,'JPL')
        if strcmp(Rlevel,'RL05')
            %JPL has 0001_0005 and 0000_0005
            datanames=ls2cell(fullfile(ddir1,'GSM*JPLEM*0005'));
            % JPL Release Level 5 has no calibrated error files
            %errornames=ls2cell(fullfile(ddir1,'GSM*0060_0004.txt'));
            %JPL has 0001_0005 and 0000_0005
            datanames_D=ls2cell(fullfile(ddir1,'GAD*0005'));
            Ldata_D=100;
        elseif strcmp(Rlevel,'RL06')
            % 10/17/2023 GR-FO only now version 6.1
            if Ldata == 60
                datanames=[ls2cell(fullfile(ddir1,'GSM*BA01_060*'))];
            elseif Ldata == 96
                datanames=[ls2cell(fullfile(ddir1,'GSM*BB01_060*'))];
            else
                error(['Solutions with requested L=' num2str(Ldata) ' not currently stored']);
            end
            datanames_D=ls2cell(fullfile(ddir1,'GAD*060*'));
            Ldata_D=180;
        else
            error('JPL RL04 solutions not currently stored');
            %elseif strcmp(Rlevel,'RL05');
            %    datanames=ls2cell(fullfile(ddir1,'GSM*JPLEM*005'));
        end
        % Know a priori what the bandwidth of the coefficients is
        %    Ldata=90;
    elseif  strcmp(Pcenters,'ITSG')
        if strcmp(Rlevel,'RL05')
            error('Not available.')
        elseif strcmp(Rlevel,'RL06')
            if Ldata == 60
                datanames=[ls2cell(fullfile(ddir1,'ITSG*n60*'))];
            elseif Ldata == 96
                datanames=[ls2cell(fullfile(ddir1,'ITSG*n96*'))];
            else
                error(['Solutions with requested L=' num2str(Ldata) ' not currently stored']);
            end
            datanames_D=ls2cell(fullfile(ddir1,'model_oceanBottom*'));
            Ldata_D=180;
        else
            error('JPL RL04 solutions not currently stored');
            %elseif strcmp(Rlevel,'RL05');
            %    datanames=ls2cell(fullfile(ddir1,'GSM*JPLEM*005'));
        end
        % Know a priori what the bandwidth of the coefficients is
        %    Ldata=90;
    elseif  strcmp(Pcenters,'COSTG')
        if strcmp(Rlevel,'RL06')
            if Ldata > 90
                error(['COSTG RL06 files currently stored here are degree 90; requested L=' num2str(Ldata)]);
            end
            datanames=ls2cell(fullfile(ddir1,'GSM-2_*COSTG*.gfc'));
            if isempty(datanames)
                error('No COSTG .gfc files were found. Expected files such as GSM-2_*COSTG*.gfc.');
            end
            ddir1_D=fullfile(getenv('IFILES'),'GRACE','Originals',Rlevel,'ITSG0324');
            datanames_D=ls2cell(fullfile(ddir1_D,'model_oceanBottom*'));
            Ldata_D=180;
        else
            error('COSTG support is currently implemented for RL06 only.');
        end
    elseif  strcmp(Pcenters,'Tongji')
        if strcmp(Rlevel,'RL06')
            if Ldata > 60
                error(['Tongji0323 files currently stored here are degree 60; requested L=' num2str(Ldata)]);
            end
            datanames=ls2cell(fullfile(ddir1,'Tongji-Grace2022_n60_*.gfc'));
            if isempty(datanames)
                error('No Tongji .gfc files were found. Expected files such as Tongji-Grace2022_n60_YYYY-MM.gfc.');
            end
            ddir1_D=fullfile(getenv('IFILES'),'GRACE','Originals',Rlevel,'ITSG0324');
            datanames_D=ls2cell(fullfile(ddir1_D,'model_oceanBottom*'));
            Ldata_D=180;
        else
            error('Tongji support is currently implemented for RL06 only.');
        end
    elseif  strcmp(Pcenters,'GRGS')
        if strcmp(Rlevel,'RL05')
            %GRGS only has TSVD_0005, here only 90 degrees are considered
            datanames=ls2cell(fullfile(ddir1,'GSM*TSVD_0500.txt'));
            %GRGS only has GAB*_0090_0005
            datanames_D=ls2cell(fullfile(ddir1,'GAB*0090_0005.txt'));
            Ldata_D=90;
            %GRGS only has GAA*_0090_0005
            datanames_IB=ls2cell(fullfile(ddir1,'GAA*0090_0005.txt'));
            Ldata_IB=90;
        elseif strcmp(Rlevel,'RL06')
            error('No RL06 for now.')
        end
    end

    % WGS84 reference SETUP
    % For now just hardcode the even zonal coefficients (J), later use
    % Frederik's GRS.m program, don't bother with the higher degrees
    j2= 0.108262982131e-2*-1.0/(2*2+1)^0.5; % will be row 4
    j4=-0.237091120053e-5*-1.0/(2*4+1)^0.5; % will be row 11
    % Also useful
    a=fralmanac('a_EGM96','Earth');

    % C20 CORRECTION SETUP

    % % Load the C(2,0) coefficients from satellite laser ranging, depending on
    % % our release level. Note, here the 'NH' part describes a no header version
    % % of the SLR datafiles.
    % if strcmp(Rlevel,'RL04')
    %     slrc20=load(fullfile(getenv('IFILES'),'SLR','C20_RL04_MA.txt'));
    % elseif strcmp(Rlevel,'RL05')
    %     slrc20=load(fullfile(getenv('IFILES'),'SLR','C20_RL05_MA.txt'));
    % elseif strcmp(Rlevel,'RL06')
    %     slrc20=load(fullfile(getenv('IFILES'),'SLR','C20_RL06_MA.txt'));
    % end
    % % The sigma error is column 4
    % slrc20_error=slrc20(:,4)*1e-10;
    % % Remove the AOD1B model which was removed from the GRACE GSM data but
    % % restored to this SLR data.  Use the raw value (column 2). See data headers.
    % slrc20=[slrc20(:,1) slrc20(:,2)-slrc20(:,5)*1e-10];
    % % Convert the dates to Matlab format
    % [n,m]=size(slrc20);
    % slrc20(:,1)=datenum([slrc20(:,1) ones(n,1) ones(n,1)]);
    % % Make slrc20 relative to the WGS84 ellipsoid
    % slrc20(:,2) = slrc20(:,2) - j2;

    % New function that returns you the coefficients for degree L = 2 (and 3
    % when available). Format is 3 columns: thedates, C20, C30. thedates are
    % the midpoint of the GRACE data month.
    % m this is also OK for 'RL05'
    [slrc20] = gracedeg2(Rlevel);

    % These are not referenced to anything, so make the 2,0 coefficient
    % relative to the WGS84 ellipsoid
    slrc20(:,2) = slrc20(:,2) - j2;


    % Degree 1 Correction Setup

    % if Rlevel=='RL04'
    %     deg1=load(fullfile(getenv('IFILES'),'GRACE','deg1_RL04_NH.txt'));
    % elseif Rlevel=='RL05'
    %     deg1=load(fullfile(getenv('IFILES'),'GRACE','deg1_RL05_NH.txt'));
    % elseif Rlevel=='RL06'
    %     %deg1=load(fullfile(getenv('IFILES'),'GRACE','deg1_RL06_new.txt'));  % Technical note 11
    %     deg1=load(fullfile(getenv('IFILES'),'GRACE','deg1_RL06_Sun2016.txt'));  % Technical note 13
    % end
    % [n,m] = size(deg1);
    % dates_str = num2str(deg1(:,1));
    % deg1dates = datenum([str2num(dates_str(:,1:4)) str2num(dates_str(:,5:6)) 15*ones(n,1)]);
    % [b,m] = unique(deg1dates);
    % deg1dates = deg1dates(m);
    % for i=1:n/2; temp = [deg1(2*i-1,2:7); deg1(2*i,2:7)]; mydeg1(i,:,:) = temp; end;

    % New function that returns you the coefficients for degree L = 1. Format
    % should be identical to old code above.
    [deg1dates,mydeg1] = gracedeg1_m(Pcenters,Rlevel,'land');

    %m New function that adds the degree 1 of GAD products for ocean
    %application
    % The only possible usage of this separate maybe in study in river
    % plume, where GAD higher than degree 1 should be removed to
    % study the residual oceanic signal (the only difference to directly
    % use TWS is the additional degree 1 of GAD).

    % if strcmp(land_or_ocean,'ocean')
    %     [deg1dates_oce,mydeg1_oce] = gracedeg1_m(Pcenter,Rlevel,land_or_ocean);
    % end

    % Initialize
    nmonths = length(datanames);
    thedates = zeros(1,nmonths);
    [dems,dels]=addmon(Ldata);
    [dems_D,dels_D]=addmon(Ldata_D);

    % Calibrated errors are normally used instead, but they are kept here for
    % completeness.

    % Last two columns here are "formal" errors
    % l m cosine sine cosine_stddev sine_stddev
    potcoffs=nan(nmonths,addmup(Ldata),6);
    % Last two columns here are "calibrated" errors
    % l m cosine sine
    cal_errors=nan(nmonths,addmup(Ldata),4);

    % Loop over the months
    % special for cases when the number of GAD is larger than the number of
    % GSM; this way we can directly get the leakage number
    index_D=0;
    for index = 1:nmonths
        disp(['Processing month number: ' num2str(index)]);

        % load geopotential coefficients
        fname1=fullfile(ddir1,datanames{index});

        if strcmp(Pcenters,'ITSG') || strcmp(Pcenters,'COSTG') || strcmp(Pcenters,'Tongji')

            % Open and scan the file (data from ITSG is 7 columns)
            fid = fopen(fname1);
            C = textscan(fid,'%s%s%s%s%s%s%s');
            fclose(fid);

            % Only grab the lines for gfc
            Carray = cat(3,C{:});
            I = strmatch('gfc',Carray(:,1,1),'exact');
            Carray = squeeze(Carray(I,1,:));

        else

            % Open and scan the file (data from all three centers is 10 columns)
            fid = fopen(fname1);
            C = textscan(fid,'%s%s%s%s%s%s%s%s%s%s');
            fclose(fid);

            % Only grab the lines for GRCOF2
            Carray = cat(3,C{:});
            I = strmatch('GRCOF2',Carray(:,1,1),'exact');
            Carray = squeeze(Carray(I,1,:));

        end

        % Only want columns 2-7, and as format double
        Carray = Carray(:,2:7);
        lmcosi_month=cellfun(@str2num,Carray);
        if strcmp(Pcenters,'COSTG') || strcmp(Pcenters,'Tongji')
            lmcosi_month = lmcosi_month(lmcosi_month(:,1) <= Ldata, :);
        end
        % This should be addmup(Ldata)
        [m,n] = size(lmcosi_month);

        if strcmp(Pcenters,'GFZ') || strcmp(Pcenters,'CSR')
            % Change the order of the coefficients so that
            % order m goes as [0 01 012 0123 ...]
            new_ordering = zeros(m,6);
            revdel=[0 Ldata:-1:0];
            i=1;
            for j=1:length(dems)
                k = dels(i)+1 + sum( revdel( (1:dems(i) + 1 ) ) );
                new_ordering(j,:) = lmcosi_month(k,:);
                i=i+1;
            end
            lmcosi_month = new_ordering;
        elseif strcmp(Pcenters,'JPL')
            % The JPL coefficients are in the order we want already; just need
            % to add the zero and one coefficients
            lmcosi_month = [0 0 0 0 0 0; 1 0 0 0 0 0; 1 1 0 0 0 0; lmcosi_month];
            % Some JPL months are only degree 60, so add in zeros to get up to
            % degree 90
            ldim=length(lmcosi_month);
            if max(lmcosi_month(:,1))==60
                [~,~,~,lmcosiE]=addmon(Ldata);
                lmcosi_month = [lmcosi_month;...
                    lmcosiE(ldim+1:end,1:2) zeros(length(lmcosiE)-ldim,4)];
            end

        elseif strcmp(Pcenters,'GRGS') || strcmp(Pcenters,'ITSG') || strcmp(Pcenters,'COSTG') || strcmp(Pcenters,'Tongji')
            % The GRGS and ITSG coefficients are in the order we want already
            ldim=length(lmcosi_month);
        else
            error('Data center error')
        end

        if  strcmp(Pcenters,'ITSG')
            % Calculate the midpoint of this data span
            monthstart = datenum([str2num(datanames{index}(end-10:end-7)) ...
                str2num(datanames{index}(end-5:end-4)) 1]);
            monthend = datenum([str2num(datanames{index}(end-10:end-7))+floor(str2num(datanames{index}(end-5:end-4))/12) ...
                mod(str2num(datanames{index}(end-5:end-4)),12)+1 1])-1;
            monthmid = (monthstart+monthend)/2;
            thedates(index) = monthmid;

        elseif strcmp(Pcenters,'COSTG')
            dateToken = regexp(datanames{index}, 'GSM-2_(\d{4})(\d{3})-(\d{4})(\d{3})', 'tokens', 'once');
            if isempty(dateToken)
                error(['Cannot parse COSTG date from file name: ' datanames{index}]);
            end
            monthstart = datenum(str2double(dateToken{1}), 1, str2double(dateToken{2}));
            monthend = datenum(str2double(dateToken{3}), 1, str2double(dateToken{4}));
            monthmid = (monthstart + monthend) / 2;
            thedates(index) = monthmid;

        elseif strcmp(Pcenters,'Tongji')
            dateToken = regexp(datanames{index}, '_(\d{4})-(\d{2})\.gfc$', 'tokens', 'once');
            if isempty(dateToken)
                error(['Cannot parse Tongji date from file name: ' datanames{index}]);
            end
            yearNum = str2double(dateToken{1});
            monthNum = str2double(dateToken{2});
            monthstart = datenum(yearNum, monthNum, 1);
            monthend = datenum(yearNum + floor(monthNum/12), mod(monthNum,12) + 1, 1) - 1;
            monthmid = (monthstart + monthend) / 2;
            thedates(index) = monthmid;

        else
            % Calculate the midpoint of this data span
            monthstart = datenum([str2num(datanames{index}(7:10))...
                1 str2num(datanames{index}(11:13))]);
            monthend = datenum([str2num(datanames{index}(15:18))...
                1 str2num(datanames{index}(19:21)) ]);
            monthmid = (monthstart+monthend)/2;
            thedates(index) = monthmid;

        end

        if  ~strcmp(Pcenters,'GRGS')
            disp(['Month midpoint date: ' datestr(thedates(index))]);
        else
            disp(['10 days midpoint date: ' datestr(thedates(index))]);
        end

        % Remove the mean value of the potential i.e. set 0,0 coff = 0
        lmcosi_month(1,3) = 0;

        % GRGS product do not need correct the degree 1 and 2
        if ~strcmp(Pcenters,'GRGS')

            % Make the geopotential relative to the WGS 84 ellipsoid
            % A bit redundant since we replace 2,0 shortly
            lmcosi_month(4,3) = lmcosi_month(4,3) - j2;
            lmcosi_month(11,3) = lmcosi_month(11,3) - j4;

            % Now replace the (2,0) and (3,0) (when available) coefficient with
            % the SLR value (referenced to WGS84 above).
            %     % NOTE: This gives a value different than if you used
            %     % (column3 - column5) from the SLR data file because that data is
            %     % referenced to an overall mean, not to the WGS 84 ellipsoid.
            %     %if index==111; keyboard; end;
            where=slrc20(:,1)>monthstart & slrc20(:,1)<monthend;
            if ~any(where)
                % If there is no SLR value within our specific interval,
                % use the closest value we have
                [~,where]=min(abs(monthmid - slrc20(:,1)));
                disp('No SLR coefficients for this month, used nearest available.')
                pause(2)
            elseif nnz(where) > 1
                % We have more than one month. Use the closest value
                [~,where]=min(abs(monthmid - slrc20(:,1)));
                disp(['More than one degree 2 value detected, likely from overlapping months. Used nearest: ' datestr(slrc20(where,1))])
                pause(1)
            end
            % else we just use "where" which only has 1 valid entry

            % Need to use slrc20(where,2)
            disp(sprintf('C20 was %12.8e now %12.8e',lmcosi_month(4,3),slrc20(where,2)))
            lmcosi_month(4,3)=slrc20(where,2);
            % Now replace C3,0 if possible
            if isfinite(slrc20(where,3))
                % Here we have a C3,0 SLR value available, so substitute.
                disp(sprintf('C30 was substituted; it was %12.8e now %12.8e',lmcosi_month(7,3),slrc20(where,3)))
                lmcosi_month(7,3)=slrc20(where,3);
            end


            % Now replace the degree 1 coefficients with those from GRACEDEG1.m and
            % references therin.

            % Find a degree 1 data point that is within the month of our GRACE data
            where1=deg1dates(:)>monthstart & deg1dates(:)<monthend;
            if ~any(where1)
                % If there is no Deg1 value within our specific interval,
                % don't change anything, because we know a few months are missing
                disp('No change to degree 1')
                pause(2)
            elseif nnz(where1) > 1
                % We have more than one month. Use the closest value
                disp('More than one degree 1 value detected')
                [~,where1]=min(abs(monthmid - deg1dates));
                lmcosi_month(2:3,1:4)=squeeze(mydeg1(where1,:,1:4));
                disp(['Deg1 value for ' datestr(deg1dates(where1)) ' used.']);
                pause(1)
            else
                disp(['Deg1 value for ' datestr(deg1dates(where1)) ' used.']);
                lmcosi_month(2:3,1:4)=squeeze(mydeg1(where1,:,1:4));
            end

        else
            disp(['Deg1 and Deg2 values do not need to be replaced for GRGS products.']);
        end

        if strcmp(land_or_ocean,'ocean')

            if isempty(datanames_D)
                error(['No atmosphere/ocean background product is available for ' Pcenters '. COSTG/Tongji have no own GAD/GAC files; check ITSG0324 model_oceanBottom files for ocean cases.']);
            end

            if strcmp(Pcenters,'COSTG') || strcmp(Pcenters,'Tongji')
                disp(['Ocean case for ' Pcenters ': using ITSG model_oceanBottom files for GAD-style restore.']);
            end
            where2=0;

            while ~where2 && index_D < length(datanames_D)
                index_D=index_D+1;

                % load geopotential coefficients of GAD
                fname3=fullfile(ddir1_D,datanames_D{index_D});

                if  strcmp(Pcenters,'ITSG') || strcmp(Pcenters,'COSTG') || strcmp(Pcenters,'Tongji')
                    % Calculate the midpoint of this data span
                    monthstart_D = datenum([str2num(datanames_D{index_D}(end-10:end-7)) ...
                        str2num(datanames_D{index_D}(end-5:end-4)) 1]);
                    monthend_D = datenum([str2num(datanames_D{index_D}(end-10:end-7))+floor(str2num(datanames_D{index_D}(end-5:end-4))/12) ...
                        mod(str2num(datanames_D{index_D}(end-5:end-4)),12)+1 1])-1;
                    monthmid_D = (monthstart_D+monthend_D)/2;
                    thedates_D(index_D) = monthmid_D;

                else

                    % Calculate the midpoint of this data span
                    monthstart_D = datenum([str2num(datanames_D{index_D}(7:10))...
                        1 str2num(datanames_D{index_D}(11:13))]);
                    monthend_D = datenum([str2num(datanames_D{index_D}(15:18))...
                        1 str2num(datanames_D{index_D}(19:21))]);
                    monthmid_D = (monthstart_D+monthend_D)/2;
                    thedates_D(index_D) = monthmid_D;

                end

                %m Add the degree (1~Lmax) of GAD for ocean usage
                where2=thedates_D(index_D)==thedates(index);

                if ~where2
                    leakage_sp=[leakage_sp; index_D thedates_D(index_D)];
                end
            end

            if ~any(where2)
                % If there is no Deg1 value within our specific interval,
                % it should be an error (since your GAD files should no less than your GSM files)
                error('No GAD (or GAB) product to corresponding GSM.')
            else

                if  strcmp(Pcenters,'ITSG') || strcmp(Pcenters,'COSTG') || strcmp(Pcenters,'Tongji')
                    % Open and scan the file (data from ITSG is 5 columns)
                    fid = fopen(fname3);
                    C = textscan(fid,'%s%s%s%s%s');
                    fclose(fid);

                    % Only grab the lines for GRCOF2
                    Carray = cat(3,C{:});
                    I = strmatch('gfc',Carray(:,1,1),'exact');
                    Carray = squeeze(Carray(I,1,:));

                    % Only want columns 2-7, and as format double
                    Carray = Carray(:,2:5);

                    % the STD for GAD of ITSG use zero [default]
                    Carray=[Carray, cellstr(num2str(zeros(size(Carray,1),1))), cellstr(num2str(zeros(size(Carray,1),1)))];

                else
                    % Open and scan the file (data from all three centers is 10 columns)
                    fid = fopen(fname3);
                    C = textscan(fid,'%s%s%s%s%s%s%s%s%s%s');
                    fclose(fid);

                    % Only grab the lines for GRCOF2
                    Carray = cat(3,C{:});
                    I = strmatch('GRCOF2',Carray(:,1,1),'exact');
                    Carray = squeeze(Carray(I,1,:));

                    % Only want columns 2-7, and as format double
                    Carray = Carray(:,2:7);
                end

                lmcosi_month_D=cellfun(@str2num,Carray);
                % This should be addmup(Ldata)
                [m,n] = size(lmcosi_month_D);

                if strcmp(Pcenters,'GFZ') || strcmp(Pcenters,'CSR')
                    % Change the order of the coefficients so that
                    % order m goes as [0 01 012 0123 ...]
                    new_ordering_D = zeros(m,6);
                    revdel=[0 Ldata_D:-1:0];
                    i=1;
                    for j=1:length(dems_D)
                        k = dels_D(i)+1 + sum( revdel( (1:dems_D(i) + 1 ) ) );
                        new_ordering_D(j,:) = lmcosi_month_D(k,:);
                        i=i+1;
                    end
                    lmcosi_month_D = new_ordering_D;
                elseif strcmp(Pcenters,'JPL')
                    % The JPL coefficients are in the order we want already;
                    % GAD products do NOT need to add the zero and one coefficients
                    %             lmcosi_month_D = [0 0 0 0 0 0; 1 0 0 0 0 0; 1 1 0 0 0 0; lmcosi_month_D];
                    % Some JPL months are only degree 60, so add in zeros to get up to
                    % degree 90

                    %actually this process is redundant because GAD for JPL is
                    %always 180 just like other institutions
                    ldim=length(lmcosi_month_D);
                    if max(lmcosi_month_D(:,1))==60
                        [~,~,~,lmcosiE]=addmon(Ldata_D);
                        lmcosi_month_D = [lmcosi_month_D;...
                            lmcosiE(ldim+1:end,1:2) zeros(length(lmcosiE)-ldim,4)];
                    end
                elseif strcmp(Pcenters,'GRGS') || strcmp(Pcenters,'ITSG') || strcmp(Pcenters,'COSTG') || strcmp(Pcenters,'Tongji')
                    % The GRGS and ITSG coefficients are in the order we want already;
                    ldim=length(lmcosi_month_D);
                else
                    error('Data center error')
                end

                % Remove the mean value of the potential i.e. set 0,0 coff = 0
                % A bit redundant since we use degrees higher than 0 of GAD
%                 lmcosi_month_D(1,3) = 0;

                disp(['Deg (1 ~ ' num2str(Ldata) ') value of GAD (or GAB) for ' ...
                    datestr(thedates_D(index_D)) ' used to restore ocean bottom pressure.']);
                lmcosi_month(2:size(lmcosi_month,1),3:4)=lmcosi_month(2:size(lmcosi_month,1),3:4)+...
                    lmcosi_month_D(2:size(lmcosi_month,1),3:4);

                % calcualte IB for institutions other than GRGS (which use GAA rather than GAD)
                if ~strcmp(Pcenters,'GRGS')
                    lmcosi_month_IB=lmcosi_month_D;Ldata_IB=Ldata_D;index_IB=index_D;
                    % Convert the geopotential coefficients into surface mass density
                    % (kg/m^2) 1 kg/m^2 /1000 kg/m^3 *1000 = 1 mm
                    % Need to make geoid first
                    lmcosi_extra_IB=plm2pot([lmcosi_month_IB(:,1:2) lmcosi_month_IB(:,5:6)*a],[],[],[],4);
                    lmcosi_month_IB=plm2pot([lmcosi_month_IB(:,1:2) lmcosi_month_IB(:,3:4)*a],[],[],[],4);
                    % Add the fornal errors back in columns 5,6
                    lmcosi_month_IB=[lmcosi_month_IB lmcosi_extra_IB(:,3:4)];

                    [IB,~,~,~,~]=plm2xyz(lmcosi_month_IB,1,c11cmn2,Ldata_IB);
                    %                 GAA(land_in)=NaN;
                    IB_C(index_IB,:)=[index_IB sum(sum(IB(~land_in)))/sum(sum(~land_in))];

                    disp(['Deg (0 ~ ' num2str(Ldata_IB) ') value of GAD for ' ...
                        datestr(thedates_D(index_IB)) ' used for correcting the inverted barometer (IB) effect.']);
                end

            end

            %         %m Add the degree 1 of GAD for ocean usage
            %         % Find a degree 1 data point that is within the month of our GRACE data
            %         where2=deg1dates_oce(:)>monthstart & deg1dates_oce(:)<monthend;
            %         if ~any(where2)
            %             % If there is no Deg1 value within our specific interval,
            %             % don't change anything, because we know a few months are missing
            %             disp('No change to degree 1')
            %         elseif nnz(where2) > 1
            %             % We have more than one month. Use the closest value
            %             disp('More than one degree 1 value for GAD detected')
            %             [~,where2]=min(abs(monthmid - deg1dates));
            %             lmcosi_month(2:3,3:4)=lmcosi_month(2:3,3:4)+squeeze(mydeg1_oce(where2,:,3:4));
            %             disp(['Deg1 value of GAD for ' datestr(deg1dates_oce(where2)) ' used.']);
            %         else
            %             disp(['Deg1 value of GAD for ' datestr(deg1dates_oce(where2)) ' used.']);
            %             lmcosi_month(2:3,3:4)=lmcosi_month(2:3,3:4)+squeeze(mydeg1_oce(where2,:,3:4));
            %         end

        end

        % Convert the geopotential coefficients into surface mass density,
        % if so desired
        if strcmp(units,'SD')
            % Need to make geoid first
            lmcosi_extra=plm2pot([lmcosi_month(:,1:2) lmcosi_month(:,5:6)*a],[],[],[],4);
            lmcosi_month=plm2pot([lmcosi_month(:,1:2) lmcosi_month(:,3:4)*a],[],[],[],4);
            % Add the fornal errors back in columns 5,6
            lmcosi_month=[lmcosi_month lmcosi_extra(:,3:4)];
        end

        % Combine into one matrix
        potcoffs(index,:,:) = lmcosi_month;

        %%%
        % CALIBRATED ERRORS
        %%%
        % We have no calibrated errors for CSR release 05, so we have to bypass
        % this section in this case.
        %     if (strcmp(Pcenters,'CSR') || strcmp(Pcenters,'JPL')) && (strcmp(Rlevel,'RL05')...
        %             || strcmp(Rlevel,'RL06'))
        %         cal_errors(index,:,:) = [lmcosi_month(:,1:2) zeros(size(lmcosi_month(:,1:2)))];
        %     else
        %         fname2=fullfile(ddir1,errornames{index});
        %
        %         % Open and scan the file (data from both Pcenters is 5 columns)
        %         fid = fopen(fname2);
        %         E = textscan(fid,'%s%s%s%s%s');
        %         fclose(fid);
        %
        %         % Only grab the lines for CALSDV
        %         Earray = cat(3,E{:});
        %         I = strmatch('CALSDV',Earray(:,1,1),'exact');
        %         Earray = squeeze(Earray(I,1,:));
        %
        %         % Only want columns 2-5, and as format double
        %         Earray = Earray(:,2:5);
        %         cal_errors_month=cellfun(@str2num,Earray);
        %         [m,n] = size(cal_errors_month);
        %
        %         % Change the order of the coefficients so that
        %         % order m goes as [0 01 012 0123 ...]
        %         revdel=[0 Ldata:-1:0];
        %         i=1;
        %         if strcmp(Pcenters,'CSR')
        %             new_ordering = zeros(m,4);
        %             demm=dems;
        %             for j=1:length(dems)
        %                 k = dels(i)+1 + sum( revdel( (1:dems(i) + 1 ) ) );
        %                 new_ordering(j,:) = cal_errors_month(k,:);
        %                 demm(j)=cal_errors_month(k,2);
        %                 i=i+1;
        %             end
        %             cal_errors_month = new_ordering;
        %         elseif strcmp('GSM-2_2006121-2006151_0028_EIGEN_G---_0004.txt',...
        %                 errornames(index)) || strcmp(Rlevel,'RL05')
        %             % for one very odd GFZ file
        %             new_ordering = zeros(m,4);
        %             i=4;
        %             for j=4:length(dems)
        %                 k = dels(i)+1 + sum( revdel( (1:dems(i) + 1 ) ) );
        %                 % This file has only 1 space for the 2,1 coefficients
        %                 if dems(i)==0
        %                     k=k-2;
        %                 else
        %                     k=k-3;
        %                 end
        %                 new_ordering(j-3,:) = cal_errors_month(k,:);
        %                 i=i+1;
        %             end
        %             cal_errors_month = new_ordering;
        %
        %         else % for the rest of GFZ, which has slightly less odd formatting
        %             new_ordering = zeros(m-1,4);
        %             i=4;
        %             for j=4:length(dems)
        %                 k = dels(i)+1 + sum( revdel( (1:dems(i) + 1 ) ) );
        %                 % These files have two spaces for the 2,1 coefficients
        %                 if j == 5
        %                     k=k-3;
        %                 else
        %                     k=k-2;
        %                 end
        %                 new_ordering(j-3,:) = cal_errors_month(k,:);
        %                 i=i+1;
        %             end
        %             cal_errors_month = new_ordering;
        %         end
        %
        %         % If from the GFZ data center, add terms for l=0 and 1
        %         if Pcenters == 'GFZ'
        %             cal_errors_month = [0 0 0 0; 1 0 0 0; 1 1 0 0; cal_errors_month];
        %         end
        %
        %         % Replace the C20 error from GRACE with the C20 error from SLR since we
        %         % used the C20 coefficient from SLR
        %         disp(sprintf('C20 error was %12.8e now %12.8e',cal_errors_month(4,3),slrc20_error(where)))
        %         cal_errors_month(4,3)=slrc20_error(where);
        %
        %         % Replace the Deg1 error from GRACE with the Deg1 error from Swenson et al.
        %         if ~any(where1)
        %             % Do nothing here
        %         else
        %             cal_errors_month(2:3,3:4)=squeeze(mydeg1(where1,:,5:6));
        %         end
        %
        %         % Convert the geopotential error coefficients into surface mass
        %         % density, if so desired
        %         if strcmp(units,'SD')
        %             % Need to make geoid first
        %             a=fralmanac('a_EGM96','Earth');
        %             cal_errors_month=plm2pot([cal_errors_month(:,1:2) cal_errors_month(:,3:4)*a],[],[],[],4);
        %         end
        %
        %         % Combine into one matrix
        %         cal_errors(index,:,:) = cal_errors_month;
        %
        %     end % We have no errors?

        disp(['Completed month number: ' num2str(index)]);
        disp(' ')
    end

    % calculate the IB correction for specific products
    % actually it is necessary for GRACE-FO (as no TN-12 is provided)
    % but we use GAA for GRGS, GAD for other institutions
    if strcmp(Pcenters,'GRGS')
        % we need to calculate the IB correction
        % The GRGS coefficients are in the order we want already;
        % GAA products do NOT need to add the zero and one coefficients
        %             lmcosi_month_D = [0 0 0 0 0 0; 1 0 0 0 0 0; 1 1 0 0 0 0; lmcosi_month_D];

        if  strcmp(Pcenters,'GRGS')
            if strcmp(Rlevel,'RL05')
                %GRGS only has GAA*_0090_0005
                datanames_IB=ls2cell(fullfile(ddir1,'GAA*0090_0005.txt'));
                Ldata_IB=90;
            elseif strcmp(Rlevel,'RL06')
                error('No RL06 for now.')
            end
        end

        for index_IB = 1:length(datanames_IB)
            % load geopotential coefficients of GAD
            fname4=fullfile(ddir1,datanames_IB{index_IB});

            % Calculate the midpoint of this data span
            monthstart_IB = datenum([str2num(datanames_IB{index_IB}(7:10))...
                1 str2num(datanames_IB{index_IB}(11:13))]);
            monthend_IB = datenum([str2num(datanames_IB{index_IB}(15:18))...
                1 str2num(datanames_IB{index_IB}(19:21))]);
            monthmid_IB = (monthstart_IB+monthend_IB)/2;
            thedates_IB(index_IB) = monthmid_IB;

            % Open and scan the file (data from all three centers is 10 columns)
            fid = fopen(fname4);
            C = textscan(fid,'%s%s%s%s%s%s%s%s%s%s');
            fclose(fid);

            % Only grab the lines for GRCOF2
            Carray = cat(3,C{:});
            I = strmatch('GRCOF2',Carray(:,1,1),'exact');
            Carray = squeeze(Carray(I,1,:));

            % Only want columns 2-7, and as format double
            Carray = Carray(:,2:7);
            lmcosi_month_IB=cellfun(@str2num,Carray);
            % This should be addmup(Ldata)
            [m,n] = size(lmcosi_month_IB);

            ldim=length(lmcosi_month_IB);

            % Convert the geopotential coefficients into surface mass density
            % (kg/m^2) 1 kg/m^2 /1000 kg/m^3 *1000 = 1 mm
            % Need to make geoid first
            lmcosi_extra_IB=plm2pot([lmcosi_month_IB(:,1:2) lmcosi_month_IB(:,5:6)*a],[],[],[],4);
            lmcosi_month_IB=plm2pot([lmcosi_month_IB(:,1:2) lmcosi_month_IB(:,3:4)*a],[],[],[],4);
            % Add the fornal errors back in columns 5,6
            lmcosi_month_IB=[lmcosi_month_IB lmcosi_extra_IB(:,3:4)];

            [IB,lon,lat,Plm,degres]=plm2xyz(lmcosi_month_IB,1,c11cmn2,Ldata_IB);
            %                 GAA(land_in)=NaN;
            IB_C(index_IB,:)=[index_IB sum(sum(IB(~land_in)))/sum(sum(~land_in))];

            disp(['Deg (0 ~ ' num2str(Ldata_IB) ') value of GAA for ' ...
                datestr(thedates_D(index_D)) ' used for correcting the inverted barometer (IB) effect.']);
        end

    end

    if strcmp(land_or_ocean,'ocean')
        save(fnpl_IB,'IB_C','thedates');
    end

    % correct GIA effect
    if ~strcmp(GIA,'No')

        if strcmp(units,'POT')
            error('Sorry, the GIA removal process does not support direct correct of geopotential (POT) yet. Please try surface density (i.e., "SD") instead.')
        end

        [thedates,GIAt]=correct4gia(thedates,GIA);

        for index = 1:nmonths
            potcoffs(index,:,3:4) = potcoffs(index,:,3:4) -  GIAt(index,1:size(potcoffs,2),3:4) ;
        end

    end

    % Save
    save(fnpl,'potcoffs','thedates','leakage_sp');

end % End if we have a save file already

% Collect output
% Here we have "thedates" twice so that we don't break older code. But in
% the future we will fix this so that we only have the two output
varns={potcoffs,thedates,thedates};
% Luckily, we can give additional information about the time epoch of
% leakage data by using this bug under specific cases
if strcmp(Pcenters,'GRGS')
    varns={potcoffs,leakage_sp,thedates};
end
varargout=varns(1:nargout);
