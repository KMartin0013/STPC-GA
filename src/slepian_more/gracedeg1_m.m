function varargout=gracedeg1_m(Pcenter,Rlevel,land_or_ocean)
% [thedates,deg1data]=GRACEDEG1(Pcenter,Rlevel,land_or_ocean)
%
% This function reads and formats the degree-1 spherical harmonic
% correction from GRACE/GRACE-FO Technical Note 13.
%
% INPUT:
%
% Pcenter     'CSR' data center at the Center for Space Research
%             'GFZ' data center at the GeoForschungsZentrum Potsdam
%             'JPL' data center at the Jet Propulsion Laboratory
% Rlevel      The release level of the solution you want.
%              Either 'RL04','RL05', or 'RL06'
% land_or_ocean
%             Whether the data was used on land or ocean. This will
%             change the process on degree 1
%
% OUTPUT:
%
% Returns these variables and saves them in a .mat file:
%    thedates       time stamps in Matlab time. This date is the midpoint
%                      of the GRACE month data span.
%    deg1data       potential coefficients for just degree l = 1
%                      [nmonths x 2 x 4]
%
% NOTES:
%
% The GRACE degree-1 correction data are distributed in Technical Note 13,
%   TN-13. There is one for each datacenter because the GRACE coefficients
%   from degree 2-60 are used in the construction of degree one terms. You
%   can find TN-13 in the docs folder in the GRACE data directory at NASA
%   PODAAC. https://podaac.jpl.nasa.gov/
%
% We check the file timestamps of both the TN-13 file and the Matlab .mat
%   file that we save/load. If the .mat file is newer, then we load it. If
%   the TN-13 document is newer, then we use it to create a new .mat file
%   and then load the .mat file.
%
% Last modified by zhongtian.ma@connect.polyu.hk, 10/15/2023
% Last modified by charig-at-arizona.edu, 11/18/2021

% Determine parameters and set defaults
defval('Pcenter','CSR')
defval('Rlevel','RL06')
defval('land_or_ocean','land')

% Where the original data files are kept
defval('ddir1',fullfile(getenv('IFILES'),'GRACE','Degree1'));
% Where you would like to save the new .mat file
defval('ddir2',fullfile(getenv('IFILES'),'GRACE','Degree1'));


% The name of the original data files

if strcmp(land_or_ocean,'land')

    if strcmp(Rlevel,'RL06')
        %used for IJAEOG
        %         fnpl1=sprintf('%s/TN-13_GEOC_%s_%s_FO.txt',ddir1,Pcenters,'RL06'); %m here changed from 'RL0601' and added 'FO' for longer records
        %         fnpl2=sprintf('%s/TN-13_GEOC_%s_%s_FO.mat',ddir2,Pcenters,'RL06'); %m here changed from 'RL0601' and added 'FO' for longer records

        % used for GRACE and GRACE-FO
        if strcmp(Pcenter,'CSR')
%             Release='RL0602';
            Release='RL0603'; % we update this
        elseif strcmp(Pcenter,'JPL') || strcmp(Pcenter,'GFZ')
%             Release='RL0601';
            Release='RL0603';  % we update this
        else
            Release='RLMean';
        end

        fnpl1=sprintf('%s/TN-13_GEOC_%s_%s.txt',ddir1,Pcenter,Release);
        fnpl2=sprintf('%s/TN-13_GEOC_%s_%s.mat',ddir2,Pcenter,Release);

    elseif strcmp(Rlevel,'RL05')
        fnpl1=sprintf('%s/deg1_RL05_MA.txt',ddir1); %m here changed
        fnpl2=sprintf('%s/deg1_RL05_MA.mat',ddir2); %m here changed
    elseif strcmp(Rlevel,'RL04')
        fnpl1=sprintf('%s/deg1_RL04_NH.txt',ddir1);
        fnpl2=sprintf('%s/deg1_RL04_NH.mat',ddir2);
    else
        error('GRACEDEG1: Wonky release level requested')
    end

elseif strcmp(land_or_ocean,'ocean')

    warning('Are you sure to use the degree 1 of GAD to substitude that from SLR? Press any key to continue.')
    pause
    % actually this is adundant to separate land and ocean according to the TN-13

    % Becuase the degree 1 values of GAD will be added to the degree 1 from SLR
    % afterwards in grace2plmt_m.m,
    % the degree 1 here should be from TN-13 rather than that of GAD.
    % The only possible usage of this separate maybe in study in river
    % plume, where GAD higher than degree 1 should be removed to
    % study the residual oceanic signal (the only difference to directly
    % use TWS is the additional degree 1 of GAD).
    % Because the paper of IJAEOG use this process, this process remained
    % here for further paper correction.

    if strcmp(Rlevel,'RL06')
        fnpl1=sprintf('%s/GAD_Degree1_GEOC_%s_%s.txt',ddir1,Pcenter,'RL06'); %m here changed from 'RL0601' and added 'FO' for longer records
        fnpl2=sprintf('%s/GAD_Degree1_GEOC_%s_%s.mat',ddir2,Pcenter,'RL06'); %m here changed from 'RL0601' and added 'FO' for longer records
    elseif strcmp(Rlevel,'RL05')
        fnpl1=sprintf('%s/GAD_Degree1_GEOC_%s_%s.txt',ddir1,Pcenter,'RL05'); %m here changed from 'RL0601' and added 'FO' for longer records
        fnpl2=sprintf('%s/GAD_Degree1_GEOC_%s_%s.mat',ddir2,Pcenter,'RL05'); %m here changed from 'RL0601' and added 'FO' for longer records
    elseif strcmp(Rlevel,'RL04')
        error('GRACEDEG1: Special Ocean Application for degree 1 is only valid for RL05 and RL06.');
        fnpl1=sprintf('%s/deg1_RL04_NH.txt',ddir1);
        fnpl2=sprintf('%s/deg1_RL04_NH.mat',ddir2);
    else
        error('GRACEDEG1: Wonky release level requested')
    end
else
    error('You must select either ocean or land.')

end

% Get the file date information
d1 = dir(fnpl1);
d2 = dir(fnpl2);



% If this file already exists, and it is more new than the source file,
% then load it.  Otherwise make a new one and return that one. This should
% automatically make and use a new .mat file if you update the source (TN-13) file.
if ~(exist(fnpl1,'file')==2)
    if strcmp(land_or_ocean,'ocean')
        disp('Please wait a bit for the creation of degree files for GAD.')
        extract_GAD_degree1(Pcenter,Rlevel)

    else

        warning(['No input degree 1 file for %s %s. We will use the averaged degree 1' ...
            ' files from CSR, JPL and GFZ instead.'],Pcenter, Rlevel)
    end
end
% Refresh file metadata because a marker file may have just been created.
d1 = dir(fnpl1);
d2 = dir(fnpl2);

% If this file already exists, and it is more new than the source file,
% then load it.  Otherwise make a new one and return that one. This should
% automatically make and use a new .mat file if you update the source (TN-13) file.
if ~(exist(fnpl2,'file')==2) || isempty(d1) || datenum(d1.date) > datenum(d2.date)
    % Here create a new save file to use.
    % This is reversed from our typical if, exist, load convention because
    % this way if the file does not exist it will skip the second if
    % condition (which would error on the datenum call) and go straight here

% % because we use ITSG data, the d1 file do not exist.
% if ~(exist(fnpl2,'file')==2)

    if ~strcmp(Pcenter,'CSR') && ~strcmp(Pcenter,'JPL') && ~strcmp(Pcenter,'GFZ')
        [thedates,deg1CSR]=gracedeg1_m('CSR',Rlevel);
        [thedates,deg1JPL]=gracedeg1_m('JPL',Rlevel);
        [thedates,deg1GFZ]=gracedeg1_m('GFZ',Rlevel);

        deg1filenum=min([size(deg1CSR,1),size(deg1JPL,1),size(deg1GFZ,1)]);
        if size(deg1CSR,1)~=deg1filenum || size(deg1JPL,1)~=deg1filenum || size(deg1GFZ,1)~=deg1filenum
            warning('The total data time spans for CSR, JPL and GFZ are different. We used the shortest time span.')
        end

        for i=1:min([size(deg1CSR,1),size(deg1JPL,1),size(deg1GFZ,1)])
            for j=1:2
                deg1Mean(i,j,1:2)=deg1CSR(i,j,1:2);
                for k=3:4
                    deg1Mean(i,j,k)=mean(deg1CSR(i,j,k)+deg1JPL(i,j,k)+deg1GFZ(i,j,k))/3;
                end
            end
        end

        thedates = thedates(1:min([size(deg1CSR,1),size(deg1JPL,1),size(deg1GFZ,1)]));
        deg1data = deg1Mean;
        
        fileID = fopen(fnpl1,'wt');
        fprintf(fileID,'Just to be consistent with other instutions and keep the process continue.');
        fclose(fileID);

    else

        if strcmp(Rlevel,'RL06')
            fid = fopen(fnpl1);
            tline = fgetl(fid);

            % Keep reading lines until we get to the 'end of header' line
            while isempty(tline) || ~strcmp(tline(1:5),'end o')
                tline = fgetl(fid);
            end

            % Now do a formatted text scan on the data
            C = textscan(fid,'%s%d%d%f%f%f%f%s%s%s');
            % Use string format for the dates so we can easily send it to datenum
            fclose(fid);

            % What we actually want is a 3 dimensional cell array with the date as
            % the first dimension and a lmcosi as dimension 2 and 3 (just like
            % grace2plmt)

            % Use the epoch start and end dates to get the midpoint date
            thedates = (datenum(C{8},'yyyymmdd') + datenum(C{9},'yyyymmdd'))./2;
            [b,m] = unique(thedates);
            thedates = thedates(m);

            % Re cast the ints as doubles to get it all in one matrix
            deg1data = [double(C{2}) double(C{3}) C{4} C{5}];

            % Now reorder into a lmcosi format for each month
            for i=1:(length(C{1}))/2
                temp = [deg1data(2*i-1,:);
                    deg1data(2*i,:)];
                deg1lmcosi(i,:,:) = temp;
            end

            deg1data = deg1lmcosi;

        end

        % m Add the usage of RL05
        if strcmp(Rlevel,'RL05')

            %              deg1=load(fnpl1);

            fileID = fopen(fnpl1,'r');
            formatSpec = '%6f%5f%5f%19s%19s%11s%11s%14f%f%[^\n\r]';
            C = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string',  'ReturnOnError', false);

            %              str=char(dataArray);
            %              data=char2num(str);

            deg1=[C{1} C{2} C{3} double(C{4}) double(C{5}) double(C{6}) double(C{7})];

            %              C = textscan(fid,'%s%d%d%f%f%f%f%s%s%s');

            [n,m] = size(deg1);
            dates_str = num2str(deg1(:,1));
            deg1dates = datenum([str2num(dates_str(:,1:4)) str2num(dates_str(:,5:6)) 15*ones(n,1)]);
            [b,m] = unique(deg1dates);
            thedates = deg1dates(m);
            for i=1:n/2
                temp = [deg1(2*i-1,2:7); deg1(2*i,2:7)];
                mydeg1(i,:,:) = temp;
            end

            deg1data = mydeg1;
        end

    end
    % Create a save file
    save(fnpl2,'thedates','deg1data')

else
    load(fnpl2)
    disp(sprintf('%s loaded by GRACEDEG1',fnpl2))
end





% Collect output
varns={thedates,deg1data};
varargout=varns(1:nargout);




