function [thedates,TWSlmcosi,Coeff,Checker_board]=MODEL_Checkerboard_Make_new_m(Dataproduct,TH,Lwindow,XY_buffer)
% Create the simulation model (Checkerboard Simulation)
% For this Checkerboard Simulation, the 'Dataproduct' consists of three
% aspects, ii and oo and rr. The ii and oo (normally oo=0) determines the 
% simulated signals and blanks in the study region. The rr denotes the
% grid intervals of the checkerboard (normally rr=3). 
% The checkerboard simulation is for regions where leak-in or leak-out
% effects are equally strong or weak. For example, river basins with few
% coastlines.

%% Coefficients and routes
defval('Dataproduct',{'2','0','5'});
defval('TH','SCSpTH');
defval('Lwindow',60);

num_in=str2num(Dataproduct{1});
num_out=str2num(Dataproduct{2});
Checker_reso=str2num(Dataproduct{3}); % resolution of the checkboard (in degree)

%     setenv('IFILES','/Users/mazhongtian/Academic/program/Slepian');
% setenv('IFILES','C:\Academic\program\Slepian');

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

fn1=fullfile(getenv('IFILES'), 'MOD');
fnmod=sprintf('%s/MODEL-I%sO%sR%s-%s-%i.mat',...
    fn1,Dataproduct{1},Dataproduct{2},Dataproduct{3},h,Lwindow);

if exist(fn1, 'dir') ~= 7
    mkdir(fn1);
end

Earth_radius = fralmanac('a_EGM96','Earth');

run(fullfile(stpc_source_root(),'MOD','Coeff_make.m'))

%% Make simulation
% c11cmn=[0 89.5 359 -89.5];
% c11cmn=[0 90 360 -90];
c11cmn=[0.5 89.5 359.5 -89.5];
lon_lim=[c11cmn(1) c11cmn(3)];
lat_lim=[c11cmn(2) c11cmn(4)];

lon=lon_lim(1):1:lon_lim(2);
lat=lat_lim(1):-1:lat_lim(2);
[lon_mesh, lat_mesh]=meshgrid(lon,lat);

Checker_board_block=zeros(Checker_reso*2,Checker_reso*2);
Checker_board_block(1:Checker_reso,1:Checker_reso)=1;
Checker_board_block(Checker_reso+1:Checker_reso*2,Checker_reso+1:Checker_reso*2)=1;

Checker_board=repmat(Checker_board_block,floor(size(lon_mesh,1)/Checker_reso/2)+1,floor(size(lon_mesh,2)/Checker_reso/2)+1);
Checker_board=Checker_board(1:size(lon_mesh,1),1:size(lon_mesh,2));

ll_bufer_in=[0,-0.5,-1,-1.5];
ll_bufer_out=[0,0.5,1,1.5];

if num_out==0
    % for land

    Coeff=[Coeffin(num_in)];

    Coeff(1).Lwindow=Lwindow;

    in_trend=Coeffin(num_in).Coef_trend;
    in_A=Coeffin(num_in).Coef_A;in_w=Coeffin(num_in).Coef_w;in_fi=Coeffin(num_in).Coef_fi;
    in_Eps=Coeffin(num_in).Coef_Eps;

    % How much will the inner study area expand? (unit: degree)
    % You should give the buffer zone.

    TH_in=TH_origin;

    % we actually keep all the mass in the Greenland
    %     XY_in=eval(sprintf('%s(%i,%f)',TH_in,0,ll_buffer_in(1)));
    XY_in=eval(sprintf('%s(%i,%f)',TH_in,0,XY_buffer));

    BasinArea_in=spharea(XY_in)*4*pi*Earth_radius^2;

    f1=figure
    subplot(121)
    plot(XY_in(:,1),XY_in(:,2),'k');

    [XY_in_in,XY_in_on] = inpolygon(lon_mesh,lat_mesh,XY_in(:,1),XY_in(:,2));
    % [XY_out_in,XY_out_on] = inpolygon(lon_mesh,lat_mesh,XY_out(:,1),XY_out(:,2));

    hold on
    plot(lon_mesh(XY_in_in & ~Checker_board),lat_mesh(XY_in_in & ~Checker_board) ,'bo');
    hold on
    plot(lon_mesh(XY_in_in & Checker_board),lat_mesh(XY_in_in & Checker_board),'ro');

    for i=2005:2015
        for j=1:12
            thedates((i-2005)*12+j)=datenum(i,j,15);
        end
    end

    in_ts=in_trend*thedates+in_A*sin(in_w*thedates+in_fi)+in_Eps*randn(1,numel(thedates)); % should be sin
    % in_ts=in_trend*thedates+in_A*cos(in_w*thedates+in_fi)+in_Eps*randn(1,numel(thedates));
    in_ts=in_ts-mean(in_ts);

    %     out_ts=out_trend*thedates+out_A*sin(out_w*thedates+out_fi)+out_Eps*randn(1,numel(thedates)); % should be sin
    % out_ts=out_trend*thedates+out_A*cos(out_w*thedates+out_fi)+out_Eps*randn(1,numel(thedates));
    %     out_ts=out_ts-mean(out_ts);
    out_ts=zeros(size(in_ts));

    subplot(122)
    plot(thedates,in_ts,'r');
    hold on
    plot(thedates,out_ts,'b')

    for i=1:numel(thedates)

        in_xyz=zeros(size(lat_mesh));
        out_xyz=zeros(size(lat_mesh));

        in_xyz(XY_in_in & Checker_board)=in_ts(i);
        out_xyz(XY_in_in & ~Checker_board)=out_ts(i);

        final_xyz=in_xyz+out_xyz;

        % Lwindow=60;

        % TWSlmcosi{i} = xyz2plm(final_xyz,Lwindow,[],lat_lim,lon_lim);
        TWSlmcosi{i} = xyz2plm(final_xyz,Lwindow);

    end

else

    Coeff=[Coeffin(num_in), Coeffout(num_out)];

    Coeff(1).Lwindow=Lwindow;Coeff(2).Lwindow=Lwindow;

    in_trend=Coeffin(num_in).Coef_trend;
    in_A=Coeffin(num_in).Coef_A;in_w=Coeffin(num_in).Coef_w;in_fi=Coeffin(num_in).Coef_fi;
    in_Eps=Coeffin(num_in).Coef_Eps;

    out_trend=Coeffout(num_out).Coef_trend;
    out_A=Coeffout(num_out).Coef_A;out_w=Coeffout(num_out).Coef_w;out_fi=Coeffout(num_out).Coef_fi;
    out_Eps=Coeffout(num_out).Coef_Eps;

    % How much will the inner study area shrink? (unit: degree)
    % default 0 (the original area of your study area)
    ll_in=1;
    % How much will the outer circle expand? (unit: degree)
    % default 1.5
    ll_out=4;

    TH_in=TH_origin;

    phi=0; theta=0; omega=0;

    XY_buffer_in=ll_bufer_in(ll_in);
    XY_buffer_out=ll_bufer_out(ll_out);

    XY_in=eval(sprintf('%s(%i,%f)',TH_in,0,XY_buffer_in));
    BasinArea_in=spharea(XY_in)*4*pi*Earth_radius^2;
    % XY_out=eval(sprintf('%s(%i,%f)',TH_out,0,XY_buffer_out));
    % BasinArea_out=spharea(XY_in)*4*pi*Earth_radius^2;

    f1=figure
    subplot(121)
    plot(XY_in(:,1),XY_in(:,2),'k');

    [XY_in_in,XY_in_on] = inpolygon(lon_mesh,lat_mesh,XY_in(:,1),XY_in(:,2));
    % [XY_out_in,XY_out_on] = inpolygon(lon_mesh,lat_mesh,XY_out(:,1),XY_out(:,2));

    hold on
    plot(lon_mesh(XY_in_in & ~Checker_board),lat_mesh(XY_in_in & ~Checker_board) ,'bo');
    hold on
    plot(lon_mesh(XY_in_in & Checker_board),lat_mesh(XY_in_in & Checker_board),'ro');

    for i=2005:2015
        for j=1:12
            thedates((i-2005)*12+j)=datenum(i,j,15);
        end
    end

    in_ts=in_trend*thedates+in_A*sin(in_w*thedates+in_fi)+in_Eps*randn(1,numel(thedates)); % should be sin
    % in_ts=in_trend*thedates+in_A*cos(in_w*thedates+in_fi)+in_Eps*randn(1,numel(thedates));
    in_ts=in_ts-mean(in_ts);

    out_ts=out_trend*thedates+out_A*sin(out_w*thedates+out_fi)+out_Eps*randn(1,numel(thedates)); % should be sin
    % out_ts=out_trend*thedates+out_A*cos(out_w*thedates+out_fi)+out_Eps*randn(1,numel(thedates));
    out_ts=out_ts-mean(out_ts);

    subplot(222)
    plot(thedates,in_ts,'r');
    hold on
    plot(thedates,out_ts,'b')

    for i=1:numel(thedates)

        in_xyz=zeros(size(lat_mesh));
        out_xyz=zeros(size(lat_mesh));

        in_xyz(XY_in_in & Checker_board)=in_ts(i);
        out_xyz(XY_in_in & ~Checker_board)=out_ts(i);

        final_xyz=in_xyz+out_xyz;

        % Lwindow=60;

        % TWSlmcosi{i} = xyz2plm(final_xyz,Lwindow,[],lat_lim,lon_lim);
        TWSlmcosi{i} = xyz2plm(final_xyz,Lwindow);

    end

end

% Save figure
tif_name = ['Model_Point_TS_Buf' num2str(XY_buffer) '.tif'];
print(f1, '-dtiff', '-r300', fullfile(fn1, tif_name));

% 1 for in, 0 for out
save(fnmod,'thedates','TWSlmcosi','Coeff','Checker_board');

end
