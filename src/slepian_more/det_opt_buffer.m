function buffer_deg=det_opt_buffer(basicInfo, slepianInfo);

arguments
    basicInfo   struct
    slepianInfo struct
end

% ---------------------------------------------------------------
% 0. Unpack basic fields
% ---------------------------------------------------------------
Code_Version = basicInfo.Code_Version;
TH_ori       = basicInfo.TH_ori;
Area         = basicInfo.Area;
Lwindow      = basicInfo.Lwindow;
c11cmn       = basicInfo.c11cmn;
Max_S        = basicInfo.Max_S;
land_or_ocean= basicInfo.land_or_ocean;
ddir1        = basicInfo.ddir1;
ddir2        = basicInfo.ddir2;
ifilesRoot   = basicInfo.ifilesRoot;

redo         = basicInfo.redo;
plotProcess  = basicInfo.plotProcess;
saveAddData  = basicInfo.saveAddData;

Radius            = slepianInfo.Radius;
group_buffer      = slepianInfo.group_buffer;

if ~redo && ~plotProcess && exist(fullfile(ddir1,'MainSlep_Optimal_Buffer.mat'), 'file')

    disp(['Slepian: Loading existing ',fullfile(ddir1,'MainSlep_Optimal_Buffer.mat') ...
        ' for plotting.']);
    load(fullfile(ddir1,'MainSlep_Optimal_Buffer.mat'),"buffer_deg");

    return

end

if strcmp(land_or_ocean,'land')
    % checkerboard

    ii=2; % This should have only 1 value.
    oo=0; % This could have more than 1 value.
    rr=3; % This should have no more than 1 value.

elseif strcmp(land_or_ocean,'ocean')
    ii=2; % This should have only 1 value.
    oo=43; % This could have more than 1 value.
    rr=[]; % This should have no more than 1 value.

elseif strcmp(land_or_ocean,'ice')
    ii=7; % This should have only 1 value.
    oo=0; % This could have more than 1 value.
    rr=[]; % This should have no more than 1 value.

end

MOD_S_selection = struct('mode', 'eigen_threshold', ...
    'value', [], 'threshold', 0.1, 'radius', []);

if isfield(basicInfo, 'sharedResultDir')
    ddir3 = fullfile(basicInfo.sharedResultDir, 'Simulate');
else
    ddir3 = fullfile(ddir1, 'Simulate');
end
if exist(ddir3, 'dir') ~= 7
    mkdir(ddir3);
end

group_buffer_num=numel(group_buffer);
[Dataproduct, FIG_Attach_Mod] = MOD_build_name_tags_and_ib_path( ...
    ii, oo, rr, Area, Lwindow);

simulationFile = fullfile(ddir3, ['Sig_' FIG_Attach_Mod '.mat']);
simulationMetadata = struct( ...
    'schema_version', 'V3-simulation-1', ...
    'TH_ori', TH_ori, ...
    'Area', Area, ...
    'land_or_ocean', land_or_ocean, ...
    'Lwindow', Lwindow, ...
    'Max_S', Max_S, ...
    'group_buffer', group_buffer, ...
    'c11cmn', c11cmn, ...
    'Radius', Radius, ...
    'complete', true);
reuseSimulation = false;
if ~redo && exist(simulationFile, 'file')
    cachedSimulation = load(simulationFile, ...
        'group_buffer_results', 'simulationMetadata');
    reuseSimulation = isfield(cachedSimulation, 'group_buffer_results') && ...
        isfield(cachedSimulation, 'simulationMetadata') && ...
        isequaln(cachedSimulation.simulationMetadata, simulationMetadata);
end

if reuseSimulation
    disp(['Slepian: Loading compatible shared simulation ', ...
        simulationFile, ' for plotting.']);
    group_buffer_results = cachedSimulation.group_buffer_results;
else
    if ~redo && exist(simulationFile, 'file')
        warning('V3:SharedSimulationMismatch', ...
            ['The shared simulation cache is incompatible with the current ', ...
             'region/groupBuffer settings and will be recomputed.']);
    end

    disp('This process may take a while. Please wait and enjoy a cup of coffee.')

    group_buffer_results=struct();
    % redo=true;
    for ll=1:group_buffer_num % normal one

        ll_buffer=group_buffer(ll);
        disp(['-------- Try buffer zone: ' num2str(ll_buffer) char(176) ...
            ' --------'])

        if ll_buffer>=0
            buffer_str=num2str(ll_buffer);
            buffer_str(buffer_str=='.')='p';
        else
            buffer_str=['neg' num2str(abs(ll_buffer))];
            buffer_str(buffer_str=='.')='p';
        end

        %% Projects model maps into the requested Slepian basis
        % 2. Decide on your choice of basis, depending on your region of interest
        % and the bandwidth you want. Using GRACE2SLEPT, project the results of
        % GRACE2PLMT into your chosen basis. We recommend that your basis is
        % chosen based on a set of synthetic experiments which estimate the
        % leakage/recovery tradeoffs.

        [XY_buffer, BasinArea_buffer, XY_ori, BasinArea_ori, Earth_radius] = ...
            compute_regions(TH_ori, ll_buffer);

        maybe_plot_region(XY_buffer, XY_ori, FIG_Attach_Mod, ll_buffer, ...
            ddir2, plotProcess);

        % mass density

        % because we don't need missing months for simulation, so it is simple
        [mod_slepcoffs,mod_dates,Coeff,TH_buffer,~,CC,V]=TWSt2slept_model(Dataproduct,...
            TH_ori,ll_buffer,Lwindow,0);

        % Example grid from first Slepian function
        [r_example, ~, ~, ~, ~] = plm2xyz(CC{1}, 1, c11cmn, Lwindow); %#ok<ASGLU>


        mod_nmonths=numel(mod_dates);
        mod_months=1:numel(mod_dates);

        [S_shannon, S, S_sig] = choose_S_from_eigenvalues(MOD_S_selection, ...
            Max_S, V, Lwindow, XY_buffer);

        mod_MASS_SS=zeros(S,length(mod_dates));
        mod_EWH_SS=zeros(S,length(mod_dates));
        %% Slepian residuals and integration
        % 3. Next run SLEPT2RESID to fit a choice of functions (e.g. lines,
        % quadratics, etc.) to the Slepian coefficients. Note: if you want to
        % remove a model of GIA then you should do that before this step, using
        % a function like CORRECT4GIA.

        fitwhat=[1 365.25 ]; %quadratic periodic residual
        %     fitwhat=[3]; %quadratic periodic residual

        [TH_use, XY, BasinArea] = MOD_select_integration_region( ...
            land_or_ocean, TH_ori, TH_buffer, XY_ori, XY_buffer, ...
            BasinArea_buffer, BasinArea_ori);

        % 4. If you want to examine the total mass trend in a region, this
        % information is calculated and returned as a result from SLEPT2RESID.
        % To summarize, each Slepian coefficient up to the Shannon truncation is
        % multiplied by the integral (found using INTEGRATEBASIS) of the corresponding
        % function over the region. This total data is then fit with TIMESERIESFIT
        % which allows fitting with a custom data variance.

        % smooth or not?
        if numel(S_sig)>1 && Max_S>S_sig(1)
            [mod_signalcoffs,mod_residcoffs,~,~,~,~,~,...
                ~,~,functionintegrals,~]=slept2resid_m(mod_slepcoffs,...
                mod_dates,fitwhat,[],[],CC,TH_use,[S_sig(1),Max_S],Radius);

            r_record = zeros(S2, size(r_example,1), size(r_example,2));
            for j = 1:S1
                [r_record(j,:,:), ~, ~, ~, ~] = ...
                    plm2xyz(CC{j}, 1, c11cmn, Lwindow);
            end
            for j = S1+1:S2
                CC_smo = plm_Gausssmooth(CC{j}, Radius);
                [r_record(j,:,:), ~, ~, ~, ~] = ...
                    plm2xyz(CC_smo, 1, c11cmn, Lwindow);
            end

        else
            [mod_signalcoffs,mod_residcoffs,~,~,~,~,~,...
                ~,~,functionintegrals,~]=slept2resid_m(mod_slepcoffs,...
                mod_dates,fitwhat,[],[],CC,TH_use,Max_S);

            r_record = zeros(S, size(r_example,1), size(r_example,2));
            for j = 1:S
                [r_record(j,:,:), ~, ~, ~, ~] = ...
                    plm2xyz(CC{j}, 1, c11cmn, Lwindow);
            end

        end

        Table2=zeros(Max_S,3);Table2_un=zeros(Max_S,4);
        mod_Grid_EWH_rms=zeros(Max_S, size(r_example,1), size(r_example,2));
        for SS=2:Max_S
            %% Individual Slepian residual (temporal)

            [Cab] = slepresid2cov(mod_residcoffs);

            % Make the coefficients with reference to some mean
            % If they already are, then this won't matter
            mod_sleptdelta = mod_slepcoffs(1:mod_nmonths,:) - repmat(mean(mod_slepcoffs(1:mod_nmonths,:),1),mod_nmonths,1);
            % no missing values necessary for model, so this is simple
            mod_signaldelta=mod_signalcoffs(1:mod_nmonths,:) - repmat(mean(mod_signalcoffs(1:mod_nmonths,:),1),mod_nmonths,1);
            mod_residdelta=mod_residcoffs(1:mod_nmonths,:) - repmat(mean(mod_residcoffs(1:mod_nmonths,:),1),mod_nmonths,1);

            eigfunINT=functionintegrals(1:SS);
            % Here do the total sum of the data
            %             total=eigfunINT*sleptdelta(:,1:N_all)';

            % Get the error
            %             thevars = diag(Cab(1:N_all,1:N_all))';
            %             alphavar = eigfunINT.^2.*thevars;
            % Now the combined error with covariance
            alphavarall = eigfunINT*Cab(1:SS,1:SS)*eigfunINT';

            [mod_MASS_CC,~,~,~,~,~]=integral_fit(SS, mod_dates, ...
                functionintegrals,mod_sleptdelta,Cab);
            mod_MASS=sum(mod_MASS_CC);

            % notion that here Cab do have actual meaning (we only want the sum-up).
            % So we do not use the error calculation here.
            [mod_MASSsig_CC,~,~,...
                ~,~,~]=integral_fit(SS,mod_dates,functionintegrals,mod_signaldelta,Cab);
            mod_MASSsig=sum(mod_MASSsig_CC);

            % notion that here Cab do have actual meaning (we only want the sum-up).
            % So we do not use the error calculation here.
            [mod_MASSres_CC,~,~,...
                ~,~,~]=integral_fit(SS,mod_dates,functionintegrals,mod_residdelta,Cab);
            mod_MASSres=sum(mod_MASSres_CC);

            MOD_MASSreconst=mod_MASSsig+mod_MASSres;
            mod_MASS_SS(SS,:)=mod_MASS;

            % change mass to EWH (transform Gt to kg/m^2 to cm: Gt * (10^9 * 10^3 / BasinArea) * (/10^3 *10^2))
            mod_EWH=mod_MASS*10^9*10^3/BasinArea/10;
            mod_EWH_SS(SS,:)=mod_EWH;

            %% Temporal statistic
            [~,~,~,~,Amplitude,Phase]=periodic_analysis_m(mod_EWH',...
                mod_dates,fitwhat); % sin(wt-fi)
            Phase=pi/2-Phase; % be consistent with previous studies for cos(wt+fi)
            Phase(Phase<0)=Phase(Phase<0)+2*pi; % -180~-180 to 0~360

            x0=[1,1,1,1];
            g=fittype(['a+b*x+a1*sin(2*pi*x/',num2str(fitwhat(2)),'+fi1)']...
                , 'coeff',{'a','b','a1','fi1'});
            [curve,goodness]=fit(mod_dates',mod_EWH',g,'StartPoint', x0);
            mod_interval=confint(curve);

            [~,~,mod_EWH_params,mod_EWH_paramserrors] = timeseriesfit([mod_dates' mod_EWH'],...
                alphavarall,1,1);
            %         mod_Mass_paramserrors = mod_Mass_paramserrors*365.25; % from Gt/day to Gt/year

            %         Total_EST_fit=EST_totalfit*10^9*10^3/BasinArea/10;
            mod_EWH_slope=mod_EWH_params(2)*365.25; % from /day to /year
            mod_EWH_signal_Amp=Amplitude;
            mod_EWH_signal_Pha=Phase;

            mod_EWH_paramerrors=mod_EWH_paramserrors(2)*365.25; % from /day to /year

            % Table PLOT for time series of different mass sea level
            Table2(SS,1)=mod_EWH_slope;
            Table2(SS,2)=mod_EWH_signal_Amp(1);
            Table2(SS,3)=mod_EWH_signal_Pha(1)/pi*180;

            Table2_un(SS,1)=mod_EWH_paramerrors;
            Table2_un(SS,2:4)=(mod_interval(2,2:end)-mod_interval(1,2:end))/2; % fit with periodic terms

            %% Spatial statistic

            mod_Grid_EWH = zeros(mod_nmonths, size(r_example,1), size(r_example,2));
            for i = 1:mod_nmonths

                sp_ewh   = 0;

                for j = 1:S
                    r = squeeze(r_record(j,:,:));
                    sp_ewh   = sp_ewh   + r * mod_sleptdelta(i,j)' / 1000 * 100;

                end
                mod_Grid_EWH(i,:,:) = sp_ewh;

            end
            mod_Grid_EWH_rms(SS,:,:)=squeeze(rms(mod_Grid_EWH,1));

        end

        Table2_un(:,[4])=Table2_un(:,[4])/pi*180; % phase by fiting with periodic term

        group_buffer_results(ll).mod_dates  = mod_dates;
        group_buffer_results(ll).Table2     = Table2;
        group_buffer_results(ll).Table2_un  = Table2_un;
        group_buffer_results(ll).S_sig      = S_sig;
        group_buffer_results(ll).mod_EWH_SS = mod_EWH_SS;
        group_buffer_results(ll).mod_Grid_EWH_rms = mod_Grid_EWH_rms;
        group_buffer_results(ll).Coeff      = Coeff;

        simulationMetadata.complete = false;
        save(simulationFile, 'group_buffer_results', 'simulationMetadata')
        simulationMetadata.complete = true;

    end
    save(simulationFile, 'group_buffer_results', 'simulationMetadata')

end
%% load useful plotting variables

run(fullfile(stpc_source_root(), 'MOD', 'Coeff_make.m')) % Coeffin and Coeffout

mod_dates       = group_buffer_results(1).mod_dates;
mod_nmonths     = numel(mod_dates);
mod_months      = 1:mod_nmonths;
Coeff           = group_buffer_results(1).Coeff;

if any(rr)
    if oo==0
        ratio=0.5;
    else
        error(['It is rather difficult to calculate the true trend, ' ...
            'amplitude or phase if the checkboard has different values.'])
    end
else
    ratio=1;
end

Coeff_pp=[Coeffin(ii)];

% in
in_trend=Coeff(1).Coef_trend;
in_A=Coeff(1).Coef_A;in_w=Coeff(1).Coef_w;in_fi=Coeff(1).Coef_fi;
in_Eps=Coeff(1).Coef_Eps;

% we will not simulate the noise here
in_ts=in_trend*mod_dates+in_A*sin(in_w*mod_dates+in_fi); % should be sin

in_ts=in_ts-mean(in_ts);

if oo~=0
    Coeff_pp=[Coeffin(ii), Coeffout(oo)];

    % out
    out_trend=Coeff(2).Coef_trend;
    out_A=Coeff(2).Coef_A;out_w=Coeff(2).Coef_w;out_fi=Coeff(2).Coef_fi;
    out_Eps=Coeff(2).Coef_Eps;

    % we will not simulate the noise here
    out_ts=out_trend*mod_dates+out_A*sin(out_w*mod_dates+out_fi); % should be sin

    out_ts=out_ts-mean(out_ts);

end

True_in_ts=in_ts*ratio;
True_trend=Coeff(1).Coef_trend*365.25/10*ratio;
True_amp=Coeff(1).Coef_A/10*ratio;
True_pha=(pi/2-Coeff(1).Coef_fi)/pi*180;


line_tre_n1=zeros();line_amp_n1=zeros();line_pha_n1=zeros();
buffer_pcc=zeros();buffer_rmse=zeros();
buffer_ts=zeros(mod_nmonths,group_buffer_num);
for ll=1:group_buffer_num

    line_tre_n1(ll)=group_buffer_results(ll).Table2(group_buffer_results(ll).S_sig,1);
    line_amp_n1(ll)=group_buffer_results(ll).Table2(group_buffer_results(ll).S_sig,2);
    line_pha_n1(ll)=group_buffer_results(ll).Table2(group_buffer_results(ll).S_sig,3);

    buffer_ts(:,ll)=group_buffer_results(ll).mod_EWH_SS(group_buffer_results(ll).S_sig,:)'*10;

    buffer_pcc(ll)=corr(True_in_ts',buffer_ts(:,ll));
    buffer_rmse(ll)=rmse(True_in_ts',buffer_ts(:,ll));
end

[Max_pcc,I_pcc] = max(buffer_pcc);
[Max_rmse,I_rmse] = min(buffer_rmse);
bufferDecision = buffer_selection_diagnostics(group_buffer, buffer_pcc, ...
    buffer_rmse, line_tre_n1, line_amp_n1, line_pha_n1, ...
    True_trend, True_amp, True_pha, I_pcc, I_rmse);

%% plot
lon_matl = c11cmn(1)-0.5 : c11cmn(3)-0.5;
lat_matl = c11cmn(2)+0.5 : -1 : c11cmn(4)+0.5;
[lonlon_matl,latlat_matl] = meshgrid(lon_matl,lat_matl);

mod_time=2005+mod_months/12-1/24;

XY_buf={};
for ll=1:group_buffer_num

    ll_buffer=group_buffer(ll);

    XY_buf{ll}=eval(sprintf('%s(%i,%f)',TH_ori,0,ll_buffer));

end

% plot the basis
F_Position=[2100,-250,1050,1200];
abc='abcdefghijklmnopqrstuvwxyz';
colorr_blind=[0,114,178;230,159,0;0,158,115;213,94,0]/255; % Blind-friendly colors
fontsize_n=12;line_width=2;

f2=figure;
set (gcf,'Position',F_Position)

F_gap1=[.07 .04];F_marg_h1=[.7 .05];F_marg_w1=[.06 .1];
h1=tight_subplot(1,3,F_gap1,F_marg_h1,F_marg_w1);

F_gap2=[.03 .05];F_marg_h2=[.35 .35];F_marg_w2=[.06 .04];
h2=tight_subplot(1,1,F_gap2,F_marg_h2,F_marg_w2);
xbou2=[2004.8,2016.2];

F_gap3=[.03 .06];F_marg_h3=[.05 .7];F_marg_w3=[.06 .04];
h3=tight_subplot(1,3,F_gap3,F_marg_h3,F_marg_w3);
xbou3=[min(group_buffer)-0.2,max(group_buffer)+0.2];

group_FaceAlpha=linspace(1,0,group_buffer_num+1);
if strcmp(land_or_ocean,'ocean')

    axes(h1(1));

    m_proj('mercator','long',[lon_matl(1), lon_matl(end)],'lat',[lat_matl(end), lat_matl(1)]);
    % m_coast('patch',[253 245 230]/255)
    m_coast('patch',[255 255 255]/255)
    m_grid('tickdir','in','yaxislocation','right',...
        'xaxislocation','bottom');

    XY_ori=eval(sprintf('%s(%i,%f)',TH_ori,0,0));
    XY_out=eval(sprintf('%s(%i,%f)',[TH_ori '_C'], 0, 1.5)); %1.5 is by default

    m_patch(XY_out(:,1), XY_out(:,2), [65 105 225]/255 )
    hold on
    m_patch(XY_ori(:,1), XY_ori(:,2), [255 255 255]/255 )
    hold on


    m_line(XY_ori(:,1), XY_ori(:,2),'LineStyle','-',...
        'Color','k','linewidth',3);
    hold on
    for ll=1:group_buffer_num

        p1=m_patch(XY_buf{ll}(:,1), XY_buf{ll}(:,2), [255 69 0]/255 );
        %     p1(ll).FaceAlpha=0.3;

        hold on

        % this is too empirical, so the buffer contour and study area my be
        % separated with 360?
        if XY_buf{ll}(1,1)>180
            m_line(XY_buf{ll}(:,1)-360, XY_buf{ll}(:,2),'LineStyle','--',...
                'Color',colorr_blind(ll,:),'linewidth',1);
        else
            m_line(XY_buf{ll}(:,1), XY_buf{ll}(:,2),'LineStyle','--',...
                'Color',colorr_blind(ll,:),'linewidth',1);
        end
        hold on

        p1.LineStyle='--';
        p1.FaceAlpha=group_FaceAlpha(ll);

    end

    % m_patch(XY_shp2(:,1), XY_shp2(:,2), [255 69 0]/255 )
    set(gca,'FontName','times new Roman','fontsize',fontsize_n)

    title({['(' abc(1) ') ' Area],'input/output-ocean simulation'},...
        'FontSize',fontsize_n+2,'FontWeight','bold');

elseif strcmp(land_or_ocean,'land')

    axes(h1(1));

    m_proj('mercator','long',[lon_matl(1), lon_matl(end)],'lat',[lat_matl(end), lat_matl(1)]);
    % m_coast('patch',[253 245 230]/255)
    m_coast('patch',[255 255 255]/255)
    m_grid('tickdir','in','yaxislocation','right',...
        'xaxislocation','bottom');

    XY_ori=eval(sprintf('%s(%i,%f)',TH_ori,0,0));
    poly_shp=polyshape(XY_ori(:,1)',XY_ori(:,2)');

    nn=0;pp=0;
    XY_c={};
    XY_ci={};
    for i=1:ceil(180/rr)
        for j=1:ceil(360/rr)
            XY_c{nn+1}=[2*rr*j-2*rr,90-2*rr*i+2*rr; 2*rr*j-rr,90-2*rr*i+2*rr; ...
                2*rr*j-rr,90-2*rr*i+rr; 2*rr*j-2*rr,90-2*rr*i+rr;];
            poly = polyshape(XY_c{nn+1}(:,1)',XY_c{nn+1}(:,2)');
            polyout_shp = intersect(poly,poly_shp);

            if ~isempty(polyout_shp.Vertices)
                XY_ci{pp+1}=polyout_shp.Vertices;
                pp=pp+1;
            end

            XY_c{nn+2}=[2*rr*j-rr,90-2*rr*i+rr; 2*rr*j,90-2*rr*i+rr; ...
                2*rr*j,90-2*rr*i; 2*rr*j-rr,90-2*rr*i;];
            poly2 = polyshape(XY_c{nn+2}(:,1)',XY_c{nn+2}(:,2)');
            polyout2_shp1 = intersect(poly2,poly_shp);

            if ~isempty(polyout2_shp1.Vertices)
                XY_ci{pp+1}=polyout2_shp1.Vertices;
                pp=pp+1;
            end

            nn=nn+2;
        end
    end

    LineStyle_ord=["--","--","--","-"];
    XY_ll_ci={};XY_buf={};
    for ll=1:group_buffer_num

        ll_buffer=group_buffer(ll);

        XY_buf{ll}=eval(sprintf('%s(%i,%f)',TH_ori,0,ll_buffer));

        poly_ll=polyshape(XY_buf{ll}(:,1)',XY_buf{ll}(:,2)');

        pp=0;
        for i=1:length(XY_c)

            polyi = polyshape(XY_c{i}(:,1)',XY_c{i}(:,2)');
            polyout_ll = intersect(polyi,poly_ll);

            if ~isempty(polyout_ll.Vertices)

                if sum(isnan(polyout_ll.Vertices(:,1)))>0

                    data=polyout_ll.Vertices(:,1);
                    % find the location of NaN
                    nanIdx = isnan(data);

                    % Use cumsum to segment-mark non-NaN portions
                    segmentIdx = cumsum(nanIdx);

                    % Remove NaN values while preserving segmentation
                    segments = {};
                    uniqueSegments = unique(segmentIdx(~nanIdx));
                    for j = 1:length(uniqueSegments)
                        segments{end+1} = polyout_ll.Vertices(segmentIdx == uniqueSegments(j) & ~nanIdx,:);
                        XY_ll_ci{ll,pp+1} = polyout_ll.Vertices(segmentIdx == uniqueSegments(j) & ~nanIdx,:);
                        pp=pp+1;
                    end

                    disp(segments);

                else

                    XY_ll_ci{ll,pp+1}=polyout_ll.Vertices;
                    pp=pp+1;

                end

            end

        end


        for i=1:pp
            p1(ll)=m_patch(XY_ll_ci{ll,i}(:,1), XY_ll_ci{ll,i}(:,2), [255 69 0]/255 );
            %     p1(ll).FaceAlpha=0.3;
            p1(ll).LineStyle=char(LineStyle_ord(ll));
            p1(ll).FaceAlpha=group_FaceAlpha(ll);
        end

        hold on

        % this is too empirical, so the buffer contour and study area my be
        % separated with 360?
        if XY_buf{ll}(1,1)>180
            l1(ll)=m_line(XY_buf{ll}(:,1)-360, XY_buf{ll}(:,2),'LineStyle','--',...
                'Color',colorr_blind(ll,:),'linewidth',1);
        else
            l1(ll)=m_line(XY_buf{ll}(:,1), XY_buf{ll}(:,2),'LineStyle','--',...
                'Color',colorr_blind(ll,:),'linewidth',1);
        end
        hold on

    end

    m_line(XY_ori(:,1), XY_ori(:,2),'LineStyle','-',...
        'Color','k','linewidth',3);

    set(gca,'FontName','times new Roman','fontsize',fontsize_n)

    title({['(' abc(1) ') ' Area],'checkerboard simulation'},'FontSize',14,'FontWeight','bold');

elseif strcmp(land_or_ocean,'ice')

    axes(h1(1));

    m_proj('stereographic','long',(c11cmn(1)+c11cmn(3))/2,'lat',(c11cmn(4)+c11cmn(2))/2,...
            'rad',21,'rec','on');
    % m_coast('patch',[253 245 230]/255)
    m_coast('patch',[255 255 255]/255)
    m_grid('tickdir','in','yaxislocation','right',...
        'xaxislocation','bottom');

    XY_ori=eval(sprintf('%s(%i,%f)',TH_ori,0,0));

    m_line(XY_ori(:,1), XY_ori(:,2),'LineStyle','-',...
        'Color','k','linewidth',3);
    hold on

    for ll=1:group_buffer_num

        p1=m_patch(XY_buf{ll}(:,1), XY_buf{ll}(:,2), [255 69 0]/255 );
        hold on
        %     p1(ll).FaceAlpha=0.3;
        % this is too empirical, so the buffer contour and study area my be
        % separated with 360?
        if XY_buf{ll}(1,1)>180
            m_line(XY_buf{ll}(:,1)-360, XY_buf{ll}(:,2),'LineStyle','--',...
                'Color',colorr_blind(ll,:),'linewidth',1);
        else
            m_line(XY_buf{ll}(:,1), XY_buf{ll}(:,2),'LineStyle','--',...
                'Color',colorr_blind(ll,:),'linewidth',1);
        end
        hold on

        p1.LineStyle='--';
        p1.FaceAlpha=group_FaceAlpha(ll);

    end

    % m_patch(XY_shp2(:,1), XY_shp2(:,2), [255 69 0]/255 )
    set(gca,'FontName','times new Roman','fontsize',fontsize_n)

    title({['(' abc(1) ') ' Area],'input/output-ocean simulation'},...
        'FontSize',fontsize_n+2,'FontWeight','bold');

end

% for plot (b)-(c)
axes(h1(2));
ll_buffer=1;

if strcmp(land_or_ocean,'ice')
    m_proj('stereographic','long',(c11cmn(1)+c11cmn(3))/2,'lat',(c11cmn(4)+c11cmn(2))/2,...
            'rad',21,'rec','on');
else
    m_proj('mercator','long',[lon_matl(1), lon_matl(end)],'lat',[lat_matl(end), lat_matl(1)]);
end

plot_map=squeeze(group_buffer_results(1).mod_Grid_EWH_rms(group_buffer_results(1).S_sig,:,:));
m_pcolor(lonlon_matl,latlat_matl,plot_map);

% m_coast('patch',[253 245 230]/255)
m_coast('patch',[255 255 255]/255)
m_grid('tickdir','in','yaxislocation','right',...
    'xaxislocation','bottom');

XY_buf_sub=XY_buf{ll_buffer};
shading flat;
colormap(mymap("coolwarm"));

hold on

m_line(XY_ori(:,1), XY_ori(:,2),'LineStyle','-',...
    'Color','k','linewidth',3);

hold on
% this is too empirical, so the buffer contour and study area my be
% separated with 360?
if XY_buf_sub(1,1)>180
    m_line(XY_buf_sub(:,1)-360, XY_buf_sub(:,2),'LineStyle','--',...
        'Color','k','linewidth',1);
else
    m_line(XY_buf_sub(:,1), XY_buf_sub(:,2),'LineStyle','--',...
        'Color','k','linewidth',1);
end

set(gca,'FontName','times new Roman','fontsize',fontsize_n)

title({['(' abc(2) ') ' Area],['buffer zone: ' num2str(group_buffer(ll_buffer)) ...
    '\circ']},'FontSize',fontsize_n+2,'FontWeight','bold');

axes(h1(3));
ll_buffer=group_buffer_num;

if strcmp(land_or_ocean,'ice')
    m_proj('stereographic','long',(c11cmn(1)+c11cmn(3))/2,'lat',(c11cmn(4)+c11cmn(2))/2,...
            'rad',21,'rec','on');
else
    m_proj('mercator','long',[lon_matl(1), lon_matl(end)],'lat',[lat_matl(end), lat_matl(1)]);
end

plot_map=squeeze(group_buffer_results(ll).mod_Grid_EWH_rms(group_buffer_results(ll).S_sig,:,:));
m_pcolor(lonlon_matl,latlat_matl,plot_map);

m_coast('patch',[255 255 255]/255)
m_grid('tickdir','in','yaxislocation','right',...
    'xaxislocation','bottom');

XY_buf_sub=XY_buf{ll_buffer};
shading flat;
colormap(mymap("coolwarm"));
cb = colorbar('h');

cb.Location = 'East';
cb.Position = [0.94 0.69 0.02 0.26];
title(cb,'RMS','fontsize',14);

% this is too empirical, so the buffer contour and study area my be
% separated with 360?
if XY_buf_sub(1,1)>180
    m_line(XY_buf_sub(:,1)-360, XY_buf_sub(:,2),'LineStyle','--',...
        'Color','k','linewidth',1);
else
    m_line(XY_buf_sub(:,1), XY_buf_sub(:,2),'LineStyle','--',...
        'Color','k','linewidth',1);
end

hold on

m_line(XY_ori(:,1), XY_ori(:,2),'LineStyle','-',...
    'Color','k','linewidth',3);

% m_patch(XY_shp2(:,1), XY_shp2(:,2), [255 69 0]/255 )
set(gca,'FontName','times new Roman','fontsize',fontsize_n)

title({['(' abc(3) ') ' Area],['buffer zone: ' ...
    num2str(group_buffer(ll_buffer)) '\circ']},...
    'FontSize',fontsize_n+2,'FontWeight','bold');

% for plot (d)-(g)

axes(h2(1));

li_n=0;

li_n=li_n+1;
p0_in=plot(mod_time,True_in_ts,'Linestyle','-', ...
    'LineWidth',line_width, 'Color','r');
le3_p=[p0_in];le3_text{li_n}='synthetic internal';
hold on

if oo~=0
    li_n=li_n+1;
    p0_out=plot(mod_time,out_ts,'Linestyle','-', ...
        'LineWidth',line_width, 'Color','b');
    hold on
    le3_p=[le3_p, p0_out];le3_text{li_n}='synthetic internal';
end

for ll=1:group_buffer_num
    li_n=li_n+1;
    p0_ts(ll)=plot(mod_time,buffer_ts(:,ll),'Linestyle','-', ...
        'LineWidth',line_width,'Color',colorr_blind(ll,:));
    hold on
    le3_p=[le3_p p0_ts(ll)];
    le3_text{li_n}=['inverse internal (buffer = ' num2str(group_buffer(ll)) '\circ)'];
end
xlim(xbou2)

legend(le3_p,le3_text,'NumColumns',3,'Location','northwest',...
    'FontSize',fontsize_n,'Box','off');
grid on

ylabel('cm','FontWeight','bold','fontsize',fontsize_n);

set(gca,'FontName','times new Roman','fontsize',fontsize_n);

title(['(' abc(4) ') synthetic and inverse internal signal time series'],...
    'FontSize',fontsize_n+2,'FontWeight','bold')

axes(h3(1));
p1_t=plot(group_buffer,line_tre_n1,'Linestyle','-', ...
    'LineWidth',line_width, 'Color','k');

hold on
p2_t=plot(group_buffer,ones(1,group_buffer_num)*True_trend,'Linestyle','-', ...
    'LineWidth',1, 'Color','r');
% text(mean(xbou3),True_trend,['b = ' num2str(True_trend)],...
%     'Color','r','fontsize',fontsize_n);
ylabel('cm/year','FontWeight','bold');
grid on
title({['(' abc(5) ') inverse trend']},'FontSize',fontsize_n+2,'FontWeight','bold');

xlabel('buffer zone range (degree)','FontWeight','bold')
l1=legend([p1_t,p2_t],{'inverse ($\lambda_S = 0.1$)','synthetic'},'Interpreter','latex','Box','off','fontsize',fontsize_n);

l1.Location='southeast';
l1.NumColumns=1;
% yticks(-0.11:0.01:-0.09)
xlim(xbou3)

set(gca,'FontName','times new Roman','fontsize',fontsize_n)

axes(h3(2));

p1_a=plot(group_buffer,line_amp_n1,'Linestyle','-', ...
    'LineWidth',line_width , 'Color','k');

hold on
p2_a=plot(group_buffer,ones(1,group_buffer_num)*True_amp,'Linestyle','-', ...
    'LineWidth',1, 'Color','r');
% text(mean(xbou3),True_amp,['A = ' num2str(True_amp)],'Color','r','fontsize',fontsize_n)
ylabel('cm','FontWeight','bold')
grid on

title({['(' abc(6) ') inverse annual amplitude']},'FontSize',fontsize_n+2,'FontWeight','bold')

xlabel('buffer zone range (degree)','FontWeight','bold')
l1=legend([p1_a,p2_a],{'inverse ($\lambda_S = 0.1$)','synthetic'},'Interpreter','latex','Box','off','fontsize',fontsize_n);

l1.Location='northwest';
l1.NumColumns=1;
xlim(xbou3)

set(gca,'FontName','times new Roman','fontsize',fontsize_n)

axes(h3(3));

p1_p=plot(group_buffer,line_pha_n1,'Linestyle','-', ...
    'LineWidth',line_width , 'Color','k');

hold on
p2_p=plot(group_buffer,ones(1,group_buffer_num)*True_pha,'Linestyle','-', ...
    'LineWidth',1, 'Color','r');
% text(mean(xbou3),True_pha,['\phi = ' num2str(True_pha)],'Color','r','fontsize',fontsize_n)
ylabel('degree','FontWeight','bold')
grid on

title({['(' abc(7) ') inverse annual phase']},'FontSize',fontsize_n+2,'FontWeight','bold')

xlabel('buffer zone range (degree)','FontWeight','bold')
l1=legend([p1_p,p2_p],{'inverse ($\lambda_S = 0.1$)','synthetic'},'Interpreter','latex','Box','off','fontsize',fontsize_n);

l1.Location='northeast';
l1.NumColumns=1;
xlim(xbou3)

set(gca,'FontName','times new Roman','fontsize',fontsize_n)

% Save figure
tif_name = ['Optimal_buffer_zone.tif'];
print(f2, '-dtiff', '-r500', fullfile(ddir2, tif_name));


%% determine the final suitable buffer zone

if redo || ~exist(fullfile(ddir1,'MainSlep_Optimal_Buffer.mat'), 'file')

    if I_pcc==I_rmse
        buffer_deg = group_buffer(I_rmse);
        buffer_selection_source = 'pcc_rmse_agreement';
        disp(['PCC and RMSE agree. The optimal buffer zone is: ', ...
            num2str(buffer_deg) char(176), '. Continue automatically.']);

    else
        disp('PCC and RMSE select different buffer zones.');
        disp(['PCC-best buffer: ', num2str(group_buffer(I_pcc)), ...
            char(176), '; RMSE-best buffer: ', ...
            num2str(group_buffer(I_rmse)), char(176), '.']);
        disp('Buffer comparison (estimated value, true value, and absolute error):');
        disp(bufferDecision.table);
        disp(bufferDecision.reason);

        defaultChoice = num2str(bufferDecision.auto_index);
        promptText = sprintf([ ...
            'Select any candidate number shown in the table above.\n' ...
            'If no valid input is received within 30 seconds, candidate %s ' ...
            '(buffer %.6g deg) will be selected automatically.\n' ...
            'Candidate [%s]: '], defaultChoice, ...
            bufferDecision.auto_buffer_deg, defaultChoice);
        [candidateAnswer, usedDefault] = timed_buffer_input( ...
            promptText, defaultChoice, 30);
        selectedIndex = str2double(candidateAnswer);

        if isfinite(selectedIndex) && selectedIndex == fix(selectedIndex) && ...
                ismember(selectedIndex, 1:numel(group_buffer))
            buffer_deg = group_buffer(selectedIndex);
            if usedDefault
                buffer_selection_source = 'secondary_metrics_timeout_auto';
            else
                buffer_selection_source = 'user_candidate_selection';
            end
        else
            error(['Invalid buffer candidate "%s". Run again and enter an ' ...
                'integer candidate number from 1 to %d.'], ...
                candidateAnswer, numel(group_buffer));
        end
    end

    close all

    buffer_selection_diagnostics_table = bufferDecision.table;
    buffer_auto_choice_deg = bufferDecision.auto_buffer_deg;
    buffer_selection_reason = bufferDecision.reason;
    save(fullfile(ddir1,'MainSlep_Optimal_Buffer.mat'), ...
        'buffer_deg', 'buffer_pcc', 'buffer_rmse', ...
        'buffer_selection_diagnostics_table', 'buffer_auto_choice_deg', ...
        'buffer_selection_reason', 'buffer_selection_source');

else

    disp(['Slepian: Loading existing ',fullfile(ddir1,'MainSlep_Optimal_Buffer.mat') ...
        ' for plotting.']);
   
    load(fullfile(ddir1,'MainSlep_Optimal_Buffer.mat'),"buffer_deg");

    disp(['Your selected buffer zone: ' num2str(buffer_deg) char(176) ])
    
    pause(1)
    
    close all
end
%%
end

function [answer, usedDefault] = timed_buffer_input( ...
        promptText, defaultAnswer, timeoutSeconds)
    answer = '';
    usedDefault = false;
    fprintf('%s', promptText);
    fprintf('Waiting %d seconds; default is "%s".\n', timeoutSeconds, defaultAnswer);

    if usejava('jvm')
        try
            inputReader = javaObject('java.io.InputStreamReader', java.lang.System.in);
            bufferedReader = javaObject('java.io.BufferedReader', inputReader);
            waitTimer = tic;
            while toc(waitTimer) < timeoutSeconds
                if bufferedReader.ready()
                    answer = strtrim(char(bufferedReader.readLine()));
                    break
                end
                pause(0.2);
            end
        catch
            pause(timeoutSeconds);
        end
    else
        pause(timeoutSeconds);
    end

    if isempty(answer)
        answer = defaultAnswer;
        usedDefault = true;
        fprintf('No input received in %d seconds. Use default "%s".\n', ...
            timeoutSeconds, defaultAnswer);
    end
end
