
% Prepare run
clear
clc
close all
rng(0,'twister');




% Global variables
global n_subpops useOnlyS useFixedLambda useFixedS fixedS fixedLambda ...
    fixedEpsilon useFixedEpsilon useLambdaEst ...
    maxNumberOfFractions useImmuneRT useRelImmune n_weeks_btwFrac ...
    useLambdadie usegamma0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start OF INPUT SECTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Path
workfolder  = pwd;

% Enter filename with data (each patient holds an individual sheet)
filenameXLS = ;


% Patients: 
% Included patients
usepat = [1,2];

% Maximum number of patients in study
nPat   = 2; 


% Immune effects
useImmuneRT    = false;
useRelImmune   = false;


% Bootstrap
varyDataPoints = true;
n_vary         = 50;


% Surviving fraction (S) settings and choice of number of subpopulations
n_subpops = 1;
useOnlyS  = true;
useFixedS = false;
fixedS    = 0.3;

% Parameter settings
useLambdadie = false;
usegamma0    = false;

% Growth rate settings
useFixedLambda    = true;
useLambdaEst      = false;
fixedLambda       = 0.0027; % In h^-1
useLambdaPreForV0 = false;

if useFixedLambda
% Used for growth rate grid search 
% allLambdas = [linspace(log(2)./(40*24),log(2)./(15*24),10),...
%               linspace(log(2)./(16*24),log(2)./(8*24),19),...
%               linspace(log(2)./(7*24),log(2)./(5*24),3)];

allLambdas = [0.0027];

else
    allLambdas = [1];
end


% Resistence rate
useFixedEpsilon = false;
fixedEpsilon    = 1.9663e-04;



% Option to exclude pateints without pretreatment data 
excludePatsWithNOPretreatmentData = false;


% How to start the run
useProcessedRun            = true;
useProcessedRunsForStart   = false;
loadParametersFromthisPath = 'my_filename';

% Plot?
doplots = true;
    

% Maximum number of fractions: 
n_maxFrac = [11]; %[5,7,9,11,13] 

% Time between fractionation
allfractions = [6];%[4,6,8,10,12] 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF INPUT SECTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop over different times between fractions
for timeBtwfrac = allfractions
n_weeks_btwFrac = timeBtwfrac;

    % Loop over all growth rates
    for il=1:length(allLambdas)
        disp(['At step ' num2str(il) 'of ' num2str(length(allLambdas))])

        % Loop over maximum number of fractions
        for mf = n_maxFrac

            % Experiment name
            maxNumberOfFractions = mf;
            expName_all = ['Bootstrap_' num2str(n_weeks_btwFrac) '_' num2str(maxNumberOfFractions)];
            if usegamma0
                expName_all = [expName_all '_gamma0free'];
            else
                expName_all = [expName_all '_gamm0eqLambda'];
            end
            uncertData  = cell(length(usepat),length(allLambdas),n_vary);

            if useFixedLambda
                fixedLambda = allLambdas(il);    
                expName_new = [expName_all '_lambda' num2str(fixedLambda)];
            else
                expName_new = [expName_all '_Fitlambda'];
            end

            if il ==1
                expName = expName_new;
                mkdir(expName)
            end

            % Folder/data containers to save to 
            allData = cell(nPat,1);

            % Loop over all patients
            for i=usepat
                
                disp(['Working on patient ', num2str(i)])
                if useProcessedRunsForStart || useProcessedRun

                    % Load reference data 
                    curdir = pwd;
                    cd(loadParametersFromthisPath)
                    load('refParameters.mat')
                    load('uncertData.mat')
                    cd(curdir)
                end



                % Load data (tumor volumes and time points)
                [num,txt,raw] = xlsread(filenameXLS,i);
                n_vols  = sum(~isnan(num(:,7)));
                volData = NaN(n_vols,2);
                inds    =  find(~isnan(num(:,7)));
                for j=1:n_vols
                    volData(j,1) = num(inds(j),1);
                    volData(j,2) = num(inds(j),7);
                end    
                allData{i}=volData;


                % Pretreatment data
                preTrVolData = volData(volData(:,1)<0,:);

                % Exceptions of patients with surgery - make sure to only use
                % post surgical volumes
                if any(diff(preTrVolData(:,2))<0)
                   neg_inds     = find(diff(preTrVolData(:,2))<0);
                   theind       = neg_inds(end);
                   preTrVolData = preTrVolData(theind+1:end,:);
                else
                    theind = 0;
                end


                % Get pre/post treatment volumes and times 
                indstart      = find(volData(:,1)>=0,1);
                postTrVolData = volData(indstart:end,:);
                data_use_orig = [preTrVolData;postTrVolData];
                time          = data_use_orig(:,1);


                % Loop over number of bootstraps (if any)
                if varyDataPoints
                    close all
                    for iv = 1:n_vary

                        % Show progress
                        disp(num2str(iv))

                        % First run with original data
                        if iv>1

                            if useProcessedRun
                                data_use      = uncertData{i,il,iv};
                            else

                                % Original
                                data_use  =  data_use_orig;


                                % Vary each data point
                                for d = 1: size(data_use,1)

                                    % Different strategy for volumina > or <
                                    % 2cmˆ3: 
                                    if data_use(d,2)>=2.

                                        % std of 0.2 mean 1 for 'large' volumes
                                        std_var  = 0.2;
                                        mean_var = 1;
                                        rand_var = std_var.*randn(size(data_use(d,2))) + mean_var;

                                        % Check last element to prevent
                                        if d == size(data_use,1)
                                            while rand_var.*data_use(d,2)<=data_use(d-1,2)
                                                rand_var = std_var.*randn(size(data_use(d,2))) + mean_var;
                                            end
                                        end

                                        % Assign
                                        data_use(d,2) = rand_var.*data_use(d,2);
                                        if data_use(d,2)<0
                                            data_use(d,2) = 0;
                                        end

                                    else
                                        % Volumina <= 2 cm^3
                                        if data_use(d,2)>0
                                            std_var = 0.5./(data_use(d,2)+0.5);
                                            mean_var = 1;
                                            rand_var = std_var.*randn(size(data_use(d,2))) + mean_var;
                                            if d == size(data_use,1)
                                                while rand_var.*data_use(d,2)<=data_use(d-1,2)
                                                    rand_var = std_var.*randn(size(data_use(d,2))) + mean_var;
                                                end
                                            end

                                            data_use(d,2) = rand_var.*data_use(d,2);
                                        else
                                            data_use(d,2) = 0;
                                        end

                                        % Ensure positive or zero values 
                                        % (negative volumia do not make sense)
                                        if data_use(d,2)<0
                                            data_use(d,2) = 0;
                                        end
                                    end                                         
                                end
                                use_vol = data_use;
                                uncertData{i,il,iv} = data_use;
                            end
                        else
                            data_use = data_use_orig;
                        end



                    % Fit the RT model to the post treatment data or load a
                    % previous run
                    use_vol  = data_use(:,2);
                    time     = data_use(:,1);
                    rt_on    = 0;
                    if useProcessedRun
                        curdir = pwd;
                        load([loadParametersFromthisPath '/refParameters.mat'])

                        if n_vary>1
                            load([loadParametersFromthisPath '/uncertData.mat'])
                        end
                        pars_out = refParameters{i,il,iv}(2*n_subpops+2:2*n_subpops+2+4);
                        init_out = refParameters{i,il,iv}(1:2*n_subpops+1);
                        cd(curdir)

                    else
                        pars_out = NaN;
                        init_out = NaN;
                    end  
                    xdata = time;
                    ydata = use_vol;



                    % Fit parameters and predict alternative schedules
                    [Init,params,rs,...
                        t_vendSPART(i,iv),...
                        t_6weeks(i,iv),...
                        v_end_standard(i,iv),...
                        t_vend_standard(i,iv), v_standard, t_standard, v_PART,t_PART,...
                        v_6weeks, t_6week,v_data,t_data,... 
                        gamma_standard, gamma_PART, gamma_6weeks,aic, r2] = ...
                        Opt_PSO_pat_RT_Bev_Pem_withLambda_fmin_use(xdata,ydata,i,rt_on,...
                         expName,pars_out, init_out);

                    % Save results                            
                    v0_RT = v_standard(t_standard==0);
                    allDataout_standard{i,il,iv}.t_standard = t_standard;
                    allDataout_standard{i,il,iv}.v_standard = v_standard;
                    allDataout_standard{i,il,iv}.t_PART     = t_PART;
                    allDataout_standard{i,il,iv}.v_PART     = v_PART;
                    allDataout_standard{i,il,iv}.t_6week    = t_6week;
                    allDataout_standard{i,il,iv}.v_6week    = v_6weeks;
                    allDataout_standard{i,il,iv}.v_data     = v_data;
                    allDataout_standard{i,il,iv}.t_data     = t_data;
                    allDataout_standard{i,il,iv}.rs         = rs;
                    allDataout_standard{i,il,iv}.v0         = v0_RT;
                    allDataout_standard{i,il,iv}.time       = time*24;
                    allDataout_standard{i,il,iv}.post_vol   = use_vol;
                    allDataout_standard{i,il,iv}.gamma_standard= gamma_standard;
                    allDataout_standard{i,il,iv}.gamma_PART = gamma_PART;
                    allDataout_standard{i,il,iv}.gamma_6week= gamma_6weeks;
                    allDataout_standard{i,il,iv}.aic        = aic;
                    allDataout_standard{i,il,iv}.r2         = r2;

                    allParams{i,il,iv} = [Init,params,rs,NaN,NaN,...
                            t_vend_standard(i),t_6weeks(i),t_vendSPART(i),v_end_standard(i)];
                    t_standard(i) = time(end)*24;
     
                    end


                    % Final plots: Tumor volume and gamma
                    if doplots
                        plot_allGrowth(allDataout_standard,i,il,expName,uncertData)
                    end
                end

            end
            refParameters = allParams;


            % Save everything  
            cd(workfolder)
            currFolder = pwd;
            cd(expName)
            save('allParams_New')
            cd(currFolder)


            try

                % Prepare KM plots
                t_6weeks(t_vend_standard==0)     = NaN;
                t_vendSPART(t_vend_standard==0)  = NaN;
                t_6weeks_l(t_vend_standard==0)   = NaN;
                t_vendSPART_l(t_vend_standard==0)= NaN;
                t_6weeks_u(t_vend_standard==0)   = NaN;
                t_vendSPART_u(t_vend_standard==0)= NaN;
                t_vend_standard_l(t_vend_standard==0)= NaN;
                t_vend_standard_u(t_vend_standard==0)= NaN;
                t_vend_standard(t_vend_standard==0)  = NaN;

                if varyDataPoints
                    all_t_pro = [t_vend_standard(:,1)';t_6weeks(:,1)';t_vendSPART(:,1)'];
                    all_t_pro_l = [nanmin(t_vend_standard');nanmin(t_6weeks');nanmin(t_vendSPART')];
                    all_t_pro_u = [nanmax(t_vend_standard');nanmax(t_6weeks');nanmax(t_vendSPART')];
                    all_t_pro_m = [nanmean(t_vend_standard');nanmean(t_6weeks');nanmean(t_vendSPART')];

                else
                    all_t_pro = [t_vend_standard;t_6weeks;t_vendSPART];
                    all_t_pro_l = [t_vend_standard_l;t_6weeks_l;t_vendSPART_l];
                    all_t_pro_u = [t_vend_standard_u;t_6weeks_u;t_vendSPART_u];
                    all_t_pro_m = all_t_pro;
                end

                % Do KM plot with bootstrap - for HFSRT, iRT+boost, iRT
                plot_multiple_KaplanMeier_new(all_t_pro_m',all_t_pro_u',all_t_pro_l')
                cd(expName)
                    filename = '_KM_uncertFig.fig';
                    savefig(gcf,filename)
                    filename = 'KM_uncert.png';
                    saveas(gcf,filename)
                    cd(currFolder)


                % KM no bootstrap but events - HFSRT, iRT,
                % iRT+boost
                timeVar  = [t_vend_standard,t_6weeks,t_vendSPART];
                groupVar = [1*ones(size(t_vend_standard)),2*ones(size(t_6weeks)),3*ones(size(t_vendSPART))];
                eventVar = ones(size(timeVar));
                timeVar  = [t_vend_standard(usepat)/24,t_6weeks(usepat)/24, t_vendSPART(usepat)/24];
                groupVar = cell(size(timeVar));
                groupVar(1:length(t_vend_standard(usepat))) = {'HFSRT'};
                groupVar(length(t_vend_standard(usepat))+1:...
                    length(t_vend_standard(usepat))+length(t_6weeks(usepat))) = {'iRT'};
                groupVar(length(t_vend_standard(usepat))+...
                    length(t_6weeks(usepat))+1:end) = {'iRT+Boost'};

                eventVar = ones(size(timeVar));
                [p, fh, stats] = MatSurv(timeVar, eventVar, groupVar);
                    cd(expName)
                    filename = 'KM_all.fig';
                    savefig(fh,filename)
                    filename = 'KM_all.png';
                    saveas(fh,filename)
                    cd(currFolder)


                % KM plot with log rank test but no bootstrap intervals - HFSRT,iRT+boost    
                timeVar = [t_vend_standard(usepat)/24, t_vendSPART(usepat)/24];
                groupVar = cell(size(timeVar));
                groupVar(1:length(t_vend_standard(usepat))) = {'HSFR'};
                groupVar(length(t_vend_standard(usepat))+1:end) = {'iRT+Boost'};
                eventVar = ones(size(timeVar));
                [p, fh, stats] = MatSurv(timeVar, eventVar, groupVar);
                    cd(expName)
                    filename = 'KM_56_PART.fig';
                    savefig(fh,filename)
                    filename = 'KM_56_PART.png';
                    saveas(fh,filename)
                    cd(currFolder)
                cd(expName)
                save('allParams_New')
                cd(currFolder)


                % Evaluate fit
                markers = ['+','o','*','x','s','d','^','v','>','<','p','h','+','o','*','x',...
                          's','d','^','v','>','<','p','h'];
                counter = 1;
                figure1 = figure;
                hold all
                for p = usepat
                    preV  = interp1(allDataout_standard{p,1,1}.t_standard, ...
                           allDataout_standard{p,1,1}.v_standard,...
                           allDataout_standard{p,1,1}.time/24);
                    if p == 1
                        allV = log(preV);
                        allVdata =log(allDataout_standard{p,1,1}.post_vol);
                    else
                        allV = [allV; log(preV)];
                        allVdata =[allVdata;log(allDataout_standard{p,1,1}.post_vol)];
                    end
                    legenduse{counter} = num2str(p);
                    h=plot(log(preV),log(allDataout_standard{p,1,1}.post_vol),'o',...
                        'MarkerSize', 5,'Marker', markers(counter));
                    try
                    set(h,'MarkerFaceColor', get(h,'Color')); 
                    catch
                    end
                    counter = counter+1;  
                end


                % Calculate the slope of v_measured vs. v_fitted (log)
                datax = allV;
                datay = allVdata;
                x = datax(~isinf(datay));
                y = datay(~isinf(datay));

                [fitres,qualityFit] = fit(x,y,'poly1');
                xplot = linspace(min(x),max(x),100);
                yplot = fitres.p1*xplot+fitres.p2;
                plot(xplot,yplot,'k-','Linewidth', 2)

                box on
                xlabel('log(measured volume)')
                ylabel('log(fitted volume)')
                set(gca,'Fontsize',16)
                annotation(figure1,'textbox',...
                    [0.144214285714286 0.830476190476191 0.1815 0.0833333333333348],...
                    'String',['R^2 = ',num2str(qualityFit.rsquare)],...
                    'FontSize',14,...
                    'FitBoxToText','off',...
                    'EdgeColor','none');

                cd(expName)
                filename = 'vsimvsVdata.fig';
                savefig(gcf,filename)
                filename = 'vsimvsVdata.png';
                saveas(gcf,filename)
                cd(currFolder)

                catch
                    disp('plots failed')

            end

        end
    end
end








