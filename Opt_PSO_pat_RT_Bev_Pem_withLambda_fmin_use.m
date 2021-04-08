function [Init,parnew,rmse,t_vendSPART,t_6weeks,v_end_standard,...
t_vend_standard,v_standard, t_standard, v_PART,t_PART,...
v_6weeks, t_6week,v_data,t_data, gamma_standard, gamma_PART,...
gamma_6weeks,aic,r2] = Opt_PSO_pat_RT_Bev_Pem_withLambda_fmin_use(xdata,...
ydata, patient,rt_on, expName,pars_out, init_out)


% Initialize all parameters
global time_firstRT n_subpops useFixedLambda ...
     fixedLambda tpre data useLambdadie usegamma0

if useFixedLambda
    lambda = fixedLambda;
end
v_data     = NaN; t_data     = NaN; 
v_standard = NaN; t_standard = NaN; 
v_PART     = NaN; t_PART     = NaN;
v_6weeks   = NaN; t_6week    = NaN;
t_vendSPART= NaN; t_6weeks   = NaN;



% Create structure with data and initial conditions (will be input to
% optimizer and solution function)
data.ydata    = ydata;
data.xdata    = xdata;
tpre          = xdata(1);
time_firstRT  = rt_on*24; 

if ~isnan(pars_out)
    % No optimization done (in case of loaded results)
    Init = init_out;   %take initial conditions
    parnew = pars_out; %take out parameters
    aic = NaN; 
else

    % Get nominal parameter values, upper and lower bound for optimization:
    %startP - start point
    %LB - parameter lower bounds
    %UB - parameter upper bounds 
    [startP,LB,UB] = getParametersAnalytical(ydata(1),patient,true); 

    tic
    
    % Fmincon optimization - custom costfunction
    A   = [];
    b   = [];
    Aeq = [];
    beq = []; 
    options = optimoptions(@fmincon,'Algorithm','sqp','MaxIterations',500,...
              'Display','iter','StepTolerance',1e-9, ...
              'OptimalityTolerance',1e-9,...
              'MaxFunctionEvaluations',5000);
    fun = @(pars_in)costfun_RT_Bev_Pem(pars_in,data);
    [curve,~,~,~] = fmincon(fun,startP,A,b,Aeq,beq,...
                    LB, UB,[],options);
    epsilon = curve(1);
    lambda  = curve(3);
    S       = curve(4);
    vpre    = curve(5);
            
    % Check options
    if useLambdadie
        lambda_die = curve(6);
    else
        lambda_die = curve(3);
    end
    if usegamma0
        gamma0 = curve(2);
    else
        gamma0 = curve(3);
    end
    opt_pars     = [0, vpre, gamma0, S, 0.2,epsilon, lambda, lambda_die];            
         
    toc
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate AIC for this patient and model %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Number of observations
    n_obs    = length(ydata);
    
    
    % Parameters: V0, epsilon, S, Lambda, gamma0
    n_params = 5; 
    if useFixedLambda
    n_params = n_params-1;
    end
    if ~usegamma0
        n_params = n_params-1;
    end
    
    % Get RSS
    ypred = analyticSol_RT_Bev_Pem_eval_withPre(data.xdata*24, epsilon,gamma0,lambda,S,vpre);
    rss = sum((ydata-ypred).^2);
    
    % Get AIC 
    if n_params<n_obs/40
        aic = n_obs*log(rss/n_obs)+2*n_params;
    else
        aic = n_obs*log(rss/n_obs)+2*n_params+2*n_params*n_obs/(n_obs-n_params-1);
    end
    if isinf(aic)
       disp('Is inf!') 
    end
    

    %Solve the system using new parameters
    Init   = opt_pars(1:2*n_subpops+1); %take initial conditions
    parnew = opt_pars(2*n_subpops+2:end); %take out parameters

end


% Evaluate alternative treatment schedules (iRT (here called "6weeks")
% and iRT+boost (here called "simplePART"))
xplot = linspace(0,500*24);
if max(data.xdata)>500*24
    xplot = linspace(0,2*max(data.xdata)*24);
end
[v_pred,yout,tout,v_end_standard,t_vend_standard,gamma_standard] = solvediscontODE_RT_Bev_Pem_eval_withPre(parnew,Init,data.xdata*24, data.ydata);
[v_6weeks,yout_6weeks,tout_6weeks, t_6weeks,~,gamma_6weeks]      = solvediscontODE_6GY6weeks_withPre(parnew,Init,data.xdata*24,v_end_standard);
[v_sPART,yout_sPART,tout_sPART, t_vendSPART,gamma_PART]          = solvediscontODE_simplePART_Vmin_withPre(parnew,Init,data.xdata*24,v_end_standard);



% Sum up contribtuions
if n_subpops>1
    v_viab  = sum(yout(n_subpops+1:2*n_subpops,:));
    v_dying = sum(yout(1:n_subpops,:));
    v       = v_viab+v_dying;
    
    v_viab_6weeks  = sum(yout_6weeks(n_subpops+1:2*n_subpops,:));
    v_dying_6weeks = sum(yout_6weeks(1:n_subpops,:));
    v_6weeks       = v_viab_6weeks+v_dying_6weeks;
    
    v_viab_sPART  = sum(yout_sPART(n_subpops+1:2*n_subpops,:));
    v_dying_sPART = sum(yout_sPART(1:n_subpops,:));
    v_sPART       = v_viab_sPART+v_dying_sPART;
else
    v_viab  = yout(n_subpops+1,:);
    v_dying = yout(1,:);
    v       = v_viab+v_dying;
    
    v_viab_6weeks  = yout_6weeks(2,:);
    v_dying_6weeks = yout_6weeks(1,:);
    v_6weeks       = v_viab_6weeks+v_dying_6weeks;
    
    v_viab_sPART  = yout_sPART(2,:);
    v_dying_sPART = yout_sPART(1,:);
    v_sPART       = v_viab_sPART+v_dying_sPART;
end


% Evaluate goodness of fit
r2   = rsquare(ydata,v_pred);
rmse = sqrt(mean((ydata-v_pred).^2));


% Output overview
v_standard = v; 
t_standard = tout/24; 
v_PART     = v_sPART;
t_PART     = tout_sPART/24;
t_6week    = tout_6weeks/24;
v_data     = ydata;
t_data     = xdata;

end