function v = analyticSol_RT_Bev_Pem_eval_withPre(x, epsilon,gamma0,lambda,s,vpre)
global time_firstRT n_subpops tpre useLambdadie usegamma0


n_subpops = 1; 
useLambdadie  = false; 
usegamma0     = false; 
time_firstRT  = 0; 



x_in = x;
if useLambdadie 
    error("Patient 22 does not have enough data for this option")
else
    lambda_die = lambda;
end
if usegamma0
    gamma0 = gamma0;
else
    gamma0 = lambda;
end

t_final_plot = 1500*24;                     % Total time frame simulated
rt_on        = [0,1,2,3,4]*24+time_firstRT; % in h - time of RT delivery



% Pretreatment
tout = [tpre*24:1:0];
yout = zeros(3,length(tout));
yout(n_subpops+1:2*n_subpops,:) = vpre*exp(-(tout(1)-tout)*lambda);
Init(n_subpops+1:2*n_subpops)   = yout(n_subpops+1:2*n_subpops,end);  



% First fraction at time 0
v0_viab       = Init(n_subpops+1:2*n_subpops);
v0_dying      = Init(1:n_subpops);
v0_dying      = v0_dying+(1-s).*v0_viab;
v0_viab       = v0_viab-(1-s).*v0_viab;
gamma0_use    = gamma0;


% Set up fraction calculation
rt_on         = rt_on(2:end);
tstart        = 0;
tfinal        = rt_on(1)-1;
numfractions  = length(rt_on);


% Take out values at time 0 (otherwise double)
if length(tout)>0
tout = tout(1:length(tout)-1);
yout = yout(:,1:length(yout)-1);
end
   

for i = 1:numfractions+1
   
   % Analytical solution 
   x       = tstart:tfinal;
   v_viab  = v0_viab*exp(lambda*(x-tstart)+gamma0_use/epsilon*(exp(-epsilon*(x-tstart))-1));
   v_dying = v0_dying*exp(-lambda_die*(x-tstart));    
   gamma   = gamma0_use*exp(-epsilon*(x-tstart));
  
   
   % Accumulate output and 
   tout = [tout x];
   yout = [yout [v_viab;v_dying;gamma]];
   
   
   % Reassign start/end times
   if i<numfractions
       tstart = rt_on(i);
       tfinal = rt_on(i+1)-1;
   else
       tstart = tfinal+1;
       tfinal = t_final_plot; 
   end
   
   
   % Set the new initial conditions
   gamma0_use   = gamma0*exp(-epsilon*tstart);
   v0_dying = v_dying(end)+(1-s).*v_viab(end);
   v0_viab  = v_viab(end)-(1-s).*v_viab(end);  
     
end
yall = sum(yout(1:2*n_subpops,:),1);

% Interpolate at x
v = interp1(tout,yall,x_in);

