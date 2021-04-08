function [ypred,yout,tout,v_end,t_vend,gammaout] = solvediscontODE_RT_Bev_Pem_eval_withPre(pars,Init,xdata,ydata)
global time_firstRT n_subpops useOnlyS

t_final_plot = 900*24; % Total time frame simulated
rt_on        = [0,1,2,3,4]*24+time_firstRT; % in h - time of RT delivery


%odeset('Events',@eventsPat,'OutputFcn',@odeplot,'OutputSel',1,...
%'Refine',refine);
options = [];   
   
   
% surviving fraction
if useOnlyS
    S = pars(1:n_subpops);
else
    error('Option useOnlyS=false is not supported!')
end


% Pretreatment
if any(xdata<0)
    lambda = pars(4*n_subpops:5*n_subpops-1);
    tout = [min(xdata)-14:1:0];
    yout = zeros(3,length(tout));
    yout(n_subpops+1:2*n_subpops,:) = Init(n_subpops+1:2*n_subpops)*exp(-(tout(1)-tout)*lambda);
    Init(n_subpops+1:2*n_subpops)=yout(n_subpops+1:2*n_subpops,end);  
else
    tout = [0];
    yout = [Init'];
end


% First fraction at time 0
v_viab       = Init(n_subpops+1:2*n_subpops);
v_dying      = Init(1:n_subpops);
v_dying      = v_dying+(1-S).*v_viab;
v_viab       = v_viab-(1-S).*v_viab;
gamma        = Init(2*n_subpops+1);
Init         = [v_dying,v_viab,gamma];
rt_on        = rt_on(2:end);
tstart       = 0;
tfinal       = rt_on(1);
numfractions = length(rt_on);
   
   

for i = 1:numfractions+1

   [t,y] = ode45(@model_ODE_pat_RT_Bev_Pem,[tstart:1:tfinal],Init,options,pars);
   y     = transpose(y);
   t     = transpose(t);
    
   % Accumulate output and set the new initial conditions
   v_viab = y(n_subpops+1:2*n_subpops,end)';
   v_dying =y(1:n_subpops,end)';
   gamma = y(1+2*n_subpops,end)';
   v_dying = v_dying+(1-S).*v_viab;
   v_viab = v_viab-(1-S).*v_viab;  
   Init=[v_dying,v_viab,gamma];
   tout = [tout t(2:end)];
   yout = [yout y(:,2:end)];
   
   % Reassign start/end times
   if i<numfractions
       tstart = rt_on(i);
       tfinal = rt_on(i+1);
   else
       tstart = tfinal;
       tfinal = t_final_plot; 
   end
     
end
yall = sum(yout(1:2*n_subpops,:),1);


% Simulate more time if needed to show the full range until time to
% regrowth
countex = 0;
while max(yall)<1.1*yall(1) && countex <10
    v_viab = y(n_subpops+1:2*n_subpops,end)';
    v_dying =y(1:n_subpops,end)';
    gamma = y(1+2*n_subpops,end)';
    Init=[v_dying,v_viab,gamma];
    tstart = tout(end);
    tfinal = tstart+6*6*7*24;
    sol = ode45(@model_ODE_pat_RT_Bev_Pem,[tstart:1:tfinal],Init,options,pars);
    [t,y] = ode45(@model_ODE_pat_RT_Bev_Pem,[tstart:1:tfinal],Init,options,pars);
    y = transpose(y);
    t = transpose(t);
    
   % Accumulate output and set the new initial conditions
   v_viab = y(n_subpops+1:2*n_subpops,end)';
   v_dying =y(1:n_subpops,end)';
   gamma = y(1+2*n_subpops,end)';
   Init=[v_dying,v_viab,gamma];
   tout = [tout t(2:end)];
   yout = [yout y(:,2:end)];
   yall = nansum(yout(1:2*n_subpops,:),1);
   countex = countex+1;
end



% Only use values which are not nan or inf
inds = ~isinf(yall)&~isnan(yall);
yall = yall(inds);
gammaout = yout(2*n_subpops+1,:); 
yout = yout(1:2*n_subpops,inds);
tout = tout(inds);
[uni_tout,uniqueinds]=unique(tout);

if length(uniqueinds)~=length(yall)
    yall=yall(uniqueinds);
    tout = uni_tout;
    yout = yout(1:2*n_subpops,uniqueinds);
    gammaout =gammaout(inds);
    gammaout =gammaout(uniqueinds);
end



% Interpolate at data points
try
    ypred = interp1(tout,yall,xdata);
catch
    disp(tout(inds));
end



% Get regrowth times - this has to be after the minimum vol has been
[~,minind] = min(yall);
v0_RT = yall(tout==0);
if ydata(end)>v0_RT
    try
        v_end = interp1(tout(minind:end),yall(minind:end),xdata(end));
    catch
        v_end = ydata(end);
    end
else
    v_end = 1.1*v0_RT;
end

try
    t_vend = interp1(yall(minind:end),tout(minind:end),v_end);
catch
    t_vend = 100000; 
end