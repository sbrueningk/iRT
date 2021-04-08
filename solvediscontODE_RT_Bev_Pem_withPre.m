function [ypred,yout,tout] = solvediscontODE_RT_Bev_Pem_withPre(pars,Init,xdata,vend)
global time_firstRT n_subpops useOnlyS


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


rt_on = [0,1,2,3,4]*24+time_firstRT; % in h
dose = 6; % in Gy

options = [];%odeset('Events',@eventsPat,'OutputFcn',@odeplot,'OutputSel',1,...
   %'Refine',refine);

% surviving fraction
if useOnlyS
    S = pars(1:n_subpops);
else
    alpha = pars(1:n_subpops);
    beta = pars(n_subpops+1:2*n_subpops);
    S = exp(-alpha*dose-beta*dose^2);
end


% First fraction at time 0
v_viab = Init(n_subpops+1:2*n_subpops);
v_dying = Init(1:n_subpops);
v_dying = v_dying+(1-S).*v_viab;
v_viab = v_viab-(1-S).*v_viab;
gamma0 = pars(4);%Init(2*n_subpops+1);
Init=[v_dying,v_viab,gamma0];
rt_on =rt_on(2:end);
tstart = 0;
tfinal = rt_on(1);
numfractions = length(rt_on);
   
   

for i = 1:numfractions+1

   sol = ode45(@model_ODE_pat_RT_Bev_Pem,[tstart tfinal],Init,options,pars);
    
   % Accumulate output and set the new initial conditions
   v_viab = sol.y(n_subpops+1:2*n_subpops,end)';
   v_dying =sol.y(1:n_subpops,end)';
   gamma = sol.y(1+2*n_subpops,end)';
   v_dying = v_dying+(1-S).*v_viab;
   v_viab = v_viab-(1-S).*v_viab;  
   Init=[v_dying,v_viab,gamma];
   tout = [tout sol.x(2:end)];
   yout = [yout sol.y(:,2:end)];
   
  
   if i<numfractions
        tstart = rt_on(i);
        tfinal = rt_on(i+1);
   else
       tstart = tfinal;
       tfinal = max(max(xdata)*1.1,tstart+6*7*24);
   end
     
end
yall = sum(yout(1:2*n_subpops,:),1);

% only use values which are not nan or inf
inds = ~isinf(yall)&~isnan(yall);
yall = yall(inds);
yout = yout(1:2*n_subpops,inds);
tout = tout(inds);
[uni_tout,uniqueinds]=unique(tout);
if length(uniqueinds)~=length(yall)
    yall=yall(uniqueinds);
    tout = uni_tout;
    yout = yout(1:2*n_subpops,uniqueinds);
end

% Interpolate at data points
try
    ypred = interp1(tout,yall,xdata);
catch
    disp(tout(inds));
end

% Get regrowth times - this has to be after the minimum vol has been
% reached
[minvol,indminvol] = min(yall);
t_vend = NaN;%interp1(yall(indminvol:end),tout(indminvol:end),v_end);
t_nadir = NaN;%interp1(yall(indminvol:end),tout(indminvol:end),1.2*minvol);

end