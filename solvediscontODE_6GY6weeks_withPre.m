function [ypred,yout,tout,t_vend,t_nadir,gammaout] = solvediscontODE_6GY6weeks_withPre(pars,Init,xdata,v_end)
global time_firstRT n_subpops useOnlyS maxNumberOfFractions useImmuneRT useRelImmune n_weeks_btwFrac
rt_on = [0:1:maxNumberOfFractions]*7*n_weeks_btwFrac*24+time_firstRT; % in h
dose = 6; % in Gy


t_final_plot  = 1500*24;
daysShift_max = 42; 

options = [];%odeset('Events',@eventsPat,'OutputFcn',@odeplot,'OutputSel',1,...
   %'Refine',refine);
   
   
   
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



% surviving fraction
if useOnlyS
    S = pars(1:n_subpops);
else
    alpha = pars(1:n_subpops);
    beta = pars(n_subpops+1:2*n_subpops);
    S = exp(-alpha*dose-beta*dose^2);
end


% First fraction at time 0
if ismember(0,rt_on)
    v_viab = Init(n_subpops+1:2*n_subpops);
    v_dying = Init(1:n_subpops);
    v_dying = v_dying+(1-S).*v_viab;
    v_viab = v_viab-(1-S).*v_viab;
    gamma = Init(2*n_subpops+1);
    Init=[v_dying,v_viab,gamma];
    rt_on =rt_on(2:end);
    inds = [0];
else
    inds = [];
end
numfractions = length(rt_on);
v_0 = sum(Init(1:2*n_subpops));
immune100 = (1-S^5)*v_0;

   
   
tstart = 0;
tfinal = rt_on(1);
all_tstart(1)=tstart;
for i = 1:numfractions

   [t,y] = ode45(@model_ODE_pat_RT_Bev_Pem,[tstart:1:tfinal],Init,options,pars);
   y = transpose(y);
   t = transpose(t);
    
   % Accumulate output and set the new initial conditions
   v_viab = y(n_subpops+1:2*n_subpops,end)';
   v_dying =y(1:n_subpops,end)';
   
   if useImmuneRT && i<numfractions
       epsilon = pars(2*n_subpops+1);
       if useRelImmune
           immune_thisFraction = (1-S).*v_viab./immune100;
           day_shift = immune_thisFraction*daysShift_max;
           if day_shift >daysShift_max
               gamma = y(1+2*n_subpops,end)'*exp(epsilon*daysShift_max*24);
           else
               gamma = y(1+2*n_subpops,end)'*exp(epsilon*day_shift*24);
           end
       else
            gamma = y(1+2*n_subpops,end)'*exp(epsilon*10*24);
       end
   else
       gamma = y(1+2*n_subpops,end)';
   end
   
   v_dying = v_dying +(1-S).*v_viab;
   v_viab = v_viab-(1-S).*v_viab;
   Init=[v_dying,v_viab,gamma];
   tout = [tout t(2:end)];
   yout = [yout y(:,2:end)];
   inds = [inds,length(tout)];
   
      
%    if i==1
%        volRTon=sum(yout(1:2,end));
%    end
   
  
   if i<numfractions-1
        tstart = rt_on(i);
        tfinal = rt_on(i+1);
   else
       tstart = rt_on(i);
       tfinal = t_final_plot; %max(max(xdata)*1.1,tstart+6*6*7*24);
   end
   all_tstart(length(all_tstart)+1)=tstart;
   
     
end
yall = nansum(yout(1:2*n_subpops,:),1);

% Potentilly simulate up to 60 more weeks (max 10x6 weeks more)
countex = 0;
while(max(yall)<v_end) && countex <10
    v_viab = y(n_subpops+1:2*n_subpops,end)';
    v_dying =y(1:n_subpops,end)';
    gamma = y(1+2*n_subpops,end)';
    Init=[v_dying,v_viab,gamma];
    tstart = tout(end);
    tfinal = tstart+6*n_weeks_btwFrac*7*24;
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
   inds = [inds,length(tout)];
   yall = nansum(yout(1:2*n_subpops,:),1);
   countex = countex+1;
end




% only use values which are not nan or inf
inds = ~isinf(yall)&~isnan(yall);
yall = yall(inds);
tout = tout(inds);
gammaout = yout(2*n_subpops+1,:); 
yout = yout(1:2*n_subpops,inds);
[uni_tout,uniqueinds]=unique(tout);
if length(uniqueinds)~=length(yall)
    yall=yall(uniqueinds);
    tout = uni_tout;
    yout = yout(1:2*n_subpops,inds);
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
% reached
minvol = yall(inds(1));
% k=2;
% while k<length(rt_on)
%     
%    minvol_new = yall(inds(k));
%    if minvol_new<minvol
%        minvol = minvol_new;
%    else
%        v_nadir = interp1(yall(inds(k-1):inds(k)),tout(inds(k-1):inds(k)),1.2*minvol);
%    end
%    
%    
%    k=k+1;
% end

if countex==10
    t_vend = 100000;
else
    try

    [a,indsRT]=ismember(all_tstart,tout);
    indsRT = indsRT(a);
    
    vol_atRTstart = yall(indsRT);
    [minVolRT,ind_minVolRT] = min(vol_atRTstart);
    
    inds = find(yall>=v_end & tout>tout(indsRT(ind_minVolRT)));
    ind = inds(1);
    tmp = find(ind>indsRT);
    try
        ind_st = indsRT(tmp(end));
        if length(indsRT)<tmp(end)+1
            ind_en = length(tout);
        else
            ind_en = indsRT(tmp(end)+1);
        end

    catch
        ind_RT =1;
        ind_st = indsRT(ind_RT);
        try
            ind_en = indsRT(ind_RT+1);
        catch
            ind_en = length(tout);
        end

    end

    try
        t_vend = interp1(yall(ind_st:ind_en),tout(ind_st:ind_en),v_end);
    catch
        t_vend = 100000;
    end
catch
     t_vend = 100000;
end
    
% inds = find(yall>=v_end);
% ind = inds(1);
% [a,indsRT]=ismember(all_tstart,tout);
% indsRT = indsRT(a);
% tmp = find(ind>indsRT);
% try
%     ind_st = indsRT(tmp(end));
% 
%     if length(indsRT)<tmp(end)+1
%         ind_en = length(tout);
%     else
%         ind_en = indsRT(tmp(end)+1);
%     end
% 
% catch
%     ind_RT =1;
%     ind_st = indsRT(ind_RT);
% %     if ind_st == 1
% %         ind_RT = ind_RT+1;
% %         ind_st = indsRT(ind_RT);
% %     end
% 
%     try
%         ind_en = indsRT(ind_RT+1);
%     catch
%         ind_en = length(tout);
%     end
% 
% end
% try
%     t_vend = interp1(yall(ind_st:ind_en),tout(ind_st:ind_en),v_end);
% catch
%     t_vend = 100000;
% end
end
t_nadir = NaN;%interp1(yall(indminvol:end),tout(indminvol:end),1.2*minvol);

end