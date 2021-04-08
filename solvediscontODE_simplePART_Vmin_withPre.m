function [ypred,yout,tout,t_vend,gammaout] = solvediscontODE_simplePART_Vmin_withPre(pars,Init,xdata,v_end)
global time_firstRT n_subpops useOnlyS maxNumberOfFractions useImmuneRT useRelImmune n_weeks_btwFrac

t_final_plot     = 1500*24;
rt_on            = time_firstRT; % in h
max_numfractions = maxNumberOfFractions;
daysShift_max    = 42;


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



% Surviving fraction
if useOnlyS
    S_in = pars(1:n_subpops);
else
    error('Option useOnlyS=false is not supported!')
end
S = S_in;
prevFractionWas0 = false;

v_0 = sum(Init(1:2*n_subpops));
immune100 = (1-S^5)*v_0;
% First fraction at time 0
if ismember(0,rt_on)
    v_viab = Init(1+n_subpops:2*n_subpops);
    v_dying = Init(1:n_subpops);
    v_dying = v_dying+(1-S).*v_viab;
    v_viab = v_viab-(1-S).*v_viab;
    gamma = Init(2*n_subpops+1);
    Init=[v_dying,v_viab,gamma];
    %rt_on =rt_on(2:end);
    tstart = 0;
    tfinal = n_weeks_btwFrac*7*24;
    i=2;
    vmin = v_viab+v_dying;
    inds = [0];
    v0 = vmin;
else
    tstart = 0;
    tfinal = rt_on(1);
    tout = [0];
    yout = [Init'];
    i=1;
    inds = [];
end


   
   
all_tstart(1)=tstart;
exit_flag = 0;
counter = 0;
while i<=max_numfractions+1 && counter<2  && i < 100  
   sol = ode45(@model_ODE_pat_RT_Bev_Pem,[tstart:1:tfinal],Init,options,pars);
   [t,y] = ode45(@model_ODE_pat_RT_Bev_Pem,[tstart:1:tfinal],Init,options,pars);
   y = transpose(y);
   t = transpose(t);
    
   % Accumulate output and set the new initial conditions
   v_viab = y(1+n_subpops:2*n_subpops,end)';
   v_dying =y(1:n_subpops,end)';
   tout = [tout t(2:end)];
   yout = [yout y(:,2:end)];
   inds = [inds,length(tout)];
   
      
   if i==1
       vmin=sum(yout(1:2,end));
       v_end6weeks = 0;
       tstart = rt_on;
       v0 = vmin;
   else
       v_end6weeks = sum(yout(1:2,end));
       tstart = tfinal;
   end
   
   % Boost if 20% above min. achieved and measured volume
   if v_end6weeks<vmin*1.2 && exit_flag==0 
       vmin = v_end6weeks;
       if i<max_numfractions
           tfinal = tstart+n_weeks_btwFrac*7*24;
       else
           tfinal = t_final_plot; 
       end
   else
       if prevFractionWas0
           vmin = v_end6weeks;
           if i<max_numfractions
               tfinal = tstart+n_weeks_btwFrac*7*24;
           else
               tfinal = t_final_plot; 
           end  
       else
       
           if exit_flag<3
               tfinal = tstart+1*24;
               exit_flag = exit_flag+1;
               prevFractionWas0 = false;
           else
               tfinal = t_final_plot; %max(max(xdata)*1.1,tstart+12*6*7*24);
               S = 1;
               counter = counter+1;
           end
       end
   end 
      
   
   % Prepare next fraction
   if v_end6weeks<-1
       S = 1;
       max_numfractions = max_numfractions+1;
       prevFractionWas0 = true;
   else
       S = S_in;
       prevFractionWas0 = false;
   end
   
   v_dying = v_dying +(1-S).*v_viab;
   v_viab = v_viab-(1-S).*v_viab;
   if useImmuneRT && S<1 %&& i<maxNumberOfFractions
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
   Init=[v_dying,v_viab,gamma];
   all_tstart(length(all_tstart)+1)=tstart;
   i=i+1;
end
yall = nansum(yout(1:2*n_subpops,:),1);
[~,minInd] = min(yall);



% Potentially simulate up to 60 more weeks to ensure covering everything
% (max 10x6 weeks more here)
countex = 0;
while(max(yall(minInd:end))<v_end) && countex <10
    v_viab = y(n_subpops+1:2*n_subpops,end)';
    v_dying =y(1:n_subpops,end)';
    gamma = y(1+2*n_subpops,end)';
    Init=[v_dying,v_viab,gamma];
    tstart = tout(end);
    tfinal = t_final_plot;%tstart+6*6*7*24;
    if tfinal-2>tstart
    sol = ode45(@model_ODE_pat_RT_Bev_Pem,[tstart:1:tfinal],Init,options,pars);
    [t,y] = ode45(@model_ODE_pat_RT_Bev_Pem,[tstart:1:tfinal],Init,options,pars);
    y = transpose(y);
    t = transpose(t);
    
   % Accumulate output and set the new initial conditions
   tout = [tout t(2:end)];
   yout = [yout y(:,2:end)];
   inds = [inds,length(tout)];
   yall = nansum(yout(1:2*n_subpops,:),1);
    else
        countex = 10; 
    end
   
   countex = countex+1;
end


% Only use values which are not nan or inf
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
try

    [a,indsRT]=ismember(all_tstart,tout);
    indsRT = indsRT(a);
    
    vol_atRTstart = yall(indsRT);
    [~,ind_minVolRT] = min(vol_atRTstart);
    
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
    disp('Could not estimate time to regrowth!!')
     t_vend = 100000;
end



