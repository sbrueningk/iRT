%--------------------------------------------------------------------------
%Costfunction used to optimize 
%--------------------------------------------------------------------------

function J = costfun_RT_Bev_Pem(pars,data)


epsilon = pars(1);
gamma0  = pars(2);
lambda  = pars(3);
s       = pars(4);
vpre    = pars(5);


% Calculate prediction
ydata = data.ydata;
ypred = analyticSol_RT_Bev_Pem_eval_withPre(data.xdata*24, epsilon,gamma0,lambda,s,vpre);


% Relative error
err = (ypred - ydata)./ydata; 

% Ignore 0 points
err(ydata==0) = (ypred(ydata==0) - ydata(ydata==0))/0.5; 
                              
% Weigh small data point more than the rest
weights = ones(size(err));% mean(ydata)./(ydata);%1./(ydata).^2;%
weights(abs(ypred - ydata)<0.5) = 0.0001;
weights(end) = 1; % Even if the last datapoint is small it should be accounted for
if ydata(end)<0.3 
    weights(end) = 2;
    if err(end) < -0.5
        weights(end) = 2000;
    end
end
if any(abs(ypred(ydata>2.5) - ydata(ydata>2.5))./ydata(ydata>2.5)>0.3)
    inds1 = find(ydata>2.5);
    inds = inds1(find(abs(ypred(ydata>2.5) - ydata(ydata>2.5))./ydata(ydata>2.5)>0.3));
    err(inds) = abs(ypred(inds) - ydata(inds));
end
if ypred(end)<0.1 
    if weights(end)<20
        weights(end) = 200;
    end
end


J = err'*(weights.*err); 














end