%--------------------------------------------------------------------------
%ODE system (takes parameters (pars) as additional input)
%--------------------------------------------------------------------------

function ydot = model_ODE_pat_RT_Bev_Pem(t,y,pars)
global n_subpops

%---------------------------------States-----------------------------------
Vd = y(1:n_subpops);
Vl  = y(n_subpops+1:2*n_subpops);
gamma = y(2*n_subpops+1);


%-------------------------------Parameters---------------------------------
%NOTE: needs to be in same order as defined on line 15 of model_pars.m
% lambda = pars(2*n_subpops+1:3*n_subpops);	
lambda = pars(2*n_subpops+2); % h^-1
lambda_die = pars(2*n_subpops+3);%lambda;%pars(31:40);	
epsilon = pars(2*n_subpops+1);


%----------------------------------ODEs------------------------------------

dgamma = -epsilon*gamma;
dVl = lambda.*Vl'-gamma.*Vl';
dVd = -lambda_die.*Vd';

   
ydot = [dVd'; dVl';dgamma'];
end