function [startP,LB,UB] = getParametersAnalytical(vpre, thispat,usePrevRes)

    global n_subpops useOnlyS useFixedLambda useFixedS fixedS fixedLambda 
    global fixedEpsilon useFixedEpsilon


    % Surviving fraction
    if useFixedS
        S_0 = fixedS;
        S_l = fixedS;
        S_h = fixedS;
    else
        S_0 = 0.293;
        S_l = 0.005;
        S_h = 1;
    end

    % Growth rate
    if useFixedLambda
        lambda_0   = fixedLambda;
        lambda_l   = fixedLambda;
        lambda_h   = fixedLambda;
    else
        lambda_0   = log(2)./(20*24)*ones(1,n_subpops);
        lambda_l   = log(2)./(50*24)*ones(1,n_subpops);
        lambda_h   = log(2)./(8*24)*ones(1,n_subpops);
    end
    
    % Death rate
    lambda_die_0   = log(2)./(20*24)*ones(1,n_subpops);
    lambda_die_l   = log(2)./(50*24)*ones(1,n_subpops);
    lambda_die_h   = log(2)./(2*24)*ones(1,n_subpops);
    
    % gamma0
    gamma0_0   = 0.0029;
    gamma0_l   = log(2)./(50*24)*ones(1,n_subpops);
    gamma0_h   = log(2)./(8*24)*ones(1,n_subpops);

    % Epsilon - resistance to immunotherapy/bev
    if useFixedEpsilon
        epsilon_0   = fixedEpsilon;
        epsilon_l   = fixedEpsilon;
        epsilon_h   = fixedEpsilon;
    else
        epsilon_0   = 0.0000291;
        epsilon_l   = 0.000000001;
        epsilon_h   = 0.1;
    end
    
    if vpre < 1
        vpre_0 = vpre;
        vpre_l = 0.0001;
        vpre_h = 2.;
    else
        vpre_0 = vpre;
        vpre_l = 0.7*vpre;
        vpre_h = max(2,1.3*vpre);
    end
    
    if usePrevRes
        pats = [1,2,4,6,7,8,10,11,12,16,17,19,21,22,25,29]; 
        
        epsilons_prev = 0.0001*[0.290510364,0.8692779,0.255032328,...
                                1.077920673,0.64477103,1.502443991,...
                                0.43286223,3.371883237,1.306126194,...
                                0.326794705,0.102876258,152.489812,...
                                1.311457956,12.06333942,3.935415676,...
                                3.504801023];
        S_prev        = [0.293268608,0.515798648,0.50975194,0.512299994,...
                         0.505388813,0.842383655,0.308630141,0.293201815,...
                         0.353382898,0.512568203,0.642897667,0.174506875,...
                         0.329158129,0.122055801,0.265754385,0.097665905]; 
        
        indpat    = find(pats ==thispat);
        epsilon_0 = epsilons_prev(indpat);
        S_0       = S_prev(indpat);

    end 


    startP = [epsilon_0,gamma0_0,lambda_0,S_0,vpre_0];
    LB     = [epsilon_l,gamma0_l,lambda_l,S_l,vpre_l];
    UB     = [epsilon_h,gamma0_h,lambda_h,S_h,vpre_h]; 
    
end