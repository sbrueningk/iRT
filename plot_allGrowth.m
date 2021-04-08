function  plot_allGrowth(allDataout_standard,i,il,expName,uncertData)

patient = i; 
n_subpops = 1;
v0 =  allDataout_standard{i,il,1}.v0;
n_vary = size(uncertData,3);
n_maxBS = 50;

if n_vary>1
    n_vary = n_maxBS;
    volData_vary=NaN(n_vary-1,length(uncertData{i,il,2}(:,2)));
end


    % Get the bounds for each model
    counter = 1;
    for iv =1:size(allDataout_standard,3)
        
        if iv ==1
            all_V_standard = allDataout_standard{i,il,iv}.v_standard;
            all_V_PART = allDataout_standard{i,il,iv}.v_PART;
            all_V_6weeks = allDataout_standard{i,il,iv}.v_6week;
            all_rs = allDataout_standard{i,il,iv}.rs;
            all_r2 = allDataout_standard{i,il,iv}.r2;
        else
            
            
            % Ensure fit has converged: 
            if allDataout_standard{i,il,iv}.r2 >0.6 
                if length(all_rs)<n_maxBS
                    
                    volData_vary(counter,:)=uncertData{i,il,iv}(:,2);
                    
                    all_V_standard = [all_V_standard;interp1(allDataout_standard{i,il,iv}.t_standard,...
                        allDataout_standard{i,il,iv}.v_standard,allDataout_standard{i,il,1}.t_standard)];
                    all_V_PART = [all_V_PART;interp1(allDataout_standard{i,il,iv}.t_PART,...
                        allDataout_standard{i,il,iv}.v_PART,allDataout_standard{i,il,1}.t_PART)];
                    all_V_6weeks = [all_V_6weeks;interp1(allDataout_standard{i,il,iv}.t_6week,...
                        allDataout_standard{i,il,iv}.v_6week,allDataout_standard{i,il,1}.t_6week)];
                    all_rs = [all_rs;allDataout_standard{i,il,iv}.rs];
                    all_r2 = [all_r2;allDataout_standard{i,il,iv}.r2];
                    counter = counter + 1;
                end
            end
        end
    end
    min_v_standard = min(all_V_standard);
    max_v_standard = max(all_V_standard);
    min_v_PART = min(all_V_PART);
    max_v_PART = max(all_V_PART);
    min_v_6weeks = min(all_V_6weeks);
    max_v_6weeks = max(all_V_6weeks);
    
  
       
    % Bootstrap Uncertainty
    if n_vary>1
        err_data = nanstd(volData_vary);
    else
        err_data = zeros(size(allDataout_standard{i,il,1}.t_data));
    end
        
    figure
    hold on
    dataplot = errorbar(allDataout_standard{i,il,1}.t_data,allDataout_standard{i,il,1}.v_data,err_data,...
        'ko','linewidth',1,'MarkerFaceColor','k', 'color', 'k', 'MarkerSize',8);
    plot(allDataout_standard{i,il,1}.t_standard,...
        allDataout_standard{i,il,1}.v_standard,'r-','linewidth',3);
    plot(allDataout_standard{i,il,1}.t_6week,...
    allDataout_standard{i,il,1}.v_6week,'b--','linewidth',3);
    plot(allDataout_standard{i,il,1}.t_PART,...
    allDataout_standard{i,il,1}.v_PART,':','linewidth',3, 'color', [0 0.5 0]);
    uistack( dataplot,'top')

    try
    x2 = [allDataout_standard{i,il,1}.t_standard, fliplr(allDataout_standard{i,il,1}.t_standard)];
    inBetween = [min_v_standard, fliplr( max_v_standard)];
    h = fill(x2, inBetween, 'r','LineStyle','none');
    set(h,'facealpha',.3)
    
        x2 = [allDataout_standard{i,il,1}.t_PART, fliplr(allDataout_standard{i,il,1}.t_PART)];
    inBetween = [min_v_PART, fliplr( max_v_PART)];
    h = fill(x2, inBetween,[0 0.5 0],'LineStyle','none');
    set(h,'facealpha',.3)
    
        x2 = [allDataout_standard{i,il,1}.t_6week, fliplr(allDataout_standard{i,il,1}.t_6week)];
    inBetween = [min_v_6weeks, fliplr( max_v_6weeks)];
    h = fill(x2, inBetween, 'b','LineStyle','none');
    set(h,'facealpha',.3)
    catch
        a=1;
    end
    

    ylabel('Volume [cm^3]')
    xlabel('Time [days]')
    legend('Data',sprintf('Fit, RMSE = %.2f',allDataout_standard{i,il,1}.rs)...
        ,'Prediction iRT','Prediction iRT+Boost',...
        'Location','best')
    set(gca,'ylim',[0,1.3*max([allDataout_standard{i,il,iv}.v_data',v0])], 'Fontsize', 20)
    set(gca,'xlim',[min(allDataout_standard{i,il,1}.t_data)-14,900], 'Fontsize', 20)
    title(['Patient ' num2str(patient)])
    box on
    legend boxoff
    set(gcf,'position',[100,100,500,450])

    currFolder = pwd;
    cd(expName)
    filename = [num2str(patient) '_ALL' num2str(n_subpops) '.fig'];
    savefig(filename)
    filename = [num2str(patient)  '_ALL' num2str(n_subpops) '.png'];
    saveas(gcf,filename)
    filename = [num2str(patient)  '_ALL' num2str(n_subpops)];
    saveas(gcf,filename,'epsc')
    cd(currFolder)




end