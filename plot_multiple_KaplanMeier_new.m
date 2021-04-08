
function plot_multiple_KaplanMeier(timesOfDeath,timesOfDeath_u, timesOfDeath_l)

numberofRuns = size(timesOfDeath,1);
numberOfSchemes = size(timesOfDeath,2);
times = 0:24*7:max(timesOfDeath(:));

figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,...
    'Position',[0.13 0.136031748847356 0.56261744966443 0.788968251152644]);
hold(axes1,'on');

set(gca,'Fontsize',16,'xLim',[0,max(timesOfDeath(:))/(7*24)])

legends = cell(1,numberOfSchemes);
percentAlive = 100*NaN(length(times),numberOfSchemes);
percentAlive_l = 100*NaN(length(times),numberOfSchemes);
percentAlive_u = 100*NaN(length(times),numberOfSchemes);
for s=1:numberOfSchemes

    for n=1:numberofRuns
        
        for t =1:length(times)
            percentAlive(t,s) = 100*length(find(timesOfDeath(:,s)>times(t)))/numberofRuns/...
                                (sum(~isnan(timesOfDeath(:,s)))/length(timesOfDeath(:,s)));
            percentAlive_l(t,s) = 100*length(find(timesOfDeath_l(:,s)>times(t)))/numberofRuns/...
                                (sum(~isnan(timesOfDeath_l(:,s)))/length(timesOfDeath_l(:,s)));
            percentAlive_u(t,s) = 100*length(find(timesOfDeath_u(:,s)>times(t)))/numberofRuns/...
                                (sum(~isnan(timesOfDeath_u(:,s)))/length(timesOfDeath_u(:,s)));
        end
    end
    
    sh_l = stairs(percentAlive_l(:,s),'LineStyle', 'none');
    sh_u = stairs(percentAlive_u(:,s),'LineStyle', 'none');
    set(get(get(sh_l,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    set(get(get(sh_u,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    curve2 = [repelem(sh_u.YData(1:end-1),2),sh_u.YData(end)];
    curve1 = [repelem(sh_l.YData(1:end-1),2),sh_l.YData(end)];
%     plot(times/(7*24), curve1, 'b', 'LineWidth', 2);
%     plot(times/(7*24), curve2, 'b', 'LineWidth', 2);
    x2 = [[sh_l.XData(1),repelem(sh_l.XData(2:end),2)],...
        fliplr([sh_u.XData(1),repelem(sh_u.XData(2:end),2)])];%[times/(7*24), fliplr(times/(7*24))];
    inBetween = [curve1, fliplr(curve2)];

if s==1
        stairs(percentAlive(:,s),'LineWidth',3,'color','r',...
            'DisplayName','HFSR');
        
    elseif s==2

        stairs(percentAlive(:,s),'LineWidth',3,'color','b',...
            'DisplayName','iRT');
    elseif s==3

        stairs(percentAlive(:,s),'LineWidth',3,'color',[0 0.5 0],...
            'DisplayName','iRT+boost');
    end
end

for s=1:numberOfSchemes

    for n=1:numberofRuns
        
        for t =1:length(times)
            percentAlive(t,s) = 100*length(find(timesOfDeath(:,s)>times(t)))/numberofRuns/...
                                (sum(~isnan(timesOfDeath(:,s)))/length(timesOfDeath(:,s)));
            percentAlive_l(t,s) = 100*length(find(timesOfDeath_l(:,s)>times(t)))/numberofRuns/...
                                (sum(~isnan(timesOfDeath_l(:,s)))/length(timesOfDeath_l(:,s)));
            percentAlive_u(t,s) = 100*length(find(timesOfDeath_u(:,s)>times(t)))/numberofRuns/...
                                (sum(~isnan(timesOfDeath_u(:,s)))/length(timesOfDeath_u(:,s)));
        end
    end
    
    sh_l = stairs(percentAlive_l(:,s),'LineStyle', 'none');
    sh_u = stairs(percentAlive_u(:,s),'LineStyle', 'none');
    set(get(get(sh_l,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    set(get(get(sh_u,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    curve2 = [repelem(sh_u.YData(1:end-1),2),sh_u.YData(end)];
    curve1 = [repelem(sh_l.YData(1:end-1),2),sh_l.YData(end)];
%     plot(times/(7*24), curve1, 'b', 'LineWidth', 2);
%     plot(times/(7*24), curve2, 'b', 'LineWidth', 2);
    x2 = [[sh_l.XData(1),repelem(sh_l.XData(2:end),2)],...
        fliplr([sh_u.XData(1),repelem(sh_u.XData(2:end),2)])];%[times/(7*24), fliplr(times/(7*24))];
    inBetween = [curve1, fliplr(curve2)];

%     

    
    if s==1
        h=fill(x2, inBetween, 'r','Edgecolor','none');
        set(h,'facealpha',.3)
        set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        uistack( h,'bottom')
        
    elseif s==2
        h=fill(x2, inBetween, 'b','Edgecolor','none');
        set(h,'facealpha',.3)
        set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        uistack( h,'bottom')

    elseif s==3
        h=fill(x2, inBetween, 'g','Edgecolor','none');
        set(h,'facealpha',.3)
        set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        uistack( h,'bottom')
    end
end



% Create ylabel
ylabel('Percent below cut-off');

% Create xlabel
xlabel('Time [days]');

set(gca,'xlim',[0,128.5])

% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes1,[0 128.571428571429]);
box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontSize',16);
% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.698756775567362 0.707487599053082 0.297651011767803 0.218809528714135]);


end


