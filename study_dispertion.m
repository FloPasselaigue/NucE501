clear all

% Get the parameters needed in this particular workspace
global ds nCells dtau s nLeg1 nCool nLeg2 ;
withPlot = false ;

nds = 30 ;
listds = linspace(.001, .02, nds) ;

dev1    = zeros(1, nds);
devHeat = zeros(1, nds);
devCool = zeros(1, nds);
dev2    = zeros(1, nds);

ds = listds(1) ;
Parameters() ;
s1 = s ;
[theta1, thetaSS1] = steady_state(withPlot) ;

dev1(1) = std(theta1(1:nLeg1) - thetaSS1(1:nLeg1));
devCool(1) = std(theta1(nLeg1+1:nLeg1+nCool) - thetaSS1(nLeg1+1:nLeg1+nCool));
dev2(1) = std(theta1(nLeg1+nCool+1:nLeg1+nCool+nLeg2) - thetaSS1(nLeg1+nCool+1:nLeg1+nCool+nLeg2));
devHeat(1) = std(theta1(nLeg1+nCool+nLeg2:end) - thetaSS1(nLeg1+nCool+nLeg2:end));

for i=2:nds-1
    disp(num2str(i) + "/" + num2str(nds))
    ds = listds(i) ;
    Parameters() ;
    [theta, thetaSS] = steady_state(withPlot) ;
    
    dev1(i) = std(theta(1:nLeg1) - thetaSS(1:nLeg1));
    devCool(i) = std(theta(nLeg1+1:nLeg1+nCool) - thetaSS(nLeg1+1:nLeg1+nCool));
    dev2(i) = std(theta(nLeg1+nCool+1:nLeg1+nCool+nLeg2) - thetaSS(nLeg1+nCool+1:nLeg1+nCool+nLeg2));
    devHeat(i) = std(theta(nLeg1+nCool+nLeg2:end) - thetaSS(nLeg1+nCool+nLeg2:end));

end

disp(num2str(i) + "/" + num2str(nds))
ds = listds(end) ;
Parameters() ;
s2 = s ;
[theta2, thetaSS2] = steady_state(withPlot) ;
dev1(end) = std(theta2(1:nLeg1) - thetaSS2(1:nLeg1));
devCool(end) = std(theta2(nLeg1+1:nLeg1+nCool) - thetaSS2(nLeg1+1:nLeg1+nCool));
dev2(end) = std(theta2(nLeg1+nCool+1:nLeg1+nCool+nLeg2) - thetaSS2(nLeg1+nCool+1:nLeg1+nCool+nLeg2));
devHeat(end) = std(theta2(nLeg1+nCool+nLeg2:end) - thetaSS2(nLeg1+nCool+nLeg2:end));

fit1 = polyval(polyfit(listds, dev1, 1), listds) ;
fitCool = polyval(polyfit(listds, devCool, 1), listds) ;
fit2 = polyval(polyfit(listds, dev2, 1), listds) ;
fitHeat = polyval(polyfit(listds, devHeat, 1), listds) ;

% prepFig(1, nCells, dtau, s, nLeg1, nCool, nLeg2)
% hold on
% plot(s1, theta1, 'LineWidth', 2, 'DisplayName', "ds = " + num2str(listds(1)));
% plot(s2, theta2, 'LineWidth', 2, 'DisplayName', "ds = " + num2str(listds(end)));
% legend('show', 'Location', 'southwest') ;

figure(3)
clf
hold on
legend show
set(gca, 'FontSize',18, 'YGrid', 'on', 'XGrid', 'on') ;
xlabel('Adimensionnal Mesh Size');
ylabel('Standard Deviation');
scatter(listds, dev1, 'DisplayName', 'Leg 1');
scatter(listds, devCool, 'DisplayName', 'Cooler');
scatter(listds, dev2, 'DisplayName', 'Leg 2');
scatter(listds, devHeat, 'DisplayName', 'Heater');


function prepFig(id, nCells, dtau, s, nLeg1, nCool, nLeg2)
    figure(id) 
    clf
    hold on
    
%     set(gcf, 'Position', [10 100 800 500]);
    set(gca, 'FontSize',18, 'YGrid', 'on') ;
    
    title("Steady State of the \theta profile - d\tau = " + num2str(dtau) );
    
    xlabel('Adimensionnal Position');
    xlim([s(1) s(end)]) ;
    
    ylabel('Adimensionnal Temperature \theta');
    ylim([-.2 1.8])
    
    plot_zones_limits(id, s, nLeg1, nCool, nLeg2)
end

function plot_zones_limits(id, s, nLeg1, nCool, nLeg2)

    figure(id)
    hold on

    hLine = xline(s(nLeg1),'k') ;
    set(get(get(hLine,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off')
    hLine = xline(s(nLeg1+nCool),'k') ;
    set(get(get(hLine,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off')
    hLine = xline(s(nLeg1+nCool+nLeg2),'k') ;
    set(get(get(hLine,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off')
    % Create textbox
    annotation(figure(id) ,'textbox',[0.21 0.85 0.082 0.06],...
        'String',{'Leg 1'},...
        'HorizontalAlignment','center',...
        'FontSize',16,...
        'EdgeColor','none');

    % Create textbox
    annotation(figure(id) ,'textbox',[0.40 0.86 0.09 0.06],...
        'String',{'Cooler'},...
        'HorizontalAlignment','center',...
        'FontSize',16,...
        'EdgeColor','none');

    % Create textbox
    annotation(figure(id) ,'textbox',[0.58 0.86 0.082 0.06],...
        'String',{'Leg 2'},...
        'HorizontalAlignment','center',...
        'FontSize',16,...
        'EdgeColor','none');
    % Create textbox
    annotation(figure(id),'textbox',[0.79 0.86 0.091 0.06],...
        'String',{'Heater'},...
        'HorizontalAlignment','center',...
        'FontSize',16,...
        'EdgeColor','none');
end







