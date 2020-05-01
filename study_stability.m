clear all

global ds dtau D_L Gr p b tau_max St_m ;
ds = 0.01 ;
Parameters() ;
withPlot = false ;

Gr = 20e10 ;

nSt = 50 ;
listSt = linspace(.5, 12, nSt);

stability = zeros(nSt, 1) ;

figure(2)
hold on
set(gca, 'FontSize',18, 'YGrid', 'on', 'XScale','lin', 'YScale','lin') ;
xlabel('St_m'); 
ylabel('\sigma evolution (%)');
ylim([-100 1500])
legend show
pause(0);
hold on

for i=1:nSt
    
    St_m = listSt(i) ;
    disp("Computing steady state " + num2str(i) + "/" + num2str(nSt) + " ...")
    [theta, thetaSS] = steady_state(false) ;
    int_theta = loop_integration_dz(thetaSS);
    Re = (2 * D_L * Gr * int_theta / p) ^ (1 / (3-b)) ; 
    disp("                           done")    
     
    disp("Computing transient...")
    evolW = solveur(Re, 1, thetaSS, false);
    disp("                           done")
        
    stability(i) = check_stability(evolW, tau_max/dtau/2) * 100 ;
        
end
    
lines = plot(listSt, stability, 'LineWidth', 2, 'DisplayName', "Gr = " + num2str(Gr/1e10) + "e10");
    
sound(sin(1:15000));















function prepFig(nCells, dtau, s, nLeg1, nCool, nLeg2)
    figure(3) 
    clf
    hold on
    
%     set(gcf, 'Position', [10 100 800 500]);
    set(gca, 'FontSize',18, 'YGrid', 'on') ;
    
    title("Steady State of the \theta profile - " +... 
        num2str(nCells) + " cells - d\tau = " + num2str(dtau) );
    
    xlabel('Adimensionnal Position');
    xlim([s(1) s(end)]) ;
    
    ylabel('Adimensionnal Temperature \theta');
    
    plot_zones_limits(3, s, nLeg1, nCool, nLeg2)
    pause(1);
end

function plot_zones_limits(idfig, s, nLeg1, nCool, nLeg2)

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
    annotation(figure(idfig) ,'textbox',[0.18625 0.134 0.081875 0.06],...
        'String',{'Leg 1'},...
        'HorizontalAlignment','center',...
        'FontSize',16,...
        'EdgeColor','none');

    % Create textbox
    annotation(figure(idfig) ,'textbox',[0.375625 0.134 0.090625 0.06],...
        'String',{'Cooler'},...
        'HorizontalAlignment','center',...
        'FontSize',16,...
        'EdgeColor','none');

    % Create textbox
    annotation(figure(idfig) ,'textbox',[0.573750000000001 0.13 0.081875 0.06],...
        'String',{'Leg 2'},...
        'HorizontalAlignment','center',...
        'FontSize',16,...
        'EdgeColor','none');
    % Create textbox
    annotation(figure(idfig),'textbox',[0.763125000000001 0.128 0.090625 0.06],...
        'String',{'Heater'},...
        'HorizontalAlignment','center',...
        'FontSize',16,...
        'EdgeColor','none');
end
