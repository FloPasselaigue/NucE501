function [theta, thetaSS] = steady_state(withPlot)
    % Compute the steady state for temperature

        % Get the parameters needed in this particular workspace 
    global s ds nLeg1 nLeg2 nCool nHeat ;

        % Mass flow rate at steady state
    w_n = 1 ;
    w_n_1 = 1 ;

        % Build resolution matrices 
    [M, N, P, ~] = build_matrices(w_n, w_n_1);
    M_inv = M \ eye(length(M)) ;
    one = ones(length(M(:,1)), 1) ;

        % Initial state
    theta = .5 + .5*(sin(s*2*pi/s(end)))' ;
        % Resolution
    theta_n_1 = M \ (N*theta + P*one) ;


       % Time loop initialization
    iter = 1 ;
    r1 = [average(abs(theta(1:nLeg1)-theta_n_1(1:nLeg1)))] ;
    r2 = [average(abs(theta(nLeg1+nCool+1:nLeg1+nCool+nLeg2)-theta_n_1(nLeg1+nCool+1:nLeg1+nCool+nLeg2)))] ;

        % Plot preparation
    if withPlot
        prepFig(1, s, nLeg1, nCool, nLeg2)
        iColor = 0 ;
        colors = [[2,56,88]; [4,90,141]; [5,112,176]; [54,144,192]; [116,169,207]; 
                [166,189,219]; [208,209,230]; [153,216,201]; [65,174,118]]/255;
        hold on
        plot(s, theta, 'LineWidth', 2, 'DisplayName', "Step 0", 'Color', colors(1,:));
        pause(0);
    end
     
        % Time loop
    while r1(end) > 5e-4 || r2(end) > 5e-4 
            % update theta
        theta = theta_n_1;

            % Resolution
        theta_n_1 = M_inv * (N*theta + P*one) ;

            % Convergence criterion
        r1 = [r1 average(abs(theta(1:nLeg1)-theta_n_1(1:nLeg1)))] ;
        r2 = [r2 average(abs(theta(nLeg1+nCool+1:nLeg1+nCool+nLeg2)-theta_n_1(nLeg1+nCool+1:nLeg1+nCool+nLeg2)))] ;

        if withPlot && rem(iter, 50) == 0
            iColor = iColor + 1 ;
            figure(1)
            hold on
            plot(s, theta, 'LineWidth', 2, 'DisplayName', "Step " + num2str(iter), 'Color', colors(1+rem(iColor, length(colors)),:));
            pause(0);
        end
            
        iter=iter+1;
    end

        % Compute smooth profile
    t1 = average(theta(1:nLeg1)) ;
    t2 = max(average(theta(nLeg1+nCool+1:nLeg1+nCool+nLeg2)), 1e-20) ;
    s1 = nLeg1 * ds ;
    sCool = nCool * ds ;
    B = - log(t1 / t2) / sCool ;
    A = t2*(t1/t2)^((s1+sCool)/sCool) ;
    thetaLeg1 = t1 * ones(nLeg1,1) ;
    thetaCool = A * exp(B * s(nLeg1+1:nLeg1+nCool)') ;
    thetaLeg2 = t2 * ones(nLeg2,1) ;
    thetaHeat = linspace(thetaLeg2(end), thetaLeg1(1), nHeat)' ;
    
    thetaSS = [thetaLeg1 ; thetaCool; thetaLeg2; thetaHeat];
    
end


%% Auxiliary functions
function prepFig(id, s, nLeg1, nCool, nLeg2)
    figure(id) 
    clf
    hold on
    set(gca, 'FontSize',18, 'YGrid', 'on') ;
    xlabel('Adimensionnal Position');
    xlim([s(1) s(end)]) ;
    ylabel('Adimensionnal Temperature \theta');
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
