function plotRGTimeEvolution(t,r,w,x1,x2,u,stateMax,stateMin, ...
    inputMax,inputMin,kappaHistory)
% Create one figure with four subplots
figure('Color','w', ...
       'Name','RG - Time Evolution and Constraints', ...
       'NumberTitle','off');

tiledlayout(4,1,'TileSpacing','compact','Padding','compact');

%% 1. Position and references
nexttile;
hold on; grid on; box on;

plot(t,r,'g--','LineWidth',1.5);
plot(t,w,'b','LineWidth',2);
plot(t,x1,'r','LineWidth',2);

yline(stateMax(1),'m--','LineWidth',1.2);
yline(stateMin(1),'m--','LineWidth',1.2);

ylabel('Position / Reference');
title('Position Tracking');

legend('Desired reference r', ...
       'Accepted reference w', ...
       'Position x_1', ...
       'Position constraints', ...
       'Location','best');

xlim([t(1) t(end)]);

%% 2. Velocity
nexttile;
hold on; grid on; box on;

plot(t,x2,'b','LineWidth',2);

yline(stateMax(2),'r--','LineWidth',1.3);
yline(stateMin(2),'r--','LineWidth',1.3);

ylabel('Velocity x_2');
title('Velocity Constraint');

legend('Velocity x_2', ...
       'Maximum velocity', ...
       'Minimum velocity', ...
       'Location','best');

xlim([t(1) t(end)]);

%% 3. Command
nexttile;
hold on; grid on; box on;

stairs(t,u,'b','LineWidth',1.8);

yline(inputMax,'r--','LineWidth',1.3);
yline(inputMin,'r--','LineWidth',1.3);

ylabel('Control input u');
title('Control Input Constraint');

legend('Command u', ...
       'Maximum input', ...
       'Minimum input', ...
       'Location','best');

xlim([t(1) t(end)]);

%% 4. Reference Governor parameter
nexttile;
hold on; grid on; box on;

stairs(t(1:end-1),kappaHistory,'b','LineWidth',1.8);

yline(1,'r--','LineWidth',1.2);
yline(0,'r--','LineWidth',1.2);

ylabel('\kappa(k)');
xlabel('Time [s]');
title('Reference Governor Parameter K');

legend('\K(k)', ...
       'Upper bound', ...
       'Lower bound', ...
       'Location','best');

ylim([-0.05 1.05]);
xlim([t(1) t(end)]);

sgtitle('Dynamic Reference Governor - Time Evolution');

end
