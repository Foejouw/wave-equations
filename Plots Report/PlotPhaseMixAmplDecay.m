load("PlotPhaseMixAmplDecay.mat")

%%
ts = linspace(tMin, tMax.*2, Nt.*2);
loglog(ts, steepAmpl, 'red', 'LineWidth',1.2)
hold on
loglog(ts, shallowAmpl, 'blue', 'LineWidth',1.2)
hold on
loglog(ts, ts.^(-3/2), 'LineWidth',1.2, 'Color', [0.9290 0.6940 0.1250])

legend(gca, {strcat('$\max_{z} (b(',num2str(xs(steep)),', z, t))$'), strcat('$\max_{z} (b(',num2str(xs(shallow)),', z, t))$'), '$t^{-3/2}$'}, "Interpreter","latex", 'Location','southwest', 'FontSize',11)
xlim([1/2 8] )
xlabel('t')
ylabel('Maximum amplitude')

%%

ts = linspace(tMin, tMax.*2, Nt.*2);
loglog(ts, steepAmpl , 'red', 'LineWidth',1.2)
hold on
loglog(ts,  0.5.*0.088.*sqrt(1./eta).*ts.^(-1.5), 'Color', [0.9290 0.6940 0.1250], 'LineWidth',1.2)
legend(gca, {strcat('$\max_{z} (b(',num2str(xs(steep)),', z, t))$'), '$\frac{1}{\sqrt{\eta}} 0.044 \cdot t^{-3/2}$'}, "Interpreter","latex", 'Location','southwest', 'FontSize',11)
xlim([1/2 8] )
xlabel('t')
ylabel('Maximum amplitude')

