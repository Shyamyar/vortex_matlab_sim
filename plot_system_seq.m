%% Plot Given UAV States

function fig = plot_system_seq(xhist)
plot_rot = [0 1 0;
            1 0 0;
            0 0 -1];

states0_plot = plot_rot * xhist(1, 1:3)';

fontsize = 16;
fig = figure('DefaultAxesFontSize', fontsize);
scatter3(states0_plot(1), states0_plot(2), states0_plot(3), 20, '*')
grid on
hold on