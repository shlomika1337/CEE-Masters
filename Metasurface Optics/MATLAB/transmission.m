close all; clear; clc;
% Coded by Shlomi Avidan
%% Experimental data plotting
interpolated_data = readmatrix("processed_transmission.csv");
tx = readmatrix("tx.csv");
ty = readmatrix("ty.csv");
phase = readmatrix("phase.csv");
phase_d = readmatrix("processed_phase.csv");

x_d = interpolated_data(:,1);
tx_d = interpolated_data(:,2);
ty_d = interpolated_data(:,3);
x_p = phase(:,1);
x_pd = phase_d(:,1);

x_tx = tx(:,1);
tx_tx = tx(:,2);
x_ty = ty(:,1);
ty_ty = ty(:,2);
p = phase(:,2);
p_d = phase_d(:,2);

figure;
hold on;
plot(x_d,tx_d, '--', 'LineWidth', 2, 'Color', 'blue');
plot(x_d,ty_d, '--', 'LineWidth', 2, 'Color', 'red');
plot(x_tx,tx_tx, '*', 'LineWidth', 5, 'MarkerSize', 12, 'Color', 'blue');
plot(x_ty,ty_ty, '*', 'LineWidth', 5, 'MarkerSize', 12, 'Color', 'red');
xline(550, '-.', 'LineWidth', 3);
xline(490, '-.', 'LineWidth', 3);
yline(0.5, '-.', 'LineWidth', 3);
yline(0.175, '-.', 'LineWidth', 3);
text(650, 0.525, "t_y = 0.5");
text(650, 0.195, "t_x = 0.175");
text(545, 0.535, "\lambda = 550 nm", 'Rotation', 90);
text(485, 0.535, "\lambda = 490 nm", 'Rotation', 90);
grid on;
title({"Experimental transmission coefficient of the nano-beam waveplate,", "for 120-nm width, 100-nm height, 200-nm period"});
xlabel("Wavelength \lambda [nm]");
ylabel("Transmission coefficient");
legend("t_x", "t_y");
fontsize(24,"points");

figure;
hold on;
plot(x_p,p, '*', 'LineWidth', 5, 'MarkerSize', 12, 'Color', 'magenta');
yline(180, '-.', 'LineWidth', 3);
xline(550, '-.', 'LineWidth', 3);
text(600, 190, "Phase = 180-degrees");
text(545, 190, "\lambda = 550 nm", 'Rotation', 90);
grid on;
title({"Experimental relative phase retardation of the nano-beam waveplate,", "for 120-nm width, 100-nm height, 200-nm period"});
xlabel("Wavelength \lambda [nm]");
ylabel("Phase [deg]");
fontsize(24,"points");

%% CoPE & CrPE
% Load the data, skipping the header row
tx_data = readmatrix('tx.csv', 'Range', 2);
ty_data = readmatrix('ty.csv', 'Range', 2);
phase_data = readmatrix('phase.csv', 'Range', 2);

% Separate columns into x and y values
x_tx = tx_data(:,1);
val_tx = tx_data(:,2);

x_ty = ty_data(:,1);
val_ty = ty_data(:,2);

x_phase = phase_data(:,1);
val_phase = phase_data(:,2);

% It's good practice to sort data before interpolation
[x_tx, sortIdx] = sort(x_tx);
val_tx = val_tx(sortIdx);

[x_ty, sortIdx] = sort(x_ty);
val_ty = val_ty(sortIdx);

[x_phase, sortIdx] = sort(x_phase);
val_phase = val_phase(sortIdx);

% --- Interpolation Step ---
% Interpolate ty and phase values onto the tx x-coordinates
ty_aligned = interp1(x_ty, val_ty, x_tx, 'linear', 'extrap');
phase_aligned = interp1(x_phase, val_phase, x_tx, 'linear', 'extrap');

% --- Calculation with Aligned Data ---
% Convert phase from degrees to radians for exp()
phase_rad = deg2rad(phase_aligned);

% Calculate CoPE and CrPE
CoPE = abs(0.5 * (val_tx + ty_aligned .* exp(1j * phase_rad))).^2;
CrPE = abs(0.5 * (val_tx - ty_aligned .* exp(1j * phase_rad))).^2;

% --- Save Results ---
% Create a table with all aligned data and results
results_table = table(x_tx, val_tx, ty_aligned, phase_aligned, CoPE, CrPE, ...
    'VariableNames', {'x', 'tx_aligned', 'ty_aligned', 'phase_aligned_deg', 'CoPE', 'CrPE'});

% Save the results to a new CSV file
writetable(results_table, 'CoPE_CrPE_aligned.csv');
% Read the aligned data from the CSV file
data_table = readtable('CoPE_CrPE_aligned.csv');

%% Create a new figure
figure('Position', [100, 100, 800, 500]);
hold on; % Hold on to plot multiple lines on the same axes

% Plot CoPE and CrPE data
plot(data_table.x, data_table.CoPE, 'LineWidth', 4, 'DisplayName', 'CoPE');
plot(data_table.x, data_table.CrPE,  'LineWidth', 4, 'DisplayName', 'CrPE');
xline(520, '--', 'LineWidth', 3);
xline(500, '--', 'LineWidth', 3);
xline(550, '--', 'LineWidth', 3);
text(495, 0.25, "\lambda = 500 nm", 'Rotation', 90);
text(515, 0.25, "\lambda = 520 nm", 'Rotation', 90);
text(545, 0.25, "\lambda = 550 nm", 'Rotation', 90);
yline(0.158, "--", 'LineWidth', 3);
yline(0.012, "--", 'LineWidth', 3);
text(650, 0.025, "y = 0.158");
text(650, 0.17, "y = 0.012");
plot(520, 0.012, '*', 'Color', 'black', 'LineWidth', 4, 'MarkerSize', 12);
plot(520, 0.158, '*', 'Color', 'black', 'LineWidth', 4, 'MarkerSize', 12);

hold off; % Release the figure

% Add title and labels
title('CoPE and CrPE vs. wavelength [nm]');
xlabel('Wavelength \lambda [nm]');
ylabel('Polarization Efficiency (a.u.)');

% Add legend and grid
legend("CoPE", "CrPE");
grid on;
fontsize(20, "points");

% Save the plot as an image
saveas(gcf, 'CoPE_CrPE_plot.png');

%% Diffraction Efficiency (eta_D) Calculation and Plot

% --- Load and Calculate Data ---
% Read the digitized experimental data (dashed line)
exp_line_data = readmatrix('diffraction_efficiency.csv');
% Read the separate experimental data points (markers)
exp_points_data = readmatrix('diffraction_experiment.csv');

% Calculate the idealized theoretical diffraction efficiency from your data
eta_D_calculated = CrPE ./ (CoPE + CrPE);

% --- Plotting Section ---
figure;
hold on;

% Plot your calculated idealized efficiency (solid black line)
plot(x_tx, eta_D_calculated * 100, 'k-', 'LineWidth', 4, 'DisplayName', 'MATLAB Calculated \eta_D (Idealized)');

% Plot the digitized experimental data from the paper (red dashed line)
plot(exp_line_data(:,1), exp_line_data(:,2) * 100, 'r-.', 'LineWidth', 4, 'DisplayName', 'Simulated Trend (from Fig. 3D)');

% Plot the separate experimental data points (blue square markers)
plot(exp_points_data(:,1), exp_points_data(:,2) * 100, 'bs', 'MarkerSize', 14, 'LineWidth', 3, 'DisplayName', 'Experimental Points (from Fig. 3D)');

xline(500, '--', 'LineWidth', 3);
text(495, 10, 'Wavelength \lambda = 500 nm', 'Rotation', 90);
yline(75, '--', 'LineWidth', 3);
text(600, 77.5, '\eta_D = 75%');
xline(520, '--', 'LineWidth', 3);
text(515, 10, 'Wavelength \lambda = 520 nm', 'Rotation', 90);
xline(550, '--', 'LineWidth', 3);
text(545, 10, 'Wavelength \lambda = 550 nm', 'Rotation', 90);

hold off;

% --- Add title, labels, and formatting ---
title('Diffraction Efficiency \eta_D');
xlabel('Wavelength \lambda [nm]');
ylabel('Diffraction Efficiency, \eta_D [%]');
grid on;
ylim([0 100]); % Set y-axis limits to 0-100%
legend('MATLAB Calculated \eta_D (Idealized)', 'Simulated Trend (from Fig. 3D)', 'Experimental Points (from Fig. 3D)');
fontsize(20, "points");

% Save the plot as an image
saveas(gcf, 'eta_D_comparison_plot_with_points.png');