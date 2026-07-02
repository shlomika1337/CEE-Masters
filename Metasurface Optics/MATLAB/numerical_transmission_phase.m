close all; clear; clc;

% MATLAB code to process and plot 2D numerical simulation data for a meta-atom.
% This script analyzes performance as a function of nano-beam width and period
% at a fixed wavelength of 550 nm.

%% --- Load and Prepare Data ---

% Load the simulation data from the .mat files
load('Tx_phix_data.mat'); % Should contain 2D matrices: 'Tx', 'phi_x'
load('Ty_phiy_data.mat'); % Should contain 2D matrices: 'Ty', 'phi_y'

% --- Define the Geometric Axes ---
width_vec = 80:2:160;   % Width: 80nm to 160nm with 2nm step
period_vec = 160:2:240;  % Period: 160nm to 240nm with 2nm step

% --- Validate Matrix Dimensions ---
[num_periods, num_widths] = size(Tx);
if num_periods ~= length(period_vec) || num_widths ~= length(width_vec)
    warning('Matrix dimensions do not match the defined axis vectors. Check the ranges.');
end


%% --- Perform Calculations (on 2D matrices) ---

% 1. Calculate phase values in degrees
phi_x_deg = rad2deg(phi_x);
phi_y_deg = rad2deg(phi_y);
phi_relative_rad = phi_y - phi_x;
phi_relative_rad_norm = unwrap(phi_relative_rad) / pi;
phi_relative_deg = rad2deg(unwrap(phi_relative_rad)); % Unwrap for continuity

% 2. Calculate efficiencies and ratios
CoPE = abs(0.5 * (Tx + Ty .* exp(1j * phi_relative_rad))).^2;
CrPE = abs(0.5 * (Tx - Ty .* exp(1j * phi_relative_rad))).^2;
diff_efficiency = CrPE ./ (CrPE + CoPE);
amplitude_ratio = Tx ./ Ty;

%% --- Find Key Performance Points and Regions ---

% 1. Find the contour where phase is ~180 degrees
C_phi = contourc(width_vec, period_vec, phi_relative_rad_norm, [1 1]);

% 2. Find the contour where Tx ~ Ty (ratio is 1)
C_ratio = contourc(width_vec, period_vec, amplitude_ratio, [1 1]);

% 3. Find the point of minimum CoPE (This is now our single "best" point)
[~, idx_min_cope] = min(CoPE(:));
[p_best, w_best] = ind2sub(size(CoPE), idx_min_cope);
best_width = width_vec(w_best);
best_period = period_vec(p_best);

% 4. Find the point of maximum CrPE
[~, idx_max_crpe] = max(CrPE(:));
[p_max_crpe, w_max_crpe] = ind2sub(size(CrPE), idx_max_crpe);

% 5. Find the point of maximum Diffraction Efficiency
[~, idx_max_eta] = max(diff_efficiency(:));
[p_max_eta, w_max_eta] = ind2sub(size(diff_efficiency), idx_max_eta);


%% --- Display Performance at Best Design Point (Min CoPE) ---
fprintf('--------------------------------------------------\n');
fprintf('Performance Metrics at the Optimal Design Point (Minimum CoPE):\n');
fprintf('  - Width:        %.1f nm\n', best_width);
fprintf('  - Period:       %.1f nm\n', best_period);
fprintf('  - Tx Amplitude: %.4f\n', Tx(p_best, w_best));
fprintf('  - Ty Amplitude: %.4f\n', Ty(p_best, w_best));
fprintf('  - Tx/Ty ratio:  %.4f\n', Tx(p_best, w_best)/Ty(p_best, w_best));
fprintf('  - Rel. Phase:   %.4f [rad/\pi]\n', phi_relative_rad_norm(p_best, w_best));
fprintf('  - CoPE:         %.4f\n', CoPE(p_best, w_best));
fprintf('  - CrPE:         %.4f\n', CrPE(p_best, w_best));
fprintf('  - Diff. Eff.:   %.4f\n', diff_efficiency(p_best, w_best));
fprintf('--------------------------------------------------\n');


%% --- Calculate Distance Between Key Points ---
width_at_max_eta = width_vec(w_max_eta);
period_at_max_eta = period_vec(p_max_eta);
width_at_min_cope = width_vec(w_best); % Updated to use 'best' variables
period_at_min_cope = period_vec(p_best);

distance_eta_cope = sqrt((width_at_max_eta - width_at_min_cope)^2 + (period_at_max_eta - period_at_min_cope)^2);
fprintf('Geometric distance between Max η_D and Min CoPE points: %.2f nm\n', distance_eta_cope);


%% --- Plotting Section ---

% --- FIGURE 1: Raw Simulation Data ---
figure('Position', [100, 100, 1200, 1000]);
ax_title = axes('Position',[0 0 1 1],'Visible','off');
text(ax_title, 0.5, 0.98, 'Figure 1: Raw Simulation Data at \lambda = 550 nm', 'HorizontalAlignment', 'center', 'FontSize', 20, 'FontWeight', 'bold');
subplot(2, 2, 1); imagesc(width_vec, period_vec, Tx); axis xy; colormap(gca, 'jet'); colorbar; title('Tx Amplitude'); ylabel('Period (nm)'); xlabel('Nano-beam Width (nm)');
subplot(2, 2, 2); imagesc(width_vec, period_vec, Ty); axis xy; colormap(gca, 'jet'); colorbar; title('Ty Amplitude'); ylabel('Period (nm)'); xlabel('Nano-beam Width (nm)');
subplot(2, 2, 3); imagesc(width_vec, period_vec, phi_x_deg); axis xy; colormap(gca, 'jet'); colorbar; title('\phi_x (degrees)'); xlabel('Nano-beam Width (nm)'); ylabel('Period (nm)');
subplot(2, 2, 4); imagesc(width_vec, period_vec, phi_y_deg); axis xy; colormap(gca, 'jet'); colorbar; title('\phi_y (degrees)'); xlabel('Nano-beam Width (nm)'); ylabel('Period (nm)');
fontsize(12, "points");
saveas(gcf, 'simulation_raw_data.png');


% --- FIGURE 2: Processed Performance Metrics with Individual Indicators ---
figure('Position', [200, 200, 1800, 1000]);
ax_title2 = axes('Position',[0 0 1 1],'Visible','off');
text(ax_title2, 0.5, 0.98, 'Figure 2: Processed Performance Metrics at \lambda = 550 nm (Individual Indicators)', 'HorizontalAlignment', 'center', 'FontSize', 20, 'FontWeight', 'bold');
subplot(2, 3, 1); imagesc(width_vec, period_vec, phi_relative_rad_norm); hold on; plot(120, 200, 'm*', 'MarkerSize', 7, 'LineWidth', 4); plot_contour(C_phi, 'w--', 2); text(85, 230, '\phi \approx \pi [rad]', 'Color', 'white', 'FontSize', 12, 'FontWeight', 'bold'); hold off; axis xy; colormap(gca, 'jet'); colorbar; title('Relative Phase (degrees)'); ylabel('Period (nm)'); xlabel('Nano-beam Width (nm)'); 
subplot(2, 3, 2); imagesc(width_vec, period_vec, CrPE); hold on; plot(120, 200, 'm*', 'MarkerSize', 7, 'LineWidth', 4); plot(width_vec(w_max_crpe), period_vec(p_max_crpe), 'r*', 'MarkerSize', 7, 'LineWidth', 2); text(width_vec(w_max_crpe)+2, period_vec(p_max_crpe), 'Max', 'Color', 'red', 'FontSize', 12, 'FontWeight', 'bold'); hold off; axis xy; colormap(gca, 'jet'); colorbar; title('CrPE (Deflected Light)'); ylabel('Period (nm)'); xlabel('Nano-beam Width (nm)');
subplot(2, 3, 3); imagesc(width_vec, period_vec, CoPE); hold on; plot(120, 200, 'm*', 'MarkerSize', 7, 'LineWidth', 4); plot(width_vec(w_best), period_vec(p_best), 'go', 'MarkerSize', 7, 'LineWidth', 2, 'MarkerFaceColor', 'g'); text(width_vec(w_best)+2, period_vec(p_best), 'Min', 'Color', 'green', 'FontSize', 12, 'FontWeight', 'bold'); hold off; axis xy; colormap(gca, 'jet'); colorbar; title('CoPE (Undeflected Light)'); ylabel('Period (nm)'); xlabel('Nano-beam Width (nm)');
subplot(2, 3, 4); imagesc(width_vec, period_vec, amplitude_ratio); hold on; plot(120, 200, 'm*', 'MarkerSize', 7, 'LineWidth', 4); plot_contour(C_ratio, 'k-', 2); text(85, 230, 'Tx \approx Ty', 'Color', 'black', 'FontSize', 12, 'FontWeight', 'bold'); hold off; axis xy; colormap(gca, 'jet'); colorbar; title('Amplitude Ratio (Tx / Ty)'); xlabel('Nano-beam Width (nm)'); ylabel('Period (nm)');
subplot(2, 3, 5); imagesc(width_vec, period_vec, diff_efficiency); hold on; plot(120, 200, 'm*', 'MarkerSize', 7, 'LineWidth', 4); plot(width_vec(w_max_eta), period_vec(p_max_eta), 'co', 'MarkerSize', 7, 'LineWidth', 2, 'MarkerFaceColor', 'c'); text(width_vec(w_max_eta)+2, period_vec(p_max_eta), 'Max', 'Color', 'cyan', 'FontSize', 12, 'FontWeight', 'bold'); hold off; axis xy; colormap(gca, 'jet'); colorbar; title('Diffraction Efficiency (\eta_D)'); xlabel('Nano-beam Width (nm)'); ylabel('Period (nm)');
fontsize(12, "points");
saveas(gcf, 'simulation_processed_data.png');


% --- FIGURE 3: Combined Indicator Analysis ---
figure('Position', [400, 400, 1800, 1000]);
ax_title4 = axes('Position',[0 0 1 1],'Visible','off');
text(ax_title4, 0.5, 0.98, 'Figure 3: Overall Performance Analysis with Optimal Point', 'HorizontalAlignment', 'center', 'FontSize', 20, 'FontWeight', 'bold');
plot_data_processed = {phi_relative_rad_norm, CrPE, CoPE, amplitude_ratio, diff_efficiency};
plot_titles_processed = {'Relative Phase (rad/\pi)', 'CrPE (Deflected Light)', 'CoPE (Undeflected Light)', 'Amplitude Ratio (Tx / Ty)', 'Diffraction Efficiency (\eta_D)'};
colormaps_processed = {'jet', 'jet', 'jet', 'jet', 'jet'};
for i = 1:5
    subplot(2, 3, i);
    imagesc(width_vec, period_vec, plot_data_processed{i});
    hold on;
    
    % Plot the individual optimal indicator for the current subplot
    switch i
        case 1 % Relative Phase plot
            plot_contour(C_phi, 'w--', 2);
            text(85, 230, '\phi \approx \pi [rad]', 'Color', 'white', 'FontSize', 12, 'FontWeight', 'bold');
        case 2 % CrPE plot
            plot(width_vec(w_max_crpe), period_vec(p_max_crpe), 'r*', 'MarkerSize', 7, 'LineWidth', 2);
            text(width_vec(w_max_crpe)+2, period_vec(p_max_crpe), 'Max', 'Color', 'red', 'FontSize', 12, 'FontWeight', 'bold');
        case 3 % CoPE plot
            plot(width_vec(w_best), period_vec(p_best), 'go', 'MarkerSize', 7, 'LineWidth', 2, 'MarkerFaceColor', 'g');
            text(width_vec(w_best)+2, period_vec(p_best), 'Min', 'Color', 'green', 'FontSize', 12, 'FontWeight', 'bold');
        case 4 % Amplitude Ratio plot
            plot_contour(C_ratio, 'k-', 2);
            text(85, 230, 'Tx \approx Ty', 'Color', 'black', 'FontSize', 12, 'FontWeight', 'bold');
        case 5 % Diffraction Efficiency plot
            plot(width_vec(w_max_eta), period_vec(p_max_eta), 'co', 'MarkerSize', 7, 'LineWidth', 2, 'MarkerFaceColor', 'c');
            text(width_vec(w_max_eta)+2, period_vec(p_max_eta), 'Max', 'Color', 'cyan', 'FontSize', 12, 'FontWeight', 'bold');
    end
    % Plot the "Optimal Point" (Min CoPE) on every subplot LAST so it's on top
    plot(best_width, best_period, 'yp', 'MarkerSize', 16, 'LineWidth', 2, 'MarkerFaceColor', 'yellow');
    
    hold off;
    axis xy; colormap(gca, colormaps_processed{i}); colorbar; title(plot_titles_processed{i});
    xlabel('Nano-beam Width (nm)');
    ylabel('Period (nm)');
end
fontsize(12, "points");
saveas(gcf, 'simulation_all_indicators_with_optimal_point.png');


% --- FIGURE 4: Side-by-Side Comparison ---
figure('Position', [500, 100, 1200, 1600]);
ax_title5 = axes('Position',[0 0 1 1],'Visible','off');
text(ax_title5, 0.5, 0.98, 'Figure 4: Side-by-Side Comparison of Optimal Indicators', 'HorizontalAlignment', 'center', 'FontSize', 20, 'FontWeight', 'bold');
for i = 1:5
    % Left subplot (Individual Indicator)
    subplot(5, 2, 2*i - 1);
    imagesc(width_vec, period_vec, plot_data_processed{i});
    hold on;
    switch i
        case 1; plot_contour(C_phi, 'w--', 2); text(85, 230, '\phi \approx \pi [rad]', 'Color', 'white', 'FontSize', 12, 'FontWeight', 'bold');
        case 2; plot(width_vec(w_max_crpe), period_vec(p_max_crpe), 'r*', 'MarkerSize', 7, 'LineWidth', 2); text(width_vec(w_max_crpe)+2, period_vec(p_max_crpe), 'Max', 'Color', 'red', 'FontSize', 12, 'FontWeight', 'bold');
        case 3; plot(width_vec(w_best), period_vec(p_best), 'go', 'MarkerSize', 7, 'LineWidth', 2, 'MarkerFaceColor', 'g'); text(width_vec(w_best)+2, period_vec(p_best), 'Min', 'Color', 'green', 'FontSize', 12, 'FontWeight', 'bold');
        case 4; plot_contour(C_ratio, 'k-', 2); text(85, 230, 'Tx \approx Ty', 'Color', 'black', 'FontSize', 12, 'FontWeight', 'bold');
        case 5; plot(width_vec(w_max_eta), period_vec(p_max_eta), 'co', 'MarkerSize', 7, 'LineWidth', 2, 'MarkerFaceColor', 'c'); text(width_vec(w_max_eta)+2, period_vec(p_max_eta), 'Max', 'Color', 'cyan', 'FontSize', 12, 'FontWeight', 'bold');
    end
    hold off;
    axis xy; colormap(gca, colormaps_processed{i}); colorbar; title(plot_titles_processed{i});
    xlabel('Nano-beam Width (nm)'); ylabel('Period (nm)');
    
    % Right subplot (Combined Indicators)
    subplot(5, 2, 2*i);
    imagesc(width_vec, period_vec, plot_data_processed{i});
    hold on;
    % Plot ALL indicators on the right side
    plot_contour(C_phi, 'w--', 2);
    plot_contour(C_ratio, 'k-', 2);
    plot(width_vec(w_max_crpe), period_vec(p_max_crpe), 'r*', 'MarkerSize', 7, 'LineWidth', 2);
    plot(width_vec(w_best), period_vec(p_best), 'go', 'MarkerSize', 7, 'LineWidth', 2, 'MarkerFaceColor', 'g');
    plot(width_vec(w_max_eta), period_vec(p_max_eta), 'co', 'MarkerSize', 7, 'LineWidth', 2, 'MarkerFaceColor', 'c');
    plot(best_width, best_period, 'mp', 'MarkerSize', 16, 'LineWidth', 2, 'MarkerFaceColor', 'magenta');
    hold off;
    axis xy; colormap(gca, colormaps_processed{i}); colorbar; title([plot_titles_processed{i}]);
    xlabel('Nano-beam Width (nm)'); ylabel('Period (nm)');
end
fontsize(12, "points");
saveas(gcf, 'simulation_side_by_side.png');


% --- Helper function to parse and plot contour data ---
function plot_contour(C, line_style, line_width)
    idx = 1;
    while idx < size(C, 2)
        %level = C(1, idx); % Level is not needed for plotting
        n_points = C(2, idx);
        x_data = C(1, idx+1 : idx+n_points);
        y_data = C(2, idx+1 : idx+n_points);
        plot(x_data, y_data, line_style, 'LineWidth', line_width);
        idx = idx + n_points + 1;
    end
end
