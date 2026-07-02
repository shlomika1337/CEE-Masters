close all; clear; clc;

% MATLAB code to process and plot 2D numerical simulation data for a meta-atom.
% Version 3: Added a dedicated column for relative phase in the output table.

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

% 1. Calculate phase values
phi_x_deg = rad2deg(phi_x); 
phi_y_deg = rad2deg(phi_y); 
phi_relative_rad = phi_y - phi_x;
phi_relative_pi = unwrap(phi_relative_rad) / pi; % Unwrap for continuity and normalize by pi

% 2. Calculate efficiencies and ratios
CoPE = abs(0.5 * (Tx + Ty .* exp(1j * phi_relative_rad))).^2;
CrPE = abs(0.5 * (Tx - Ty .* exp(1j * phi_relative_rad))).^2;
diff_efficiency = CrPE ./ (CrPE + CoPE);
amplitude_ratio = Tx ./ Ty;

%% --- Find Key Performance Points and Regions ---

% 1. Find the contour where relative phase is pi (level = 1 in phi_relative_pi)
C_phi_pi = contourc(width_vec, period_vec, phi_relative_pi, [1 1]);

% 2. Find the contour where Tx ~ Ty (ratio is 1)
C_ratio = contourc(width_vec, period_vec, amplitude_ratio, [1 1]);

% 3. Find the point of minimum CoPE (This is the primary "best" point)
[~, idx_min_cope] = min(CoPE(:));
[p_best, w_best] = ind2sub(size(CoPE), idx_min_cope);
best_width = width_vec(w_best);
best_period = period_vec(p_best);

% 4. Find the point of maximum Diffraction Efficiency
[~, idx_max_eta] = max(diff_efficiency(:));
[p_max_eta, w_max_eta] = ind2sub(size(diff_efficiency), idx_max_eta);

% 5. Define and analyze the data point from the paper
paper_width = 120; % nm
paper_period = 200; % nm
% Find the closest indices in our vectors for the paper's data point
[~, w_idx_paper] = min(abs(width_vec - paper_width));
[~, p_idx_paper] = min(abs(period_vec - paper_period));


%% --- Display Performance at Best Design Point (Min CoPE) ---
fprintf('--------------------------------------------------\n');
fprintf('Performance Metrics at the Optimal Design Point (Minimum CoPE):\n');
fprintf('  - Width:        %.1f nm\n', best_width);
fprintf('  - Period:       %.1f nm\n', best_period);
fprintf('  - Tx Amplitude: %.4f\n', Tx(p_best, w_best));
fprintf('  - Ty Amplitude: %.4f\n', Ty(p_best, w_best));
fprintf('  - Tx/Ty ratio:  %.4f\n', amplitude_ratio(p_best, w_best));
fprintf('  - Rel. Phase:   %.4f (pi rad)\n', phi_relative_pi(p_best, w_best));
fprintf('  - CoPE:         %.4f\n', CoPE(p_best, w_best));
fprintf('  - CrPE:         %.4f\n', CrPE(p_best, w_best));
fprintf('  - Diff. Eff.:   %.4f\n', diff_efficiency(p_best, w_best));
fprintf('--------------------------------------------------\n');


%% --- REVISION: Update Table to Include Relative Phase ---

fprintf('\n--- Table of Design Points Where Relative Phase ≈ π ---\n\n');
% Added 'Rel. Phase' column to the header and data rows
fprintf('%-10s | %-10s | %-15s | %-10s | %-10s | %-15s\n', 'Width (nm)', 'Period (nm)', 'Rel. Phase (π)', 'Tx', 'Ty', 'Diff. Eff. (%)');
fprintf('------------------------------------------------------------------------------------------\n');

% Extract and sample points from the pi-phase contour
idx = 1;
all_contour_widths = [];
all_contour_periods = [];
while idx < size(C_phi_pi, 2)
    n_points = C_phi_pi(2, idx);
    x_data = C_phi_pi(1, idx+1 : idx+n_points);
    y_data = C_phi_pi(2, idx+1 : idx+n_points);
    all_contour_widths = [all_contour_widths, x_data];
    all_contour_periods = [all_contour_periods, y_data];
    idx = idx + n_points + 1;
end

% Select a few (e.g., 5) evenly spaced points from the contour for the table
if ~isempty(all_contour_widths)
    num_samples = 5;
    sample_indices = round(linspace(1, length(all_contour_widths), num_samples));
    
    for i = sample_indices
        w = all_contour_widths(i);
        p = all_contour_periods(i);
        % Interpolate the metrics at this specific (w, p)
        tx_val = interp2(width_vec, period_vec, Tx, w, p);
        ty_val = interp2(width_vec, period_vec, Ty, w, p);
        eta_val = interp2(width_vec, period_vec, diff_efficiency, w, p);
        
        % For these points, the phase is by definition 1.00 * pi
        fprintf('%-10.1f | %-11.1f | %-15.2f | %-10.4f | %-10.4f | %-15.2f\n', w, p, 1.00, tx_val, ty_val, eta_val * 100);
    end
end

fprintf('------------------------------------------------------------------------------------------\n');
% Add the data point from the paper to the table
paper_tx = Tx(p_idx_paper, w_idx_paper);
paper_ty = Ty(p_idx_paper, w_idx_paper);
paper_eta = diff_efficiency(p_idx_paper, w_idx_paper);
paper_phi = phi_relative_pi(p_idx_paper, w_idx_paper);
fprintf('%-10.1f | %-11.1f | %-15.2f | %-10.4f | %-10.4f | %-15.2f  <-- Paper''s Point\n', ...
        width_vec(w_idx_paper), period_vec(p_idx_paper), paper_phi, paper_tx, paper_ty, paper_eta * 100);
fprintf('------------------------------------------------------------------------------------------\n');
fprintf('%-10.1f | %-11.1f | %-15.2f | %-10.4f | %-10.4f | %-15.2f  <-- Our Point\n', ...
        width_vec(w_best), period_vec(p_best), phi_relative_pi(p_best, w_best), Tx(p_best, w_best), Ty(p_best, w_best),100*diff_efficiency(p_best, w_best));
fprintf('------------------------------------------------------------------------------------------\n\n');

% Max CrPE
[~, idx_max_crpe] = max(CrPE(:));
[row_max_crpe, col_max_crpe] = ind2sub(size(CrPE), idx_max_crpe);
% Min CoPE
[~, idx_min_cope] = min(CoPE(:));
[row_min_cope, col_min_cope] = ind2sub(size(CoPE), idx_min_cope);
%% --- Plotting Section ---
% Note: Plotting section remains the same as v2

% --- FIGURE 1: Raw Simulation Data ---
figure('Position', [100, 100, 1200, 1000]);
ax_title = axes('Position',[0 0 1 1],'Visible','off');
text(ax_title, 0.5, 0.98, 'Figure 1: Raw Simulation Data at \lambda = 550 nm', 'HorizontalAlignment', 'center', 'FontSize', 20, 'FontWeight', 'bold');
subplot(2, 2, 1); imagesc(width_vec, period_vec, Tx); axis xy; colormap(gca, 'jet'); colorbar; title('Tx Amplitude'); ylabel('Period (nm)'); xlabel('Nano-beam Width (nm)');
subplot(2, 2, 2); imagesc(width_vec, period_vec, Ty); axis xy; colormap(gca, 'jet'); colorbar; title('Ty Amplitude'); ylabel('Period (nm)'); xlabel('Nano-beam Width (nm)');
subplot(2, 2, 3); imagesc(width_vec, period_vec, phi_x_deg); axis xy; colormap(gca, 'jet'); colorbar; title('\phi_x (degrees)'); xlabel('Nano-beam Width (nm)'); ylabel('Period (nm)');
subplot(2, 2, 4); imagesc(width_vec, period_vec, phi_y_deg); axis xy; colormap(gca, 'jet'); colorbar; title('\phi_y (degrees)'); xlabel('Nano-beam Width (nm)'); ylabel('Period (nm)');
fontsize(12, "points");
saveas(gcf, 'simulation_raw_data_v3.png');
%% --- FIGURE 2: Overall Analysis ---
figure('Position', [400, 400, 1800, 1000]);
ax_title4 = axes('Position',[0 0 1 1],'Visible','off');
text(ax_title4, 0.5, 0.98, 'Figure 2: Overall Performance Analysis with Key Design Points', 'HorizontalAlignment', 'center', 'FontSize', 20, 'FontWeight', 'bold');
plot_data_processed = {phi_relative_pi, CrPE, CoPE, amplitude_ratio, diff_efficiency};
plot_titles_processed = {'Relative Phase (\pi rad)', 'CrPE (Deflected Light)', 'CoPE (Undeflected Light)', 'Amplitude Ratio (Tx / Ty)', 'Diffraction Efficiency (\eta_D)'};

for i = 1:5
    subplot(2, 3, i);
    imagesc(width_vec, period_vec, plot_data_processed{i});
    hold on;
        
    hold off;
    axis xy; colormap(gca, 'jet'); colorbar; title(plot_titles_processed{i});
    xlabel('Nano-beam Width (nm)');
    ylabel('Period (nm)');
end
fontsize(12, "points");

%% --- FIGURE 2: Overall Analysis ---
figure('Position', [400, 400, 1800, 1000]);
ax_title4 = axes('Position',[0 0 1 1],'Visible','off');
text(ax_title4, 0.5, 0.98, 'Figure 2: Overall Performance Analysis with Key Design Points', 'HorizontalAlignment', 'center', 'FontSize', 20, 'FontWeight', 'bold');
plot_data_processed = {phi_relative_pi, CrPE, CoPE, amplitude_ratio, diff_efficiency};
plot_titles_processed = {'Relative Phase (\pi rad)', 'CrPE (Deflected Light)', 'CoPE (Undeflected Light)', 'Amplitude Ratio (Tx / Ty)', 'Diffraction Efficiency (\eta_D)'};

for i = 1:5
    subplot(2, 3, i);
    imagesc(width_vec, period_vec, plot_data_processed{i});
    hold on;

    % Plot the primary Optimal Point and the Paper's Point on every subplot
    plot(width_vec(w_idx_paper), period_vec(p_idx_paper), 'mp', 'MarkerSize', 12, 'LineWidth', 2.5, 'MarkerFaceColor', 'm'); 
        
    hold off;
    axis xy; colormap(gca, 'jet'); colorbar; title(plot_titles_processed{i});
    xlabel('Nano-beam Width (nm)');
    ylabel('Period (nm)');
end
fontsize(12, "points");
subplot(2,3,6);
axis off;
hold on;
p2 = plot(NaN, NaN, 'mp', 'MarkerSize', 12, 'LineWidth', 2.5, 'MarkerEdgeColor', 'm');
legend(p2, ...
    {
     sprintf('Paper''s Point (W: %dnm, P: %dnm)', paper_width, paper_period)}, ...
    'Location', 'west');
hold off;


%% --- FIGURE 3: Overall Performance Analysis ---
figure('Position', [400, 400, 1800, 1000]);
ax_title4 = axes('Position',[0 0 1 1],'Visible','off');
text(ax_title4, 0.5, 0.98, 'Figure 2: Overall Performance Analysis with Key Design Points', 'HorizontalAlignment', 'center', 'FontSize', 20, 'FontWeight', 'bold');
plot_data_processed = {phi_relative_pi, CrPE, CoPE, amplitude_ratio, diff_efficiency};
plot_titles_processed = {'Relative Phase (\pi rad)', 'CrPE (Deflected Light)', 'CoPE (Undeflected Light)', 'Amplitude Ratio (Tx / Ty)', 'Diffraction Efficiency (\eta_D)'};

for i = 1:5
    subplot(2, 3, i);
    imagesc(width_vec, period_vec, plot_data_processed{i});
    hold on;
    
    switch i
        case 1 
            plot_contour(C_phi_pi, 'b--', 2);
            text(85, 230, '\phi \approx \pi', 'Color', 'white', 'FontSize', 12, 'FontWeight', 'bold');
        case 2
            plot(width_vec(col_max_crpe), period_vec(row_max_crpe), 'co', 'MarkerSize', 8, 'LineWidth', 2, 'MarkerFaceColor', 'c');
            text(width_vec(col_max_crpe)-2, period_vec(row_max_crpe)+5, 'Max', 'Color', 'cyan', 'FontSize', 12, 'FontWeight', 'bold');
        case 3
            plot(width_vec(col_min_cope), period_vec(row_min_cope), 'co', 'MarkerSize', 8, 'LineWidth', 2, 'MarkerFaceColor', 'c');
            text(width_vec(col_min_cope)+2, period_vec(row_min_cope), 'Min', 'Color', 'cyan', 'FontSize', 12, 'FontWeight', 'bold');            
        case 4
            plot_contour(C_ratio, 'k-', 2);
            text(85, 230, 'Tx \approx Ty', 'Color', 'black', 'FontSize', 12, 'FontWeight', 'bold');
        case 5 
            plot(width_vec(w_max_eta), period_vec(p_max_eta), 'co', 'MarkerSize', 8, 'LineWidth', 2, 'MarkerFaceColor', 'c');
            text(width_vec(w_max_eta)+2, period_vec(p_max_eta), 'Max', 'Color', 'cyan', 'FontSize', 12, 'FontWeight', 'bold');
    end

    % Plot the primary Optimal Point and the Paper's Point on every subplot
    plot(width_vec(w_idx_paper), period_vec(p_idx_paper), 'mp', 'MarkerSize', 12, 'LineWidth', 2.5, 'MarkerFaceColor', 'm'); 
        
    hold off;
    axis xy; colormap(gca, 'jet'); colorbar; title(plot_titles_processed{i});
    xlabel('Nano-beam Width (nm)');
    ylabel('Period (nm)');
end

% Add a legend to the final subplot
subplot(2,3,6);
axis off;
hold on;
p2 = plot(NaN, NaN, 'mp', 'MarkerSize', 12, 'LineWidth', 2.5, 'MarkerEdgeColor', 'm');
p3 = plot(NaN, NaN, 'co', 'MarkerSize', 8, 'LineWidth', 2, 'MarkerFaceColor', 'c', 'MarkerEdgeColor', 'k');
p4 = plot(NaN, NaN, 'b--', 'LineWidth', 2);
p5 = plot(NaN, NaN, 'k-', 'LineWidth', 2);
legend([p2, p3, p4, p5], ...
    {
     sprintf('Paper''s Point (W: %dnm, P: %dnm)', paper_width, paper_period), ...
     'Max Efficiency Point', ...
     'Phase Retardation = \pi', ...
     'Amplitude Ratio = 1'}, ...
    'Location', 'west');
hold off;

fontsize(12, "points");
saveas(gcf, 'simulation_combined_analysis_v3.png');


%% --- FIGURE 3: Overall Performance Analysis ---
figure('Position', [400, 400, 1800, 1000]);
ax_title4 = axes('Position',[0 0 1 1],'Visible','off');
text(ax_title4, 0.5, 0.98, 'Figure 2: Overall Performance Analysis with Key Design Points', 'HorizontalAlignment', 'center', 'FontSize', 20, 'FontWeight', 'bold');
plot_data_processed = {phi_relative_pi, CrPE, CoPE, amplitude_ratio, diff_efficiency};
plot_titles_processed = {'Relative Phase (\pi rad)', 'CrPE (Deflected Light)', 'CoPE (Undeflected Light)', 'Amplitude Ratio (Tx / Ty)', 'Diffraction Efficiency (\eta_D)'};

for i = 1:5
    subplot(2, 3, i);
    imagesc(width_vec, period_vec, plot_data_processed{i});
    hold on;
    
    switch i
        case 1 
            plot_contour(C_phi_pi, 'b--', 2);
            text(85, 230, '\phi \approx \pi', 'Color', 'white', 'FontSize', 12, 'FontWeight', 'bold');
        case 2
            plot(width_vec(col_max_crpe), period_vec(row_max_crpe), 'co', 'MarkerSize', 8, 'LineWidth', 2, 'MarkerFaceColor', 'c');
            text(width_vec(col_max_crpe)-2, period_vec(row_max_crpe)+5, 'Max', 'Color', 'cyan', 'FontSize', 12, 'FontWeight', 'bold');
        case 3
            plot(width_vec(col_min_cope), period_vec(row_min_cope), 'co', 'MarkerSize', 8, 'LineWidth', 2, 'MarkerFaceColor', 'c');
            text(width_vec(col_min_cope)+2, period_vec(row_min_cope), 'Min', 'Color', 'cyan', 'FontSize', 12, 'FontWeight', 'bold');            
        case 4
            plot_contour(C_ratio, 'k-', 2);
            text(85, 230, 'Tx \approx Ty', 'Color', 'black', 'FontSize', 12, 'FontWeight', 'bold');
        case 5 
            plot(width_vec(w_max_eta), period_vec(p_max_eta), 'co', 'MarkerSize', 8, 'LineWidth', 2, 'MarkerFaceColor', 'c');
            text(width_vec(w_max_eta)+2, period_vec(p_max_eta), 'Max', 'Color', 'cyan', 'FontSize', 12, 'FontWeight', 'bold');
    end

    % Plot the primary Optimal Point and the Paper's Point on every subplot
    plot(best_width, best_period, 'yp', 'MarkerSize', 16, 'LineWidth', 2, 'MarkerFaceColor', 'yellow');
    plot(width_vec(w_idx_paper), period_vec(p_idx_paper), 'mp', 'MarkerSize', 12, 'LineWidth', 2.5, 'MarkerFaceColor', 'm'); 
        
    hold off;
    axis xy; colormap(gca, 'jet'); colorbar; title(plot_titles_processed{i});
    xlabel('Nano-beam Width (nm)');
    ylabel('Period (nm)');
end

% Add a legend to the final subplot
subplot(2,3,6);
axis off;
hold on;
p1 = plot(NaN, NaN, 'yp', 'MarkerSize', 16, 'LineWidth', 2, 'MarkerFaceColor', 'yellow', 'MarkerEdgeColor', 'k');
p2 = plot(NaN, NaN, 'mp', 'MarkerSize', 12, 'LineWidth', 2.5, 'MarkerEdgeColor', 'm');
p3 = plot(NaN, NaN, 'co', 'MarkerSize', 8, 'LineWidth', 2, 'MarkerFaceColor', 'c', 'MarkerEdgeColor', 'k');
p4 = plot(NaN, NaN, 'b--', 'LineWidth', 2);
p5 = plot(NaN, NaN, 'k-', 'LineWidth', 2);
legend([p1, p2, p3, p4, p5], ...
    {sprintf('Optimal Point (W: %.1fnm, P: %.1fnm)', best_width, best_period), ...
     sprintf('Paper''s Point (W: %dnm, P: %dnm)', paper_width, paper_period), ...
     'Max Efficiency Point', ...
     'Phase Retardation = \pi', ...
     'Amplitude Ratio = 1'}, ...
    'Location', 'west');
hold off;

fontsize(12, "points");
saveas(gcf, 'simulation_combined_analysis_v3.png');

%% --- FOR PDF ---
figure('Position', [400, 400, 1800, 1000]);
ax_title4 = axes('Position',[0 0 1 1],'Visible','off');
text(ax_title4, 0.5, 0.98, 'Figure 2: Overall Performance Analysis with Key Design Points', 'HorizontalAlignment', 'center', 'FontSize', 20, 'FontWeight', 'bold');
plot_data_processed = {phi_relative_pi, CrPE, CoPE, amplitude_ratio, diff_efficiency};
plot_titles_processed = {'Relative Phase (\pi rad)', 'CrPE (Deflected Light)', 'CoPE (Undeflected Light)', 'Amplitude Ratio (Tx / Ty)', 'Diffraction Efficiency (\eta_D)'};

for i = 1:5
    subplot(2, 3, i);
    imagesc(width_vec, period_vec, plot_data_processed{i});
    hold on;
    
    switch i
        case 1 
            plot_contour(C_phi_pi, 'b--', 2);
            text(85, 230, '\phi \approx \pi', 'Color', 'white', 'FontSize', 12, 'FontWeight', 'bold');
        case 2
            plot(width_vec(col_max_crpe), period_vec(row_max_crpe), 'co', 'MarkerSize', 8, 'LineWidth', 2, 'MarkerFaceColor', 'c');
            text(width_vec(col_max_crpe)-2, period_vec(row_max_crpe)+5, 'Max', 'Color', 'cyan', 'FontSize', 12, 'FontWeight', 'bold');
        case 3
            plot(width_vec(col_min_cope), period_vec(row_min_cope), 'co', 'MarkerSize', 8, 'LineWidth', 2, 'MarkerFaceColor', 'c');
            text(width_vec(col_min_cope)+2, period_vec(row_min_cope), 'Min', 'Color', 'cyan', 'FontSize', 12, 'FontWeight', 'bold');            
        case 4
            plot_contour(C_ratio, 'k-', 2);
            text(85, 230, 'Tx \approx Ty', 'Color', 'black', 'FontSize', 12, 'FontWeight', 'bold');
        case 5 
            plot(width_vec(w_max_eta), period_vec(p_max_eta), 'co', 'MarkerSize', 8, 'LineWidth', 2, 'MarkerFaceColor', 'c');
            text(width_vec(w_max_eta)+2, period_vec(p_max_eta), 'Max', 'Color', 'cyan', 'FontSize', 12, 'FontWeight', 'bold');
    end

    % Plot the primary Optimal Point and the Paper's Point on every subplot
    plot(best_width, best_period, 'yp', 'MarkerSize', 16, 'LineWidth', 2, 'MarkerFaceColor', 'yellow');
    plot(width_vec(w_idx_paper), period_vec(p_idx_paper), 'mp', 'MarkerSize', 12, 'LineWidth', 2.5, 'MarkerFaceColor', 'm'); 
        
    hold off;
    axis xy; colormap(gca, 'jet'); colorbar; title(plot_titles_processed{i});
    xlabel('Nano-beam Width (nm)');
    ylabel('Period (nm)');
end

% Add a legend to the final subplot
subplot(2,3,6);
axis off;
hold on;
p1 = plot(NaN, NaN, 'yp', 'MarkerSize', 16, 'LineWidth', 2, 'MarkerFaceColor', 'yellow', 'MarkerEdgeColor', 'k');
p2 = plot(NaN, NaN, 'mp', 'MarkerSize', 12, 'LineWidth', 2.5, 'MarkerEdgeColor', 'm');
p3 = plot(NaN, NaN, 'co', 'MarkerSize', 8, 'LineWidth', 2, 'MarkerFaceColor', 'c', 'MarkerEdgeColor', 'k');
p4 = plot(NaN, NaN, 'b--', 'LineWidth', 2);
p5 = plot(NaN, NaN, 'k-', 'LineWidth', 2);
legend([p1, p2, p3, p4, p5], ...
    {sprintf('Optimal Point (W: %.1fnm, P: %.1fnm)', best_width, best_period), ...
     sprintf('Paper''s Point (W: %dnm, P: %dnm)', paper_width, paper_period), ...
     'Max Efficiency Point', ...
     'Phase Retardation = \pi', ...
     'Amplitude Ratio = 1'}, ...
    'Location', 'west');
hold off;

fontsize(16, "points");

%% 2. Define the condition to find values close to pi.
% We are looking for where the normalized phase is within 3% of 1 (pi).
condition = abs(phi_relative_pi - 1) < 0.03;

% 3. Find the linear indices of the elements that satisfy the condition.
linear_indices = find(condition);

% 4. Check if any values were found.
if isempty(linear_indices)
    disp('No data points found that satisfy the condition.');
    disp('You may want to relax the tolerance (e.g., set it to 0.05).');
else
    % 5. Retrieve all corresponding values using the linear indices.
    [row_indices, col_indices] = find(condition);
    
    found_periods = period_vec(row_indices)'; % Transpose to make it a column vector
    found_widths = width_vec(col_indices)';   % Transpose to make it a column vector
    
    found_phi_rad = phi_relative_rad(linear_indices);
    found_Tx = Tx(linear_indices);
    found_Ty = Ty(linear_indices);

    % 6. Calculate the additional requested metrics.
    amplitude_ratio = found_Tx ./ found_Ty;
    
    % Efficiencies are calculated based on the formulas from your term paper.
    % CoPE = |1/2 * (tx + ty*exp(i*phi))|^2
    CoPE = 0.25 * abs(found_Tx + found_Ty .* exp(1i * found_phi_rad)).^2;
    
    % CrPE = |1/2 * (tx - ty*exp(i*phi))|^2
    CrPE = 0.25 * abs(found_Tx - found_Ty .* exp(1i * found_phi_rad)).^2;
    
    % Diffraction Efficiency is the ratio of deflected power to total transmitted power
    diff_efficiency = CrPE ./ (CoPE + CrPE);

    % 7. Create and display a table with all the results.
    results_table = table(found_widths, found_periods, found_phi_rad./pi, found_Tx, found_Ty, ...
        amplitude_ratio, CoPE, CrPE, diff_efficiency * 100, ...
        'VariableNames', {'Width [nm]', 'Period [nm]', 'Phase difference [rad/pi]' 'Tx', 'Ty', ...
        'Amplitude Ratio (Tx/Ty)', 'CoPE', 'CrPE', 'Diffraction Efficiency [%]'});

    fprintf('Found %d data points where the relative phase is approximately pi:\n\n', height(results_table));
    disp(results_table);
end


%% --- Helper function to parse and plot contour data ---
function plot_contour(C, line_style, line_width)
    if isempty(C)
        return; % Do nothing if the contour matrix is empty
    end
    idx = 1;
    while idx < size(C, 2)
        n_points = C(2, idx);
        x_data = C(1, idx+1 : idx+n_points);
        y_data = C(2, idx+1 : idx+n_points);
        plot(x_data, y_data, line_style, 'LineWidth', line_width);
        idx = idx + n_points + 1;
    end
end