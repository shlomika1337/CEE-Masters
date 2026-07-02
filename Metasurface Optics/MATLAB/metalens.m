close all; clear; clc;

% MATLAB code to simulate a metalens DGMOE
% Coded by Shlomi Avidan

%% --- Constants and Parameters ---
lambda = 550e-9;                % Wavelength (550 nm)
k = 2 * pi / lambda;            % Wavenumber

% Lens-specific parameters from the paper
f_lens = 100e-6;                % Focal length (100 um)
R_lens = 48e-6;                 % Radius of the metasurface (48 um)

%% --- Define the Metasurface (Source) Grid ---
N_source = 1000;                 % Increased grid size for larger lens
L_source = 100e-6;              % Grid width to contain the 96um lens
x0_vec = linspace(-L_source/2, L_source/2, N_source);
y0_vec = linspace(-L_source/2, L_source/2, N_source);
dx0 = x0_vec(2) - x0_vec(1);
dy0 = y0_vec(2) - y0_vec(1);
[X0, Y0] = meshgrid(x0_vec, y0_vec);

%% --- Define the Field at the Metasurface (z=0) ---
% This is the section that has been modified for the lens
r0 = sqrt(X0.^2 + Y0.^2);

% 1. Incident Wave: Uniform illumination (plane wave)
incident_wave = 1;

% 2. Aperture: Circular aperture with the lens radius
aperture = double(r0 <= R_lens);

% 3. Phase Profile: Hyperboloidal phase profile for a lens
phase_profile_lens = (2*pi/lambda) * (f_lens - sqrt(r0.^2 + f_lens^2));
lens_phase = exp(1j * phase_profile_lens);

% Combine to get the final field at the metasurface plane
U0 = incident_wave .* aperture .* lens_phase;

%% --- Define the Observation Grid (XZ Plane) ---
Nz = 500;   % Adjusted for reasonable calculation time
Nx = 200;   % Use an odd number to have a perfect center pixel
z_vec = linspace(75e-6, 125e-6, Nz); % Observe from 0 to 200um
x_vec = linspace(-10e-6, 10e-6, Nx);  % Observe a smaller transverse region
U_xz = zeros(Nx, Nz);

%% --- Huygens-Fresnel Numerical Integration ---
fprintf('Calculating XZ plane propagation...\n');
for i = 1:Nz
    z = z_vec(i);
    % Only calculate for the central xz-plane slice (y=0)
    for j = 1:Nx
        x = x_vec(j);
        y = 0; % Keep y=0 for the xz-plane
        r01 = sqrt((X0 - x).^2 + (Y0 - y).^2 + z^2);
        integrand = U0 .* (z ./ (r01.^2)) .* (1 + 1j ./ (k.*r01)) .* exp(1j .* k .* r01);
        U_xz(j, i) = sum(integrand(:));
    end
    fprintf('Progress: %d/%d\n', i, Nz);
end
U_xz = U_xz * dx0 * dy0 / (1j * lambda);
I_xz = abs(U_xz).^2;
I_xz_norm = I_xz / max(I_xz(:));
fprintf('\nCalculations complete. Now plotting...\n');
%% Numerical XY plane calculation
fprintf('\nCalculating XY focal plane propagation (this will take longer)...\n');
z_focus = f_lens; % Set observation plane at the focal length
Nx_xy = 200; % Use a smaller grid for faster XY calculation
x_vec_xy = linspace(-2e-6, 2e-6, Nx_xy); % Zoom in to see the focal spot
y_vec_xy = x_vec_xy;
U_xy = zeros(Nx_xy, Nx_xy);

for i = 1:Nx_xy
    for j = 1:Nx_xy
        x = x_vec_xy(i);
        y = y_vec_xy(j);
        r01 = sqrt((X0 - x).^2 + (Y0 - y).^2 + z_focus^2);
        integrand = U0 .* (z_focus ./ (r01.^2)) .* (1 + 1j ./ (k.*r01)) .* exp(1j .* k .* r01);
        U_xy(j, i) = sum(integrand(:));
    end
    fprintf('Progress: %d/%d\n', i, Nx_xy);
end
U_xy = U_xy * dx0 * dy0 / (1j * lambda);
I_xy_norm = abs(U_xy).^2 / max(abs(U_xy(:)).^2);
fprintf('\nXY plane calculation complete.\n');

fprintf('\nDone.\n');
%% --- PLOTTING SECTION (Simplified for Lens) ---

% Plot 1: Intensity in the XZ plane to show focusing
figure('Position', [100, 100, 1000, 400]);
hold on;
imagesc(z_vec * 1e6, x_vec * 1e6, I_xz_norm);
colormap('hot');
colorbar;
axis xy;
axis equal;
% Find the coordinates of the maximum intensity (the focal point)
[~, linear_idx] = max(I_xz_norm(:)); % Find linear index of max value
[row_idx, col_idx] = ind2sub(size(I_xz_norm), linear_idx); % Convert to 2D index
z_at_focus = z_vec(col_idx); % Z-coordinate of focus
x_at_focus = x_vec(row_idx); % X-coordinate of focus

% Plot a marker at the focal point
plot(z_at_focus * 1e6, x_at_focus * 1e6, 'b+', 'MarkerSize', 20, 'LineWidth', 3);
text(z_at_focus * 1e6 + 2, x_at_focus * 1e6 - 2.5, {'Numerically calculated', 'focal point (~0.05\mum shift)'}, 'Color', 'blue');

% Add a text label showing the coordinates
text_str = sprintf('Focus at z = %.1f \\mum', z_at_focus*1e6);
text(z_at_focus*1e6 - 60, x_at_focus*1e6 + 2, text_str, 'Color', 'white', 'FontSize', 14, 'FontWeight', 'bold');

xlim([75 125]);
ylim([-10 10]);
xlabel('Propagation distance z (\mum)');
ylabel('Transverse distance x (\mum)');
title('Simulated Intensity Profile (XZ Plane)');
xline(f_lens * 1e6, 'w--', 'LineWidth', 2, 'Label', 'z=100\mum');
fontsize(16,"points");

% Plot 2: Intensity at the focal plane
% We can get this from our XZ calculation
[~, focal_plane_idx] = min(abs(z_vec - f_lens));
intensity_at_focus = I_xz_norm(:, focal_plane_idx);
FWHM = (sum(intensity_at_focus >= 0.5) * (x_vec(2)-x_vec(1))) * 1e9; % FWHM in nm

figure('Position', [200, 200, 800, 600]);
hold on;
plot(x_vec * 1e6, intensity_at_focus, 'r-', 'LineWidth', 3);
grid on;
% Find the maximum intensity and its index
[max_intensity, max_idx] = max(intensity_at_focus);
% Get the x-position of the maximum intensity
x_at_max = x_vec(max_idx);

% Plot a marker at the actual peak point
plot(x_at_max * 1e6, max_intensity, 'ko', 'MarkerSize', 12, 'LineWidth', 3, 'MarkerFaceColor', 'yellow');

% Plot a vertical dashed line at the peak's x-position
xline(x_at_max * 1e6, '--', 'Color', 'red', 'LineWidth', 2, 'Label', sprintf('Peak at x = %.2f um', x_at_max * 1e6));

title(['Intensity at Focal Plane (z = ' num2str(f_lens*1e6) '\mum)']);
xlabel('Transverse distance x (\mum)');
ylabel('Normalized Intensity');
text(0, 0.6, sprintf('FWHM ≈ %.0f nm', FWHM), 'HorizontalAlignment', 'center', 'FontSize', 14, 'BackgroundColor', 'white');
fontsize(16,"points");


%% --- PLOTTING SECTION - RECREATING THE MULTI-PANEL FIGURE ---
fprintf('\nPlotting multi-panel figure for the metalens...\n');

% Create the main figure
figure('Position', [100, 100, 1200, 800]);
% --- Main XZ Intensity Plot ---
% This is the large background image showing the beam focusing
ax_main = axes('Position', [0.15 0.1 0.75 0.8]);
hold on;
imagesc(z_vec * 1e6, x_vec * 1e6, I_xz_norm);
colormap(ax_main, 'hot');
c = colorbar;
c.Label.String = 'Normalized Intensity';
axis equal;
xline(f_lens * 1e6, 'w--', 'LineWidth', 2, 'Label', 'z=100\mum');
% Plot a marker at the focal point
plot(z_at_focus * 1e6, x_at_focus * 1e6, 'b+', 'MarkerSize', 20, 'LineWidth', 3);
text(z_at_focus * 1e6 + 2, x_at_focus * 1e6 - 2.5, {'Numerically calculated', 'focal point (~0.05\mum shift)'}, 'Color', 'blue');
hold off;
xlabel('Propagation distance z (\mum)');
ylabel('Transverse distance x (\mum)');
set(gca, 'YAxisLocation', 'right');
xlim([75 125]);
ylim([-10 10]);


% --- On-Axis Intensity Plot (Top) ---
% This shows the intensity along the z-axis (at x=0)
ax_top = axes('Position', [0.15 0.72 0.668 0.15]);
x_center_idx = round(Nx / 2);
I_on_axis = I_xz_norm(x_center_idx, :);
plot(z_vec * 1e6, I_on_axis, 'b-', 'LineWidth', 2);
set(ax_top, 'XAxisLocation', 'top');
set(ax_top, 'XTickLabel', []);
title('MATLAB Simulated Metalens Intensity Profile for \lambda = 550 nm', 'FontSize', 18);
ylabel('I (a.u.)');
grid on;
xlim([75 125]);
ylim([0 1.05]);

% --- Transverse Intensity Plot (Left) ---
% This shows the intensity profile at the focal plane
ax_left = axes('Position', [0.035 0.3 0.1 0.40]);
% Find the z-index closest to the design focal length
[~, focal_plane_idx] = min(abs(z_vec - f_lens));
I_transverse_at_focus = I_xz_norm(:, focal_plane_idx);
% Plot Intensity on the X-axis and Position on the Y-axis
plot(I_transverse_at_focus, x_vec * 1e6, 'r-', 'LineWidth', 2);
set(ax_left, 'XDir', 'reverse');
set(ax_left, 'YTickLabel', []);
xlabel('I (a.u.)');
grid on;
ylim([-10 10]); % Match the y-limits of the main plot
xlim([0 1.05]);

% Set the font size for the whole figure
fontsize(24,"points");

fprintf('\nDone.\n');

%% PLOT 2: XY cross-section at the focal plane
figure('Position', [200, 200, 800, 800]);
imagesc(x_vec_xy * 1e6, y_vec_xy * 1e6, I_xy_norm');
axis xy; axis equal; colormap('hot');
c = colorbar; c.Label.String = 'Normalized Intensity';
title(['Intensity at Focal Plane (z = ' num2str(z_focus*1e6) '\mum)']);
xlabel('x (\mum)'); ylabel('y (\mum)');
fontsize(24,"points");
xlim([-2 2]);
ylim([-2 2]);
fprintf('\nDone.\n');

%% Calculate FWHM
[~, focal_plane_idx] = min(abs(z_vec - f_lens));
intensity_at_focus = I_xz_norm(:, focal_plane_idx);
FWHM = (sum(intensity_at_focus >= 0.5) * (x_vec(2)-x_vec(1))) * 1e9; % FWHM in nm
fprintf("FHWM: %e\n", FWHM .* 1e-9);
fprintf("FHWM: ~%.1f-um\n", FWHM .* 1e-3);