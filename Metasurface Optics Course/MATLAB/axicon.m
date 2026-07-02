close all; clear; clc;

% MATLAB code to simulate a Bessel beam from an axicon DGMOE
% Coded by Shlomi Avidan

%% --- Constants and Parameters ---
lambda = 550e-9;                % Wavelength (550 nm)
k = 2 * pi / lambda;            % Wavenumber
Lambda_axicon = 3.2e-6;         % Axicon period (3.2 um)
kr = 2 * pi / Lambda_axicon;    % Radial wavevector from axicon phase
R = 32e-6;                      % Radius of the metasurface (32 um)

% --- Define the Metasurface (Source) Grid ---
N_source = 400; 
L_source = 80e-6; 
x0_vec = linspace(-L_source/2, L_source/2, N_source);
y0_vec = linspace(-L_source/2, L_source/2, N_source);
dx0 = x0_vec(2) - x0_vec(1);
dy0 = y0_vec(2) - y0_vec(1);
[X0, Y0] = meshgrid(x0_vec, y0_vec);

% --- Define the Field at the Metasurface (z=0) ---
w0 = R / 2; 
r0 = sqrt(X0.^2 + Y0.^2);
incident_wave = exp(-r0.^2 / w0^2); 
aperture = double(r0 <= R); 
axicon_phase = exp(1j * (-kr * r0)); 
U0 = incident_wave .* aperture .* axicon_phase; 

% --- Define the Observation Grid (XZ Plane) ---
Nz = 2000; 
Nx = 300; 
z_vec = linspace(0.1e-6, 200e-6, Nz);
x_vec = linspace(-15e-6, 15e-6, Nx);
U_xz = zeros(Nx, Nz); 

%% --- Huygens-Fresnel Numerical Integration ---
fprintf('Calculating XZ plane propagation...\n');
for i = 1:Nz
    for j = 1:Nx
        z = z_vec(i);
        x = x_vec(j);
        y = 0;
        r01 = sqrt((X0 - x).^2 + (Y0 - y).^2 + z^2);
        integrand = U0 .* (z ./ (r01.^2)) .* (1 + 1j ./ (k.*r01)) .* exp(1j .* k .* r01);
        U_xz(j, i) = sum(integrand(:));
    end
    fprintf('Calculating XZ plane propagation... %d/%d\n', i, Nz);
end
U_xz = U_xz * dx0 * dy0 / (1j * lambda);
I_xz = abs(U_xz).^2;
I_xz_norm = I_xz / max(I_xz(:));

% --- Calculation for Transverse Plane at z=50um ---
z_focus = 50e-6; 
U_xy = zeros(Nx, Nx); 
fprintf('\nCalculating transverse XY plane at z=50um...\n');
for i = 1:Nx 
    for j = 1:Nx 
        x = x_vec(i);
        y = x_vec(j);
        r01 = sqrt((X0 - x).^2 + (Y0 - y).^2 + z_focus^2);
        integrand = U0 .* (z_focus ./ (r01.^2)) .* (1 + 1j ./ (k.*r01)) .* exp(1j .* k .* r01);
        U_xy(i, j) = sum(integrand(:));
    end
    fprintf('Calculating XY plane propagation... %d/%d\n', i, Nx);
end
U_xy = U_xy * dx0 * dy0 / (1j * lambda);
I_xy = abs(U_xy).^2;
I_xy_norm = I_xy / max(I_xy(:));

fprintf('\nCalculations complete. Now plotting...\n');

%% --- PLOTTING SECTION - FIG 1E RECREATION ---
figure('Position', [100, 100, 2560, 1440]);

% --- Main XZ Intensity Plot ---
ax_main = axes('Position', [0.1 0.01 0.88 0.85]);
imagesc(z_vec * 1e6, x_vec * 1e6, I_xz_norm);
colormap(ax_main, 'hot');
colorbar;
axis equal;
xlabel('Propagation distance z (\mum)');
ylabel('Transverse distance x (\mum)');
set(gca, 'YAxisLocation', 'right');
ylim([-15 15]);
xlim([0 200]);


% --- On-Axis Intensity Plot (Top) ---
ax_top = axes('Position', [0.1 0.55 0.805 0.15]);
x_center_idx = round(Nx / 2);
I_on_axis = I_xz_norm(x_center_idx, :);
plot(z_vec * 1e6, I_on_axis, 'b-', 'LineWidth', 1.5);
set(ax_top, 'XAxisLocation', 'top');
set(ax_top, 'XTickLabel', []); 
title('Simulated (MATLAB) Bessel Beam Intensity Profile for \lambda = 550 nm');
ylabel('I (a.u.)');
grid on;
xlim([0 200]);
ylim([0 1.01])

% --- Transverse Intensity Plot (Left) ---
ax_left = axes('Position', [0.04 0.33 0.052 0.21]);
I_transverse_x_cut = I_xy_norm(:, x_center_idx);

% MODIFICATION: Plot Intensity on the X-axis and Position on the Y-axis
plot(I_transverse_x_cut, x_vec * 1e6, 'r-', 'LineWidth', 1.5);

% Flip the plot direction to match the article
set(ax_left, 'XDir', 'reverse'); 
set(ax_left, 'YTickLabel', []);
xlabel('I (a.u.)');
grid on;
ylim([-15 15]);
xlim([0 1.01]);

fontsize(16,"points");

%% --- PLOTTING SECTION - XY PLANE ---
figure;

% Plot the 2D intensity profile as a colormap
imagesc(x_vec * 1e6, x_vec * 1e6, I_xy_norm'); 
axis xy; % Put origin at bottom-left
axis equal; % Ensure correct aspect ratio (no distortion)
colormap('hot');
c = colorbar;
c.Label.String = 'Normalized Intensity';
title(['XY-plane cross-section of intensity profile at z = ' num2str(z_focus*1e6) ' \mum']);
xlabel('x (\mum)');
ylabel('y (\mum)');
xlim([-7.5 7.5]);
ylim([-7.5 7.5]);
fontsize(21,"points");
fprintf('\nDone.\n');

%% --- PLOTTING SECTION - XZ PLANE ---
figure;

% Find the index for y = 0
y_center_idx = round(Nx / 2);

% Extract the intensity data along the x-axis at y = 0
intensity_cross_section = I_xy_norm(:, y_center_idx);

% Create the line plot
hold on;
plot(x_vec * 1e6, intensity_cross_section, 'r-', 'LineWidth', 5);
yline(0.154, '--', 'LineWidth', 4);
text(1.65, 0.7, 'x = 1.2\mum', 'Rotation', 90);
text(3.35, 0.7, 'x = 2.85\mum', 'Rotation', 90);
text(5, 0.7, 'x = 4.45\mum', 'Rotation', 90);
text(5, 0.185, 'f(x) = 0.154');
title('Intensity profile at xy-plane cross-section z=50\mum, for y=0');
xlabel('Transverse distance x [\mum]');
ylabel('Normalized intensity (a.u.)');
xline(1.2, '--', 'LineWidth', 4);
xline(2.85, '--', 'LineWidth', 4);
xline(4.45, '--', 'LineWidth', 4);
grid on;
xlim([-15 15]);
ylim([0 1.05]); % Set y-limit for better visualization
fontsize(21,"points");
fprintf('\nDone.\n');

%% --- Integration for inner spot and outer ring ---

% Find the indices of the data points within the first range
idx1 = x_vec >= 0 & x_vec <= 1.2e-6;
idx2 = x_vec >= 1.2e-6 & x_vec <= 2.85e-6;
idx3 = x_vec >= 2.85e-6 & x_vec <= 4.45e-6;

% Extract the x and y data for this range
x_range1 = x_vec(idx1);
x_range2 = x_vec(idx2);
x_range3 = x_vec(idx3);
y_range1 = intensity_cross_section(idx1);
y_range2 = intensity_cross_section(idx2);
y_range3 = intensity_cross_section(idx3);

% Perform the numerical integration
integral_value1 = trapz(x_range1, y_range1);
integral_value2 = trapz(x_range2, y_range2);
integral_value3 = trapz(x_range3, y_range3);

fprintf('The integral for inner spot (from x=0um to x=1.2um) is: %e\n', integral_value1);
fprintf('The integral for first outer ring (from x=1.2um to x=2.85um) is: %e\n', integral_value2);
fprintf('The integral for second outer ring (from x=2.85um to x=4.45um) is: %e\n', integral_value3);