close all; clear; clc;

data = readmatrix("refractive_index_interpolated.csv");

x = data(:,1);
y1 = data(:,2);
y2 = data(:,3);

figure;
hold on;
plot(x,y1,'LineWidth', 5);
plot(x,y2,'LineWidth', 5);
xline(550, '-.', 'LineWidth', 3);
yline(0.161,'-.', 'LineWidth', 3);
yline(4.305,'-.', 'LineWidth', 3);
plot(550,0.161, '*', 'LineWidth', 12, 'Color', 'black', 'MarkerSize', 12);
plot(550,4.305, '*', 'LineWidth', 12, 'Color', 'black', 'MarkerSize', 12);
text(535,1, 'Wavelength \lambda=550nm', 'Rotation', 90);
text(555,4.5, '(550, 4.305)');
text(555,0.36, '(550, 0.161)');
grid on;
title("Refractive Index of deposited poly-Si nano-antenna");
xlabel("Wavelength \lambda [nm]");
ylabel("Refractive Index n+ik (a.u.)");
legend("Real", "Imaginary");
fontsize(24,"points")

