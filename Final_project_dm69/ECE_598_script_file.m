%% Ece 598 final project script

clc
clear all
close all

%%  Inputs for 3-level inverter (all types)

fsw = 1000;                                                               % Switching frequency (Hz)
fe = 100;                                                                   % sine frequency (Hz)
Vm = 200;                                                                  % input voltage peak (V)
Vdc = 400;                                                                 % Input DC voltage (V)



%% Motor parameters


%% Setting up the simulink model



% % Input parameters
% Vpk = 220*sqrt(2/3);
% fe = 60;
p = 4;


% Circuit parameters
Lsl = 9.6e-3;
Lrl_dash = 9.6e-3;
M = 165e-3;
turns_ratio = 0.682;
Rr = 4.23;
Rs = 3.57;

% Mechanical parameters
J = 0.01;
tau_L = 2; %1.1e-2;
B = 0*0.0025;


% Derived quantities

Lm = M;
Ls_dash = M + Lsl;
Lr_dash = M + Lrl_dash;

sig = 1 - Lm^2/(Ls_dash*Lr_dash);


% Inductance matrix

L = [Ls_dash 0 Lm 0; 0 Ls_dash 0 Lm; Lm 0 Lr_dash 0; 0 Lm 0 Lr_dash];

L_inv = inv(L);


%% Anlysing spectrum of v_az for DCML inverter


% Extract time and voltage from simulation output
time = out.tout;
voltage = out.v_az;

% Compute sampling time and frequency
Ts = time(2) - time(1);    % Sampling interval
Fs = 1 / Ts;               % Sampling frequency

% Optional: Trim the steady-state portion (e.g., after 0.2 seconds)
steady_idx = find(time >= 0.01);
v_ss = voltage(steady_idx);
t_ss = time(steady_idx);

% ========== THD Calculation using Signal Processing Toolbox ==========
% Requires: thd() function
thd_dB = thd(v_ss, Fs);               % THD in dB
thd_percent = 100 * 10^(thd_dB/20);   % Convert dB to %
fprintf('THD = %.2f %%\n', thd_percent);

% ========== FFT Spectrum Plot ==========
L = length(v_ss);
Y = fft(v_ss);
P2 = abs(Y / L);                      % Two-sided spectrum
P1 = P2(1:floor(L/2)+1);              % One-sided spectrum
P1(2:end-1) = 2*P1(2:end-1);
f = Fs * (0:(L/2)) / L;

% % Plot the frequency spectrum
% figure;
% plot(f, P1);
% xlabel('Frequency (Hz)');
% ylabel('Magnitude');
% title('Frequency Spectrum of V_{AZ}');
% grid on;
% xlim([0, 2000]);  % adjust as needed based on fundamental frequency


figure;
stem(f/fe, P1./Vdc, 'k', ...
    'Marker', 'none', ...
    'LineWidth', 1);     % Thin lines to mimic bars

xlabel('$n$', 'Interpreter', 'latex');
ylabel('$V_{AZ}/V_{dc}$', 'Interpreter', 'latex');
xlim([-2 65]);
ylim([0 0.6]);
grid on;
title('Frequency spectrum of phase voltage of CHB')

x_range = xlim;
y_range = ylim;

x_pos = x_range(2) - 0.05 * diff(x_range);  % 5% from right edge
y_pos = y_range(2) - 0.05 * diff(y_range);  % 5% from top edge

text(x_pos, y_pos, sprintf('THD = %.2f%%', thd_percent), ...
    'FontSize', 10, ...
    'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'top');

set(gca, 'FontName', 'Times New Roman', 'FontSize', 10);


%% Anlysing spectrum of v_ab for DCML inverter


% Extract time and voltage from simulation output
time = out.tout;
voltage = out.v_ab;

% Compute sampling time and frequency
Ts = time(2) - time(1);    % Sampling interval
Fs = 1 / Ts;               % Sampling frequency

% Optional: Trim the steady-state portion (e.g., after 0.2 seconds)
steady_idx = find(time >= 0.01);
v_ss = voltage(steady_idx);
t_ss = time(steady_idx);

% ========== THD Calculation using Signal Processing Toolbox ==========
% Requires: thd() function
thd_dB = thd(v_ss, Fs);               % THD in dB
thd_percent = 100 * 10^(thd_dB/20);   % Convert dB to %
fprintf('THD = %.2f %%\n', thd_percent);

% ========== FFT Spectrum Plot ==========
L = length(v_ss);
Y = fft(v_ss);
P2 = abs(Y / L);                      % Two-sided spectrum
P1 = P2(1:floor(L/2)+1);              % One-sided spectrum
P1(2:end-1) = 2*P1(2:end-1);
f = Fs * (0:(L/2)) / L;

% % Plot the frequency spectrum
% figure;
% plot(f, P1);
% xlabel('Frequency (Hz)');
% ylabel('Magnitude');
% title('Frequency Spectrum of V_{AZ}');
% grid on;
% xlim([0, 2000]);  % adjust as needed based on fundamental frequency



figure;
stem(f/fe, P1./Vdc, 'k', ...
    'Marker', 'none', ...
    'LineWidth', 1);     % Thin lines to mimic bars

xlabel('$n$', 'Interpreter', 'latex');
ylabel('$V_{AB}/V_{dc}$', 'Interpreter', 'latex');
xlim([-2 65]);
ylim([0 0.6]);
grid on;
title('Frequency spectrum of line-to-line voltage of CHB')

x_range = xlim;
y_range = ylim;

x_pos = x_range(2) - 0.05 * diff(x_range);  % 5% from right edge
y_pos = y_range(2) - 0.05 * diff(y_range);  % 5% from top edge

text(x_pos, y_pos, sprintf('THD = %.2f%%', thd_percent), ...
    'FontSize', 10, ...
    'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'top');

set(gca, 'FontName', 'Times New Roman', 'FontSize', 10);
