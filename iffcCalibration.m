%% IFFC calibration from runs 1 and 2
% Extracts the phasors at w0, identifies the transfer functions and saves the
% controller parameters used by run 3.

lx = 250;                       % Same as in runSolver.m
actuator_profile = 'dipole';    % Same as in runSolver.m

dt = 1e-2;
omega0 = 0.15;
N_periods = 10;     % Phasor window (last N periods)
R_penal = 0;        % Control effort penalty

%% Load sensor histories

out_dir = fullfile('results', sprintf('lx%d', lx));
dir_run1 = fullfile(out_dir, 'run1_open_loop');
dir_run2 = fullfile(out_dir, ['run2_actuator_test_' actuator_profile]);

h1 = loadHistory(dir_run1);
h2 = loadHistory(dir_run2);

% Sensor positions on the mesh
mesh = load(fullfile(dir_run1, sprintf('nt_%010d.mat', 0)), 'X');
[~, i_ref] = min(abs(mesh.X - 150));
[~, i_target] = min(abs(mesh.X - 250));
Delta_x = mesh.X(i_target) - mesh.X(i_ref);

%% Phasors and transfer functions

Y_hat     = phasor(h1.t, h1.p_ref,    omega0, N_periods);   % Reference, run 1
Z_un_hat  = phasor(h1.t, h1.p_target, omega0, N_periods);   % Target, run 1
U0_hat    = phasor(h2.t, h2.u,        omega0, N_periods);   % Test tone, run 2
Y_act_hat = phasor(h2.t, h2.p_ref,    omega0, N_periods);   % Reference, run 2
Z_act_hat = phasor(h2.t, h2.p_target, omega0, N_periods);   % Target, run 2

G_yz = Z_un_hat / Y_hat;        % T-S path, reference -> target
G_uz = Z_act_hat / U0_hat;      % Control path, actuator -> target
G_uy = Y_act_hat / U0_hat;      % Parasitic path, actuator -> reference

%% Controller synthesis

if R_penal == 0
    K_hat = -G_yz / G_uz;
else
    K_hat = -conj(G_uz) * G_yz / (abs(G_uz)^2 + R_penal);
end
K_mag = abs(K_hat);
phi_K = angle(K_hat);

% Phasor command (cosine convention)
U_star = -Z_un_hat / G_uz;
A_ctrl = abs(U_star);
phi_ctrl = angle(U_star);

% Causal delay equivalent to the phase of K at w0
dt_n = (2*pi*ceil(phi_K/(2*pi)) - phi_K) / omega0;
n_delay = round(dt_n / dt);

% Upstream leakage and parasitic loop
eps_leak = abs(Y_act_hat) / abs(Y_hat);
L = K_hat * G_uy;

% Echo compensation: gain |G_uy| and causal delay with the phase of G_uy
Guy_mag = abs(G_uy);
dt_uy = mod(-angle(G_uy), 2*pi) / omega0;
n_uy = round(dt_uy / dt);

%% Report

deg = 180/pi;
fprintf('\nIFFC calibration, lx = %d, %s actuator\n', lx, actuator_profile)
fprintf('  Y_hat     = %.4e  %+8.2f deg\n', abs(Y_hat), deg*angle(Y_hat))
fprintf('  Z_un_hat  = %.4e  %+8.2f deg\n', abs(Z_un_hat), deg*angle(Z_un_hat))
fprintf('  G_yz      = %.4f      %+8.2f deg\n', abs(G_yz), deg*angle(G_yz))
fprintf('  G_uz      = %.4f      %+8.2f deg\n', abs(G_uz), deg*angle(G_uz))
fprintf('  G_uy      = %.4f      %+8.2f deg\n', abs(G_uy), deg*angle(G_uy))
fprintf('  K         = %.4f      %+8.2f deg\n', K_mag, deg*phi_K)
fprintf('  U*        = %.4e  %+8.2f deg\n', A_ctrl, deg*phi_ctrl)
fprintf('  dt_n      = %.2f (%d steps)\n', dt_n, n_delay)
fprintf('  eps_leak  = %.2f %%\n', 100*eps_leak)
fprintf('  |L|       = %.3f\n', abs(L))
fprintf('  dt_uy     = %.2f (%d steps)\n', dt_uy, n_uy)

% Phase unwrapping: T-S wave has phase speed between 0.30 and 0.40
fprintf('\n  n   alpha    c/U_inf\n')
for n = 0:8
    alpha_n = (abs(angle(G_yz)) + 2*pi*n) / Delta_x;
    c_n = omega0 / alpha_n;
    mark = '';
    if c_n >= 0.30 && c_n <= 0.40
        mark = '  <- T-S';
    end
    fprintf('  %d   %.4f   %.4f%s\n', n, alpha_n, c_n, mark)
end

%% Save

params_file = fullfile(out_dir, ['controller_params_' actuator_profile '.mat']);
save(params_file, 'omega0', 'dt', 'N_periods', 'R_penal', ...
    'Y_hat', 'Z_un_hat', 'Y_act_hat', 'Z_act_hat', 'G_yz', 'G_uz', 'G_uy', ...
    'K_hat', 'K_mag', 'phi_K', 'U_star', 'A_ctrl', 'phi_ctrl', ...
    'dt_n', 'n_delay', 'eps_leak', 'L', 'Guy_mag', 'dt_uy', 'n_uy')
fprintf('\nSaved %s\n', params_file)

%% Local functions

function h = loadHistory(case_dir)
% Sensor and command histories up to the last completed step
d = load(fullfile(case_dir, 'history.mat'));
h.t = d.hist_t(1:d.i);
h.p_ref = d.hist_p_ref(1:d.i);
h.p_target = d.hist_p_target(1:d.i);
h.u = d.hist_u(1:d.i);
end

function X_hat = phasor(t, x, omega0, N_periods)
% Coherent correlation over the last N_periods periods
T_win = N_periods * 2*pi/omega0;
if t(end) - t(1) < T_win
    error('Run is shorter than %d periods', N_periods)
end
k = t >= t(end) - T_win - 0.5*(t(2)-t(1));
tw = t(k);
X_hat = 2/(tw(end)-tw(1)) * trapz(tw, x(k) .* exp(-1i*omega0*tw));
end
