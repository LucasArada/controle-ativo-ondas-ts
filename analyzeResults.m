%% Attenuation and reference contamination of a closed-loop run
% Phasors at w0 by least squares (cosine, sine and constant) over the last
% N_periods periods. Run 1 is the open-loop reference.

lx = 250;
actuator_profile = 'dipole';
control_mode = 'IFFC_buffer_comp';

omega0 = 0.15;
N_periods = 10;     % 10 for runs ending at t = 2000, 5 for runs ending at t = 1500

out_dir = fullfile('results', sprintf('lx%d', lx));
h1 = load(fullfile(out_dir, 'run1_open_loop', 'history.mat'));
h3 = load(fullfile(out_dir, ['run3_' control_mode '_' actuator_profile], 'history.mat'));

% Same window for both runs
n = min(h1.i, h3.i);
t = h1.hist_t(1:n);
win = t >= t(n) - N_periods*2*pi/omega0 - 1e-9;
tw = t(win);

Z_open = lsPhasor(tw, h1.hist_p_target(win), omega0);
Z_ctrl = lsPhasor(tw, h3.hist_p_target(win), omega0);
D_f = 20*log10(abs(Z_open) / abs(Z_ctrl));
D_rms = 20*log10(std(h1.hist_p_target(win)) / std(h3.hist_p_target(win)));

% Reading used by the controller (echo removed in the compensated mode)
y_read = h3.hist_p_ref(1:n);
if strcmp(control_mode, 'IFFC_buffer_comp')
    p = load(fullfile(out_dir, ['controller_params_' actuator_profile '.mat']));
    u_del = [zeros(1, p.n_uy) h3.hist_u(1:n-p.n_uy)];
    y_read = y_read - p.Guy_mag * u_del;
end
Y_open = lsPhasor(tw, h1.hist_p_ref(win), omega0);
Y_read = lsPhasor(tw, y_read(win), omega0);
contamination = abs(Y_read - Y_open) / abs(Y_open);

fprintf('%s, %s actuator, t in [%.0f, %.0f]\n', control_mode, actuator_profile, tw(1), tw(end))
fprintf('  Target attenuation (phasor): %.1f dB\n', D_f)
fprintf('  Target attenuation (rms):    %.1f dB\n', D_rms)
fprintf('  Reference contamination:     %.3f %%\n', 100*contamination)

function X_hat = lsPhasor(t, x, omega0)
% x ~ a cos(w0 t) + b sin(w0 t) + c, phasor a - jb
A = [cos(omega0*t(:)) sin(omega0*t(:)) ones(numel(t),1)];
c = A \ x(:);
X_hat = c(1) - 1i*c(2);
end
