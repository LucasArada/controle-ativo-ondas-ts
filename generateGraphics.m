%% Figures of a closed-loop run
% Compares run 3 with run 1 (open loop). Run 2 is used for the superposition test.
% Attenuation metrics are computed by analyzeResults.m.

lx = 250;
actuator_profile = 'dipole';
control_mode = 'IFFC_buffer_comp';

gamma = 1.3;
omega0 = 0.15;
N_periods = 10;     % Window for the test tone phasor of run 2
y_wave = 0.43;      % Height of the u' maximum
t_start = 450;      % Start of the time series plots
t0 = 850;           % Control start

out_dir = fullfile('results', sprintf('lx%d', lx));
dir_run1 = fullfile(out_dir, 'run1_open_loop');
dir_run2 = fullfile(out_dir, ['run2_actuator_test_' actuator_profile]);
dir_run3 = fullfile(out_dir, ['run3_' control_mode '_' actuator_profile]);
fig_dir = fullfile('figures', sprintf('lx%d', lx));
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir)
end

T0 = 2*pi/omega0;
base = load(fullfile('data', sprintf('baseflow-%d.mat', lx)));
X = base.X(:);
Y = base.Y(:);
[~, j_wave] = min(abs(Y - y_wave));
phys = X >= 0 & X <= lx;
x_marks = [150.22 200 250];

%% Time series at the sensors

h1 = load(fullfile(dir_run1, 'history.mat'));
h3 = load(fullfile(dir_run3, 'history.mat'));
n = min(h1.i, h3.i);
t = h1.hist_t(1:n);
z_open = h1.hist_p_target(1:n);
z_ctrl = h3.hist_p_target(1:n);
y_ctrl = h3.hist_p_ref(1:n);

k = t >= t_start;
[tc, a_ref] = periodAmplitude(t(k), y_ctrl(k), T0);
[~, a_target] = periodAmplitude(t(k), z_ctrl(k), T0);
[~, a_open] = periodAmplitude(t(k), z_open(k), T0);

fig = figure('Color', 'w', 'Position', [100 100 900 620]);
subplot(2,1,1)
plot(t(k), z_open(k) - mean(z_open(k)), 'Color', [0.6 0.6 0.6])
hold on
plot(t(k), y_ctrl(k) - mean(y_ctrl(k)))
plot(t(k), z_ctrl(k) - mean(z_ctrl(k)))
xline(t0, 'k--');
hold off
xlim([t_start t(end)])
ylabel('p''')
legend('Alvo sem controle', 'Referência com controle', 'Alvo com controle', ...
    'Início do controle', 'Location', 'northoutside', 'Orientation', 'horizontal')
title('(a) Pressão nos sensores')
grid on

subplot(2,1,2)
semilogy(tc, a_open, 'Color', [0.6 0.6 0.6])
hold on
semilogy(tc, a_ref, 'o-', tc, a_target, 's-')
xline(t0, 'k--');
hold off
xlim([t_start t(end)])
xlabel('t')
ylabel('Amplitude por período')
title('(b) Amplitude por período')
grid on
print(fig, fullfile(fig_dir, ['sinais_' control_mode '_' actuator_profile '.png']), '-dpng', '-r300')

%% Wall pressure and u' along the plate

[p_open, u_open] = spatialPhasor(dir_run1, base, gamma, omega0, j_wave);
[p_ctrl, u_ctrl] = spatialPhasor(dir_run3, base, gamma, omega0, j_wave);

fig = figure('Color', 'w', 'Position', [100 100 900 620]);
subplot(2,1,1)
plot(X(phys), abs(p_open(phys)), X(phys), abs(p_ctrl(phys)))
leg = {'Sem controle', 'Com controle'};
if exist(fullfile(dir_run2, 'history.mat'), 'file')
    % Linear superposition of runs 1 and 2 with the optimal command
    h2 = load(fullfile(dir_run2, 'history.mat'));
    t2 = h2.hist_t(1:h2.i);
    w2 = t2 >= t2(end) - N_periods*T0;
    U0_hat = 2/(t2(find(w2,1,'last'))-t2(find(w2,1))) * trapz(t2(w2), h2.hist_u(w2) .* exp(-1i*omega0*t2(w2)));
    prm = load(fullfile(out_dir, ['controller_params_' actuator_profile '.mat']));
    p_act = spatialPhasor(dir_run2, base, gamma, omega0, j_wave);
    p_sup = p_open + prm.U_star/U0_hat * p_act;
    hold on
    xs = find(phys);
    xs = xs(1:3:end);
    plot(X(xs), abs(p_sup(xs)), 'o')
    hold off
    leg{end+1} = 'Superposição (1+2)';
end
for xm = x_marks
    xline(xm, 'k:');
end
xlim([0 lx])
ylabel('|p''(x, y = 0)|')
legend(leg, 'Location', 'northoutside', 'Orientation', 'horizontal')
title('(a) Pressão na parede em \omega_0')
grid on

subplot(2,1,2)
plot(X(phys), 20*log10(abs(p_open(phys)) ./ abs(p_ctrl(phys))))
yline(0, 'k:');
for xm = x_marks
    xline(xm, 'k:');
end
xlim([0 lx])
xlabel('Posição longitudinal x')
ylabel('\Delta(x) [dB]')
title('(b) Atenuação local')
grid on
print(fig, fullfile(fig_dir, ['pressao_parede_' control_mode '_' actuator_profile '.png']), '-dpng', '-r300')

fig = figure('Color', 'w', 'Position', [100 100 900 400]);
plot(X(phys), abs(u_open(phys)), X(phys), abs(u_ctrl(phys)), '--')
for xm = x_marks
    xline(xm, 'k:');
end
xlim([90 lx])
xlabel('Posição longitudinal x')
ylabel(sprintf('|u''(x, y = %.2f)|', Y(j_wave)))
legend('Sem controle', 'Com controle', 'Location', 'northoutside', 'Orientation', 'horizontal')
grid on
print(fig, fullfile(fig_dir, ['amplitude_u_' control_mode '_' actuator_profile '.png']), '-dpng', '-r300')

%% Local functions

function [p_hat, u_hat] = spatialPhasor(case_dir, base, gamma, omega0, j_wave)
% Phasors of the wall pressure and of u'(y_wave) over the last two periods
files = dir(fullfile(case_dir, 'nt_*.mat'));
names = sort({files.name});
last = load(fullfile(case_dir, names{end}), 't');
t_min = last.t - 2*2*pi/omega0;
P_base = (gamma-1) * base.R(:,1) .* base.E(:,1);
tk = [];
p = [];
u = [];
for k = numel(names):-1:1
    d = load(fullfile(case_dir, names{k}), 'U', 'R', 'E', 't');
    if d.t < t_min
        break
    end
    tk(end+1) = d.t;
    p(:,end+1) = (gamma-1) * d.R(:,1) .* d.E(:,1) - P_base;
    u(:,end+1) = d.U(:,j_wave) - base.U(:,j_wave);
end
[tk, s] = sort(tk);
w = exp(-1i*omega0*tk);
p_hat = 2/(tk(end)-tk(1)) * trapz(tk, p(:,s) .* w, 2);
u_hat = 2/(tk(end)-tk(1)) * trapz(tk, u(:,s) .* w, 2);
end

function [tc, a] = periodAmplitude(t, x, T0)
% Amplitude of each period, sqrt(2) times the rms
np = floor((t(end)-t(1))/T0);
tc = zeros(1,np);
a = zeros(1,np);
for k = 1:np
    m = t >= t(1)+(k-1)*T0 & t < t(1)+k*T0;
    tc(k) = t(1) + (k-0.5)*T0;
    a(k) = sqrt(2) * std(x(m), 1);
end
end
