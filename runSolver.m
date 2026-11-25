%% Case selection
% Run order: 1 (open loop) -> 2 (actuator test) -> iffcCalibration.m -> 3 (closed loop)
% Run 0 regenerates the base flow from the Blasius profile (optional).

run_number = 3;                     % 0, 1, 2 or 3
lx = 250;                           % Physical domain length: 250 or 500
actuator_profile = 'dipole';        % 'uniform' or 'dipole'
control_mode = 'IFFC_buffer_comp';  % Run 3: 'IFFC_phasor', 'IFFC_buffer' or 'IFFC_buffer_comp'

%% Define parameters

Re = 1000;
Ma = 0.5;
Pr = 0.72;
gamma = 1.3;

dt = 1e-2;
tf = 2000;
nsave = 100;
nplot = 0;      % Live plot interval in steps (0 disables)

ly = 10;

% T-S oscillator
osc.x0 = 90;
osc.length = 20;
osc.omega = 0.15;
osc.amp = 1e-4;

% Wall pressure sensors
x_ref = 150;
x_target = 250;

% Actuator: one weight per node (width 10 covers 5 nodes)
x_act = 200;
act_width = 10;
switch actuator_profile
    case 'uniform'
        w_act = [1; 1; 1; 1; 1];
    case 'dipole'
        w_act = [1; 1; 0; -1; -1];
end

% Control
ctrl.omega0 = osc.omega;
ctrl.U0 = 1e-4;                         % Run 2 test tone amplitude
ctrl.t0 = 850;                          % Control start
ctrl.T_ramp = 4 * 2*pi/ctrl.omega0;     % Smooth start duration (4 periods)

% Mesh
mesh.x.d0 = 2;
mesh.y.d0 = 1;
mesh.x.type = 'attractors';
mesh.x.attractorPoints = [];
mesh.x.attractorStrength = [];
mesh.x.attractorSize = [];
mesh.x.attractorRegions = [-30 -10 10 30 1];
mesh.x.matchFixed = 2;
mesh.x.periodic = false;
mesh.x.fixPeriodicDomainSize = false;
mesh.x.extraRefinement = 0;

mesh.y.type = 'attractors';
mesh.y.attractorPoints = [];
mesh.y.attractorStrength = [];
mesh.y.attractorSize = [];
mesh.y.attractorRegions = [-10 -0.5 0.5 10 12];
mesh.y.matchFixed = 2;
mesh.y.periodic = false;
mesh.y.fixPeriodicDomainSize = false;
mesh.y.extraRefinement = 0;

mesh.x.buffer.i.n = 40;
mesh.x.buffer.i.type = 'sigmoid';
mesh.x.buffer.i.stretching = 30;
mesh.x.buffer.i.transition = 0.2;

mesh.x.buffer.f.n = 30;
mesh.x.buffer.f.type = 'sigmoid';
mesh.x.buffer.f.stretching = 15;
mesh.x.buffer.f.transition = 0.2;

mesh.y.buffer.i.n = 0;

mesh.y.buffer.f.n = 40;
mesh.y.buffer.f.type = 'sigmoid';
mesh.y.buffer.f.stretching = 30;
mesh.y.buffer.f.transition = 0.2;

% Sponge: relaxation towards the base flow inside the buffer zones
sponge_strength = 1;

%% Case setup

out_dir = fullfile('results', sprintf('lx%d', lx));
base_file = fullfile('data', sprintf('baseflow-%d.mat', lx));
params_file = fullfile(out_dir, ['controller_params_' actuator_profile '.mat']);

switch run_number
    case 0
        case_name = 'run0_baseflow';
        osc.active = false;
        ctrl.mode = 'off';
    case 1
        case_name = 'run1_open_loop';
        osc.active = true;
        ctrl.mode = 'off';
    case 2
        case_name = ['run2_actuator_test_' actuator_profile];
        osc.active = false;
        ctrl.mode = 'test_tone';
    case 3
        case_name = ['run3_' control_mode '_' actuator_profile];
        osc.active = true;
        ctrl.mode = control_mode;
        if ~exist(params_file, 'file')
            % Parameters identified in the monograph, if runs 1 and 2 were not done
            params_file = fullfile('data', sprintf('controller-params-%d-%s.mat', lx, actuator_profile));
            if ~exist(params_file, 'file')
                error('No controller parameters found. Run iffcCalibration.m first.')
            end
            fprintf('Using %s\n', params_file)
        end
        p = load(params_file);
        if abs(p.dt - dt) > eps
            error('Controller parameters were identified with dt = %g', p.dt)
        end
        ctrl.U_star = p.U_star;
        ctrl.K_mag = p.K_mag;
        ctrl.n_delay = p.n_delay;
        ctrl.Guy_mag = p.Guy_mag;
        ctrl.n_uy = p.n_uy;
end

case_dir = fullfile(out_dir, case_name);
if ~exist(case_dir, 'dir')
    mkdir(case_dir)
end

%% Generate mesh

[mesh.X, mesh.x, mesh.nx] = generateMesh(0,lx,mesh.x,'X');
[mesh.Y, mesh.y, mesh.ny] = generateMesh(0,ly,mesh.y,'Y');

X = mesh.X;
Y = mesh.Y;

T = 0:dt:tf;

nt = length(T);
nx = length(X);
ny = length(Y);

%% Initialize flow

E0 = 1/((gamma^2-gamma)*Ma^2);
if run_number == 0
    [U, V] = blasiusFlow(X,Y,Re);
    R = ones(nx,ny);
    E = E0*ones(nx,ny);
else
    load(base_file,'U','V','R','E')
end

% Reference state for the sponge and for the pressure perturbation
base.U = U;
base.V = V;
base.R = R;
base.E = E;

t = 0;

file_name = fullfile(case_dir, sprintf('nt_%010d.mat', 0));
save(file_name,'U','V','R','E','X','Y','t')

%% Initialize the derivative and filter functions
[Dx, Dy, F] = initDerivatives(X,Y);

%% Sensors and actuator

[~, i_ref] = min(abs(X - x_ref));
[~, i_target] = min(abs(X - x_target));
P_base = (gamma-1) * base.R(:,1) .* base.E(:,1);

idx_act = find(X >= x_act - act_width/2 & X <= x_act + act_width/2);
if numel(idx_act) ~= numel(w_act)
    error('Actuator has %d nodes but w_act has %d weights', numel(idx_act), numel(w_act))
end

hist_t = zeros(1, nt);
hist_p_ref = zeros(1, nt);
hist_p_target = zeros(1, nt);
hist_u = zeros(1, nt);

%% Sponge mask

sx = ones(nx,1);
sy = ones(1,ny);
sx(1:mesh.x.buffer.i.n) = 0.5+0.5*sin(linspace(-pi/2,pi/2,mesh.x.buffer.i.n));
sx(nx-mesh.x.buffer.f.n+1:end) = 0.5+0.5*sin(linspace(pi/2,-pi/2,mesh.x.buffer.f.n));
sy(ny-mesh.y.buffer.f.n+1:end) = 0.5+0.5*sin(linspace(pi/2,-pi/2,mesh.y.buffer.f.n));
sponge = sponge_strength * (1 - sx.*sy);

%% Apply boundary conditions
[U,V,R,E] = boundary(U,V,R,E,X,Y,T(1),gamma,E0,osc,idx_act,0);

%% Time stepping
for i = 1:nt

    t = T(i);

    % Command is held constant over the four RK4 stages
    u = controller(t, i, hist_p_ref, hist_u, ctrl);
    hist_u(i) = u;
    v_act = w_act * u;

    [K1U, K1V, K1R, K1E] = nseq(U,V,R,E,Re,Ma,Pr,gamma,Dx,Dy);
    Unew = U + dt/2*K1U;
    Vnew = V + dt/2*K1V;
    Rnew = R + dt/2*K1R;
    Enew = E + dt/2*K1E;
    [Unew,Vnew,Rnew,Enew] = boundary(Unew,Vnew,Rnew,Enew,X,Y,t+dt/2,gamma,E0,osc,idx_act,v_act);

    [K2U, K2V, K2R, K2E] = nseq(Unew,Vnew,Rnew,Enew,Re,Ma,Pr,gamma,Dx,Dy);
    Unew = U + dt/2*K2U;
    Vnew = V + dt/2*K2V;
    Rnew = R + dt/2*K2R;
    Enew = E + dt/2*K2E;
    [Unew,Vnew,Rnew,Enew] = boundary(Unew,Vnew,Rnew,Enew,X,Y,t+dt/2,gamma,E0,osc,idx_act,v_act);

    [K3U, K3V, K3R, K3E] = nseq(Unew,Vnew,Rnew,Enew,Re,Ma,Pr,gamma,Dx,Dy);
    Unew = U + dt*K3U;
    Vnew = V + dt*K3V;
    Rnew = R + dt*K3R;
    Enew = E + dt*K3E;
    [Unew,Vnew,Rnew,Enew] = boundary(Unew,Vnew,Rnew,Enew,X,Y,t+dt,gamma,E0,osc,idx_act,v_act);

    [K4U, K4V, K4R, K4E] = nseq(Unew,Vnew,Rnew,Enew,Re,Ma,Pr,gamma,Dx,Dy);

    Unew = U + dt/6 * (K1U + 2*K2U + 2*K3U + K4U);
    Vnew = V + dt/6 * (K1V + 2*K2V + 2*K3V + K4V);
    Rnew = R + dt/6 * (K1R + 2*K2R + 2*K3R + K4R);
    Enew = E + dt/6 * (K1E + 2*K2E + 2*K3E + K4E);

    [Unew,Vnew,Rnew,Enew] = boundary(Unew,Vnew,Rnew,Enew,X,Y,t+dt,gamma,E0,osc,idx_act,v_act);

    Unew = F(Unew);
    Vnew = F(Vnew);
    Rnew = F(Rnew);
    Enew = F(Enew);

    Unew = Unew - dt*sponge.*(Unew-base.U);
    Vnew = Vnew - dt*sponge.*(Vnew-base.V);
    Rnew = Rnew - dt*sponge.*(Rnew-base.R);
    Enew = Enew - dt*sponge.*(Enew-base.E);

    [Unew,Vnew,Rnew,Enew] = boundary(Unew,Vnew,Rnew,Enew,X,Y,t+dt,gamma,E0,osc,idx_act,v_act);

    U = Unew;
    V = Vnew;
    R = Rnew;
    E = Enew;

    % Wall pressure perturbation at the sensors
    P_wall = (gamma-1) * R(:,1) .* E(:,1) - P_base;
    hist_t(i) = t;
    hist_p_ref(i) = P_wall(i_ref);
    hist_p_target(i) = P_wall(i_target);

    if nplot > 0 && mod(i,nplot) == 0
        figure(1)
        subplot(2,1,1)
        pcolor(X,Y,(U-base.U)')
        axis([0 lx 0 ly])
        shading interp
        colorbar
        title(sprintf('u'' em t = %.1f', t))
        subplot(2,1,2)
        plot(hist_t(1:i), hist_p_ref(1:i), hist_t(1:i), hist_p_target(1:i))
        legend('Referência', 'Alvo')
        xlabel('t')
        ylabel('p'' na parede')
        drawnow
    end

    if mod(i,nsave) == 0
        file_name = fullfile(case_dir, sprintf('nt_%010d.mat', i));
        save(file_name,'U','V','R','E','X','Y','t')
        fprintf('%s  t = %.2f\n', case_name, t)
    end

    failed = ~isreal(U) || any(~isfinite(U(:)));
    if failed || mod(i,100*nsave) == 0 || i == nt
        save(fullfile(case_dir, 'history.mat'), 'hist_t', 'hist_p_ref', 'hist_p_target', 'hist_u', 'i')
    end
    if failed
        warning('Simulation diverged at t = %.2f', t)
        break
    end

end

if run_number == 0
    save(fullfile(case_dir, 'baseflow.mat'), 'U','V','R','E','X','Y')
    fprintf('Base flow saved in %s. Copy it to %s to use it.\n', case_dir, base_file)
end
