function u = controller(t, i, p_ref, u_hist, ctrl)
% Actuator command at time step i.
% p_ref: reference sensor history, u_hist: command history, ctrl: control parameters

switch ctrl.mode
    case 'off'
        u = 0;

    case 'test_tone'
        % Run 2: pure tone used to identify the actuator paths
        u = ctrl.U0 * sin(ctrl.omega0 * t);

    case 'IFFC_phasor'
        % Optimal command synchronized with the oscillator: Re{U* e^(j w0 t)}
        u = envelope(t, ctrl) * real(ctrl.U_star * exp(1i * ctrl.omega0 * t));

    case 'IFFC_buffer'
        % Real time: u(t) = s(t) |K| y(t - dt_n)
        k = i - ctrl.n_delay;
        if k >= 1
            u = envelope(t, ctrl) * ctrl.K_mag * p_ref(k);
        else
            u = 0;
        end

    case 'IFFC_buffer_comp'
        % Real time with echo compensation:
        % y_c(t) = y(t) - |G_uy| u(t - dt_uy),  u(t) = s(t) |K| y_c(t - dt_n)
        k = i - ctrl.n_delay;
        if k >= 1
            k_u = k - ctrl.n_uy;
            if k_u >= 1
                echo = ctrl.Guy_mag * u_hist(k_u);
            else
                echo = 0;
            end
            y_c = p_ref(k) - echo;
            u = envelope(t, ctrl) * ctrl.K_mag * y_c;
        else
            u = 0;
        end

    otherwise
        error('Unknown control mode: %s', ctrl.mode)
end
end

function s = envelope(t, ctrl)
% Smooth start between t0 and t0 + T_ramp
if t < ctrl.t0
    s = 0;
elseif t <= ctrl.t0 + ctrl.T_ramp
    s = 0.5 * (1 - cos(pi * (t - ctrl.t0) / ctrl.T_ramp));
else
    s = 1;
end
end
