function [U,V,R,E] = boundary(U,V,R,E,X,Y,t,gamma,E0,osc,idx_act,v_act)
% Boundary conditions, T-S oscillator and control actuator.
% osc: oscillator struct (active, x0, length, omega, amp)
% idx_act, v_act: actuator nodes and wall-normal velocity imposed on them

P = (gamma-1) * R .* E;

% Inlet
U(1,:) = 1;
V(1,:) = 0;
E(1,:) = E0;
P(1,:) = P(2,:);

% Outlet
U(end,:) = neumann(U,2);
V(end,:) = neumann(V,2);
E(end,:) = neumann(E,2);
P(end,:) = (gamma-1)*E0;

% Outer flow
U(:,end) = neumann(U,4);
V(:,end) = neumann(V,4);
E(:,end) = neumann(E,4);
P(:,end) = neumann(P,4);

% Flat plate
U(X<0) = neumann(U(X<0,:),3);
U(X>=0,1) = 0;
V(:,1) = 0;
E(:,1) = E0;
P(:,1) = neumann(P,3);

% T-S oscillator (wall blowing and suction)
if osc.active
    idx_osc = X >= osc.x0 & X <= osc.x0 + osc.length;
    V(idx_osc,1) = osc.amp * sin(osc.omega * t);
end

% Control actuator
V(idx_act,1) = v_act;

R = P ./ ((gamma-1)*E);
end

function u = neumann(U,dir)
    NC1 = {[-1, 1], [-3/2, 2, -1/2], [-11/6, 3, -3/2, 1/3], [-25/12, 4, -3, 4/3, -1/4], [-137/60, 5, -5, 10/3, -5/4, 1/5], [-49/20, 6, -15/2, 20/3, -15/4, 6/5, -1/6]};
    order = 4;

    coeffs = -NC1{order(1)}(2:end)/NC1{order(1)}(1);

    switch dir
        case 1
            u = sum(coeffs'.*U(2:order+1,:),1);
        case 2
            u = sum(coeffs(end:-1:1)'.*U(end-order:end-1,:),1);
        case 3
            u = sum(coeffs.*U(:,2:order+1),2);
        case 4
            u = sum(coeffs(end:-1:1).*U(:,end-order:end-1),2);
    end
end
