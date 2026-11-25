function [Ut, Vt, Rt, Et] = nseq(U,V,R,E,Re,Ma,Pr,gamma,Dx,Dy)

    % Pressure and derivatives
    P = (gamma-1)*R.*E;
    
    Ux = Dx(U);
    Uy = Dy(U);
    Vx = Dx(V);
    Vy = Dy(V);
    Rx = Dx(R);
    Ry = Dy(R);
    Ex = Dx(E);
    Ey = Dy(E);
    Px = Dx(P);
    Py = Dy(P);


    % Viscous term
    suther = 110/300;

    mu1 = (1+suther)./((gamma^2 - gamma)*Ma^2*E+suther);
    mu2 = (gamma^2 - gamma)*Ma^2*E;
    mu = mu1.*mu2.*sqrt(mu2);
    
    tauxx = mu/Re .* (2*Ux-2/3*(Ux+Vy));
    tauyy = mu/Re .* (2*Vy-2/3*(Ux+Vy));
    tauxy = mu/Re .* (Uy+Vx);

    vtx = Dx(tauxx) + Dy(tauxy);
    vty = Dx(tauxy) + Dy(tauyy);

    % Heat term
    qx = -gamma/Re/Pr*mu.*Ex;
    qy = -gamma/Re/Pr*mu.*Ey;
    
    htx = Dx(qx);
    hty = Dy(qy);

    % Navier Stokes
    Rt = -R.*(Ux+Vy)-Rx.*U-Ry.*V;
    Ut = -(Ux.*U + Uy.*V) + (vtx - Px)./R;
    Vt = -(Vx.*U + Vy.*V) + (vty - Py)./R;
    Et = -(Ex.*U+Ey.*V) + (-P.*(Ux+Vy) + (tauxx.*Ux+tauyy.*Vy) + tauxy.*(Vx+Uy)-(htx+hty))./R;

end