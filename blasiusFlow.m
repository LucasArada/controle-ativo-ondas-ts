function [U, V] = blasiusFlow(X,Y,Re)

ybl = 0:0.0001:10;
ubl = blasius(ybl);
thetabl = 0.664155332943009;

nx = length(X);
ny = length(Y);

U = ones(nx,ny);
V = zeros(nx,ny);

% If blasiusFit defined, the blasius boundary layer will fit
% to the wall.
% The value of blasiusFit defines the maximum slope of Y0 for
% the blasius profile

for i = 1:nx
    xi = X(i);
    if xi > 0
        theta = 0.664 * sqrt(xi/Re);
        U(i,:) = interp1(ybl*theta/thetabl,ubl,Y,'spline')';
    end
end

U(U>1) = 1;
U(isnan(U)) = 1;

end

function [u] = blasius(y)

f0 = [0 0 0.33204312];

[~,u] = ode45(@blasiusEq,y,f0);

u = u(:,2);

end

function d = blasiusEq(~,f)
    d(3) = -0.5*f(1)*f(3);
    d(2) = f(3);
    d(1) = f(2);
    d = d';
end