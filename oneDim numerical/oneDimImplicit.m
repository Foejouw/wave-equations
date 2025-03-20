function w = oneDimImplicit(g, h, V_a, xMin, xMax, Dx, tMin, tMax, Dt, xss, gt)

Nx = round((xMax-xMin)/Dx);
Nt = round((tMax-tMin)/Dt);
%xs = zeros(1, 2*Nx, 1);
%xs(1:2:(2*Nx-1)) = xss(:);
%xs(2:2:2*Nx) = xss(:);



u = zeros(2*Nx, Nt);
u(1:2:2*Nx, 1) = g;
u(2:2:2*Nx, 1) = h;
for it=1:Nt
        %u(1, it) = gt(it);
    end
b = zeros(2*Nx, 1);
b(1) = -1;
b(end) = -1;
r = V_a * Dt/Dx
A = sparse(diag(2*ones(1,2*Nx)) + diag(-1*ones(1,2*Nx-1),1) + diag(-1*ones(1,2*Nx-1),-1));

B = sparse(inv(4/r^2*eye(2*Nx) + A));


for it=2:(Nt-1)
    u(:, it+1) = B*(2*(4/r^2*eye(2*Nx)  -A)*u(:, it) - 4*b) - u(:, it-1);
    g(1);
    u(2, it+1) = 0;%h(1);
    u(2*Nx-1, it+1) = 0;%g(Nx);
    u(2*Nx, it+1) = 0;%h(Nx);

end

w = u(1:2:2*Nx, :);