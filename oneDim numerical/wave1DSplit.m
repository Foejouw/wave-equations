
%%
% Initialization

%Nx = 1000;
%Nt = 100;
xMin = 0;
xMax = 10;
tMin = 0;
tMax = 40;
%Dx = (xMax-xMin)/Nx;
%Dt = (tMax-tMin)/Nt;
Dt = 0.001;
Dx = 0.1;
Nx = round((xMax-xMin)/Dx);
Nt = round((tMax-tMin)/Dt);
xs = linspace(xMin, xMax, Nx);
ts = linspace(tMin, tMax, Nt);

V = 1;


vNow = zeros(1, Nx);
vNext = zeros(1, Nx);
UNow = zeros(1, Nx);
UNext = zeros(1, Nx);


%%
% Boundary conditions

v_x_t0 = @(x) exp(-(x - (xMax-xMin)/2).^2); % v(x, 0)
v_x0_t = @(t) 0; % v(xMin, t)
v_x1_t = @(t) 0; % v(xMax, t)
%dv_x_t0 = @(x) cos(sin(x)); % v_t(x, 0)
dv_x_t0 = @(x) 0; % v_t(x, 0)
dv_x0_t = @(t) 0; % v_t(xMin, t)
dv_x1_t = @(t) 0; % v_t(xMax, t)

vNow(:) = v_x_t0(xs);
vNow(1) = v_x0_t(ts(1));
vNow(Nx) = v_x1_t(ts(1));
UNow(:) = dv_x_t0(xs);
UNow(1) = dv_x0_t(ts(1));
UNow(Nx) = dv_x1_t(ts(1));



%%

(V./Dx).^2 .* Dt
handle_line = plot(xs,vNow(:),'LineWidth',1);
axis([0,xMax-xMin,-4,4]);
xlabel('x'); ylabel('v');
title('Wave equation');

for it=2:Nt

    vNext = vNow + Dt.*UNow;
    UNext(2:Nx-1) = UNow(2:Nx-1) + (V./Dx).^2 .* Dt .* (vNow(3:Nx) - 2.*vNow(2:Nx-1) + vNow(1:Nx-2));

    vNow = vNext;
    UNow = UNext;

    % BC:
    vNow(1) = v_x0_t(ts(it));
    vNow(Nx) = v_x1_t(ts(it));
    UNow(1) = dv_x0_t(ts(it));
    UNow(Nx) = dv_x1_t(ts(it));

    handle_line.YData = vNow(:);
    drawnow;

end