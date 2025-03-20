clear all
%%

xMin = 0;
xMax = 10;
tMin = 0;
tMax = 3;
Dx = 0.1;
Dt = 0.0005;
Nx = round((xMax-xMin)/Dx);
Nt = round((tMax-tMin)/Dt);
xs = linspace(xMin, xMax, round((xMax-xMin)/Dx));
ts = linspace(tMin, tMax, round((tMax-tMin)/Dt));

%%
% D'Alembert solution
% Solving u_tt = V_a^2 * v_zz
% with V_a = B/sqrt(mu*rho)

% BC: u(x, 0) = g(x)
g = @(x) exp(-(x-5).^2);
%g = @(x) sin(x)/2;
%g = @(x) 3.*exp(-(x-6).^2);
%g = @(x) exp(-144.*(x-6).^2);
%g = @(x) 0;

% BC: u_t(x, 0) = h(x)
%h = @(x) exp(-(x-5).^2);
h = @(x) 0;
%h = @(x) -2.*(x).*exp(-(x-5).^2);

V_a = 100;

%uALem = @(x, t) 1/2 .* (g(x-V_a.*t) + g(x+V_a.*t)) ;
%+ 1/(2.*V_a).*integral(h, x-V_a.*t, x+V_a.*t);


uALem = @(x, t) 1/2 .* (piece(xMin, xMax, x-V_a.*t, g(x-V_a.*t)) + piece(xMin, xMax, x+V_a.*t, g(x+V_a.*t)));

%%
% Numerical analysis


uNum = oneDimImplicit(g(xs), h(xs), V_a, xMin, xMax, Dx, tMin, tMax, Dt, xs, g(ts));


%%
% Plot Numerical


s = surf(ts, xs, uNum);
set(s, 'LineStyle', 'none');
%zlim([-2, 2])
xlabel("t")
ylabel("x")

%%

% Plot animation

f = figure;
ax = gca;
xlabel('z')
ylabel('v')
ylim([-50 50]);
xlim([2.5, 7.5]);
ax.NextPlot = 'replaceChildren';
M(Nt) = struct('cdata',[],'colormap',[]);
f.Visible = 'on';
v = VideoWriter("Test13", "MPEG-4");
open(v);
for it=1:1:Nt
    plot(xs, uNum(:, it));
    drawnow limitrate;
    M(it) = getframe(f);
    %writeVideo(v, M(it));
end
%open(v)
%writeVideo(v, M(:));
close(v)

%%


it = round(Nt/3);
it = 1;
alam = zeros(Nx, 1);
for ix=1:Nx
    alam(ix) = uALem(xs(ix), ts(it));
end

plot(xs, alam, "red")
hold on
plot(xs, uNum(:, it), "blue")
ylim([-2, 2]);


%%

am = zeros(Nx, Nt);
for it=1:1:Nt/2
    for ix=1:Nx
        am(ix, it) = uALem(xs(ix), ts(it));
    end
end

%%
% Plot animation with D'Alembert

f = figure;
ax = gca;
xlabel('z')
ylabel('v')
ylim([-1 1]);
ax.NextPlot = 'replaceChildren';
M(Nt) = struct('cdata',[],'colormap',[]);
f.Visible = 'on';
v = VideoWriter("Test13", "MPEG-4");
open(v);

for it=1:1:Nt
    %plot(xs, uNum(:, it));
    %hold on
    plot(xs, am(:, it));
    %hold off
    drawnow limitrate;
    M(it) = getframe(f);
    %writeVideo(v, M(it));
end
%open(v)
%writeVideo(v, M(:));
close(v)


%%%
%%


function y = piece(min, max, x, y)
    period = 2*(max-min);
    if mod(x, period) < min/2
        y = y;
    end
end

