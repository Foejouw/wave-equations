clear all



%%
%Nx = 1000;
%Nt = 100;
xMin = 0;
xMax = 10;
tMin = 0;
tMax = 1.5;
%Dx = (xMax-xMin)/Nx;
%Dt = (tMax-tMin)/Nt;
Dt = 0.005;
Dx = 0.1;
Nx = round((xMax-xMin)/Dx);
Nt = round((tMax-tMin)/Dt);
xs = linspace(xMin, xMax, Nx);
ts = linspace(tMin, tMax, Nt);
uPrev = zeros( Nx,1);
uNow = zeros( Nx, 1);
uNext = zeros( Nx,1);
A = -2*eye(Nx, Nx) + diag(ones(1, Nx-1), -1) + diag(ones(1, Nx-1), 1);

V = 2;


%u0 = @(x) exp(-(x-(xMax-xMin)./2).^2);
u0 = @(x) sin(pi.*x) ; 
g = @(x) 0;
uALem = @(x, t) 1/2 .* (u0(mod(x-V.*t, xMax)) + u0(mod(x+V.*t, xMax)));


%% WAVE Numerical
uPrev(:) = u0(xs);
uNow(:) = u0(xs);
uPrev(1) = g(xs(1));
uPrev(Nx) = g(xs(Nx));
uNow(1) = g(xs(1));
uNow(Nx) = g(xs(Nx));
u(1, :) = uNow;


r = Dt*V/Dx
f1 = figure;
figure(f1);
handle_line = plot(xs,uNow(:),'LineWidth',1);
axis([0,xMax-xMin,-2,2]);
xlabel('x'); ylabel('u');
title('Wave equation');





for it=2:Nt
    uNext(:) = r^2 * A * uNow + 2.*uNow - uPrev;
    uPrev = uNow;
    uNow = uNext;
    uNow(1) = g(xs(1));
    uNow(Nx) = g(xs(Nx));
    uPrev(1) = g(xs(1));
    uPrev(Nx) = g(xs(Nx));
    u(it, :) = uNow;
    handle_line.YData = uNow(:);

    drawnow;
end
%%
figure(1);
surf(xs, ts, u, 'LineStyle','none')

test = @(x,t) sin(pi.*x).*cos(2.*pi.*t);
figure(2);
surf(xs, ts, test(xs, ts'))
%%

u - test(xs, ts');
surf(xs, ts, u - test(xs, ts'), 'LineStyle','none')
%%

%%


test = @(x,t) sin(pi.*x).*cos(2.*pi.*t);

semilogy(ts, max(abs(u-test(xs, ts')), [], 2), Color='red', LineWidth=1.2)
hold on
yline(0.04, 'Color', [0.9290 0.6940 0.1250], 'LineWidth',1.2)
ylim([10^(-4) 1])
xlabel('t')
ylabel('error')
legend({'$\max_x \left| v_{num} - \tilde{v} \right|$', '$y=0.04$'}, 'Interpreter','latex', 'FontSize',11)

%%


%% ERROR in time, animation

r = Dt*V/Dx

uPrev(:) = u0(xs);
uNow(:) = u0(xs);
uPrev(1) = g(xs(1));
uPrev(Nx) = g(xs(Nx));
uNow(1) = g(xs(1));
uNow(Nx) = g(xs(Nx));


handle_line2 = plot(xs, uNow(:) - transpose(uALem(xs, ts(1))),'LineWidth',1);
%axis([0,xMax-xMin,-0.01,0.01]);
xlabel('x'); ylabel('error');
title('Wave equation');
u = zeros(1, Nx);

for it=2:Nt-1
    uNext(:) = r^2 * A * uNow + 2.*uNow - uPrev;
    uPrev = uNow;
    uNow = uNext;
    uNow(1) = g(xs(1));
    uNow(Nx) = g(xs(Nx));
    uPrev(1) = g(xs(1));
    uPrev(Nx) = g(xs(Nx));
    handle_line2.YData = uNow(:) - transpose(uALem(xs, ts(it)));
    drawnow;
end

%% ERROR in time, logplot

r = Dt*V/Dx

uPrev(:) = u0(xs);
uNow(:) = u0(xs);
uPrev(1) = g(xs(1));
uPrev(Nx) = g(xs(Nx));
uNow(1) = g(xs(1));
uNow(Nx) = g(xs(Nx));
    

maxErr = zeros(1, Nt);

for it=2:Nt
    uNext(:) = r^2 * A * uNow + 2.*uNow - uPrev;
    uPrev = uNow;
    uNow = uNext;
    uNow(1) = g(xs(1));
    uNow(Nx) = g(xs(Nx));
    uPrev(1) = g(xs(1));
    uPrev(Nx) = g(xs(Nx));
    maxErr(it) = max(abs(uNow(:) - transpose(uALem(xs, ts(it)))));
end

%%

semilogy(ts, maxErr)
xlabel('t')
ylabel('error')
title(strcat('Error numerical vs theoretical (r = ', num2str(r), ')'))

