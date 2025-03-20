% wave.m

clear all;
close all;
x0 = 0;
xf = 1;
umin = -1.5;
umax = 1.5;
t0 = 0;
tf = 10;
f = 2;
omega = 2*pi*f;
g = @(t) sin(omega*t);
c = 5;
Nx = 100;
x = linspace(x0,xf,Nx);
dx = (xf-x0)/(Nx-1);
dt = dx/c/1.1246; % stable
%dt = dx/c/1.1245; % unstable
t = t0;
u = zeros(size(x));
plot(x,u);
axis([x0 xf umin umax]);
drawnow;
v = zeros(size(x));
t = t+dt;
v(1) = (4*v(2)-v(3)-2*dx*g(t))/3;
plot(x,v);
axis([x0 xf umin umax]);
drawnow;
while t<tf,
    u(end) = 0;
    u(1) = 2*v(1)-u(1)+(c*dt/dx)^2*(8*v(2)-v(3)-7*v(1)-6*dx*g(t))/2;
    u(2:end-1) = 2*v(2:end-1)-u(2:end-1)+(c*dt/dx)^2*(v(1:end-2)-2*v(2:end-1)+v(3:end));
    plot(x,u);
    axis([x0 xf umin umax]);
    t = t+dt;
    drawnow;
    v(end) = 0;
    v(1) = 2*u(1)-v(1)+(c*dt/dx)^2*(8*u(2)-u(3)-7*u(1)-6*dx*g(t))/2;
    v(2:end-1) = 2*u(2:end-1)-v(2:end-1)+(c*dt/dx)^2*(u(1:end-2)-2*u(2:end-1)+u(3:end));
    plot(x,v);
    axis([x0 xf umin umax]);
    t = t+dt;
    drawnow;
end