



Dx = 0.1;
Dt = 0.01;
Nx = 100;
Nt = 100;
a = 1;
freq = 10;

xMin = 0;
xMax = Nx*Dx;
tMin = 0;
tMax = Nt*Dt;

xs = linspace(xMin, xMax, Nx);
ts = linspace(tMin, tMax, Nt);
A = zeros(Nt, Nx);

% Initial condition: u(t=0,x) = 0
% Already satisfied by initialisation of A
for i=1:Nx
    A(1, i) = exp(-(xs(i) - xMax/2)^2);
    %A(1, i) = 0;
end

% Initial condition: u_t(t=0,x) = 0
for i=2:Nx-1
    A(2,i) = 2*A(1,i)+ 1/2 * (a * Dt/Dx)^2*(A(1,i+1) - 2*A(1,i) + A(1,i-1))-2*(xs(i)-xMax/2)*exp(-(xs(i)-xMax/2)^2) ;
    %-4*(xs(i)-xMax/2)*exp(-(xs(i)-xMax/2)^2)
end

% Boundary condition: u(t, x=xMin) = ...
for j=1:Nt
    A(j,1) = 0;
    %A(j, 1) = 3*sin(2*pi*ts(j)*freq)*exp(-ts(j));
end

% Boundary condition: u(t, x=xMax) = 0
for j=1:Nt
    A(j,Nx) = 0;
end

(a * Dt/Dx)^2


% Simulation
for j=2:(Nt-1)
    for i=2:Nx-1
        %A(j+1,i) = 2*A(j,i) - A(makeTimePos(j-1),i) + (a*Dt/Dx)^2*(A(j,makeNeg(i+1, xMax)) - 2*A(j,i) + A(j, makePos(i-1, xMax)));
        %A(j+1,i) = A(j, i) - Dt/Dx*a*(A(j,makeNeg(i+1, xMax)) - A(j,makePos(i-1, xMax)));
        A(j+1,i) = 2*A(j,i) - A(j-1,i) + (a*Dt/Dx)^2*(A(j,i+1) - 2*A(j,i) + A(j, i-1));
    end
end


for j=1:5:Nt
    plot(xs, A(j,:), '.-')
    hold on
end

%%

h = figure;
ax = gca;
xlabel('z')
ylabel('v')
ylim([-15 15]);
ax.NextPlot = 'replaceChildren';
M(Nt) = struct('cdata',[],'colormap',[]);
h.Visible = 'on';
v = VideoWriter("Test13", "MPEG-4");
open(v);
for it=1:Nt
    plot(xs, A(it, :));
    drawnow limitrate;
    M(it) = getframe(h);
    %writeVideo(v, M(it));
end
%open(v)
%writeVideo(v, M(:));
close(v)


