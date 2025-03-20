clear all

tic

Dx = 0.05;
Dz = 0.01;
Dt = 0.0001;
Nx = 2/Dx;
Nz = 2/Dz;
Nt = 5000;
a = 1/2;
freq = 5;
eta = 5*10^(-4);
n = 10;

xMin = 0;
xMax = Nx*Dx;
zMin=0;
zMax = Nz*Dz;
tMin = 0;
tMax = Nt*Dt;

xs = linspace(xMin, xMax, Nx);
zs = linspace(zMin, zMax, Nz);
ts = linspace(tMin, tMax, Nt);
%U = zeros(Nx, Nz, Nt);
b = (zeros(Nx, Nz, Nt));
Va = @(x) 1/2*cos(pi*x) + 1; 





% BC: b(x, z=0,t) = 0, ok

% Initial condition: b(x,z,t=0)
for iz=1:Nz
    if 0 < zs(iz) && zs(iz) < 2/n
        b(:, iz, 1) = -1 - cos(n*pi*(zs(iz) - 1/n));
    elseif 2/n <= zs(iz) && zs(iz) < 4/n
        b(:, iz, 1) = 1 + cos(n*pi*(zs(iz) - 3/n));
    end
end

% IC: -V * b_z(x,z,t=0)
for iz=1:Nz-1
    if 0 < zs(iz) && zs(iz) < 2/n
        for ix=1:Nx
            b(ix,iz+1,1) = ( - Dz * n*pi* Va(xs(ix)) * sin(n*pi*(zs(iz)-1/n)) + b(ix, iz,1));
        end
    elseif 2/n <= zs(iz) && zs(iz) < 4/n
        for ix=1:Nx
            b(ix,iz+1,1) = ( Dz * n*pi* Va(xs(ix)) * sin(n*pi*(zs(iz)-3/n)) + b(ix, iz,1));
        end
    else
        b(:, iz+1, 1) = 0;
    end
end

% IC: b_t(x, z, t=0)
for iz=1:Nz
    if 0 < zs(iz) && zs(iz) < 2/n
        for ix=1:Nx
            b(ix,iz,2) = (-Dt*n*pi*Va(xs(ix))*sin(n*pi*(zs(iz)-1/n)) + b(ix,iz,1));
        end
    elseif 2/n <= zs(iz) && zs(iz) < 4/n
        for ix=1:Nx
            b(ix,iz,2) = (Dt*n*pi*Va(xs(ix))*sin(n*pi*(zs(iz)-3/n)) + b(ix,iz,1));
        end
    else
        b(:, iz, 2) = 0;
    end
end




r = 1/(1 - 3*eta*Dt/Dx^2)
r2 = (max(Va(xs)) * Dt/Dz)^2
r3 = eta*Dt/Dx^2

% Simulation
for it=2:(Nt-1)
    for ix=1:(Nx)
        if ix==1 % BC: b_x(x=0, z, t) = 0
            for iz=2:(Nz-1)
                b(ix,iz,it+1) = (r * ...
                (2*b(ix,iz,it) - b(ix,iz,it-1) + Va(xs(ix))^2*(Dt/Dz)^2 * ...
                (b(ix,iz+1,it) + b(ix,iz-1,it) - 2*b(ix,iz,it)    )   + ...
                eta*Dt/Dx^2 * ( b(ix+1,iz,it) + b(ix+1,iz,it) -2*b(ix,iz,it)  -b(ix+1,iz,it-1) - ...
                b(ix+1,iz,it-1) - b(ix,iz,it-1)  )    ));
            end
        elseif ix==Nx % BC: n_x(x=1, z, t) = 0
            for iz=2:(Nz-1)
                b(ix,iz,it+1) = (r * ...
                (2*b(ix,iz,it) - b(ix,iz,it-1) + Va(xs(ix))^2*(Dt/Dz)^2 * ...
                (b(ix,iz+1,it) + b(ix,iz-1,it) - 2*b(ix,iz,it)    )   + ...
                eta*Dt/Dx^2 * ( 2*b(ix-1,iz,it)  -2*b(ix,iz,it)  -2*b(ix-1,iz,it-1) + ...
                 - b(ix,iz,it-1)  )    ));
            end
        else
            for iz=2:(Nz-1)
                b(ix,iz,it+1) = (r * ...
                (2*b(ix,iz,it) - b(ix,iz,it-1) + Va(xs(ix))^2*(Dt/Dz)^2 * ...
                (b(ix,iz+1,it) + b(ix,iz-1,it) - 2*b(ix,iz,it)    )   + ...
                eta*Dt/Dx^2 * ( b(ix+1,iz,it) + b(ix-1,iz,it) -2*b(ix,iz,it)  -b(ix+1,iz,it-1) - ...
                b(ix-1,iz,it-1) - b(ix,iz,it-1)  )    ));
            end
        end
    end
end

toc

%%
%save("magn", "b");

%%

ax = axes();
hold(ax, 'on')
for it=Nt/2:1000:Nt
    plot(zs, b(Nx/2, :, it), "DisplayName",[num2str(ts(it))])
    hold on
end
%hold(ax, 'off')
legend()


%%

surf(zs(1:100), xs, b(:,1:100,Nt))

xlabel('z')
ylabel('x')
zlabel('b')

%%

h = figure;
ax = gca;
xlabel('z')
ylabel('x')
zlabel('b')
%zlim([-3 3]);
ax.NextPlot = 'replaceChildren';
M(Nt) = struct('cdata',[],'colormap',[]);
h.Visible = 'on';
v = VideoWriter("Test13", "MPEG-4");
open(v);
for it=1:30:Nt
    surf(zs(1:round(Nz/4)), xs, b(:,1:round(Nz/4),it), 'LineStyle', 'none');
    drawnow limitrate;
    M(it) = getframe(h);
    %writeVideo(v, M(it));
end
%open(v)
%writeVideo(v, M(:));
close(v)

%%

h.Visible = 'on';
movie(M);

%%

reader = VideoReader("Test.mp4");
while hasFrame(reader)
    imshow(readFrame(reader));
end

clear reader