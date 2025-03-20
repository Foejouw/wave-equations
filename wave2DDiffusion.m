clear all
close all force

%%

xMin = -1;
xMax = 1.4;
zMin = -1;
zMax = 2;
tMin = 0;
tMax = 1;

Dx = 0.01;
Dz = 0.05;
Dt = 0.0001;


Nx = round((xMax-xMin)/Dx);
Nz = round((zMax-zMin)/Dz);
Nt = round((tMax-tMin)/Dt);
xs = gpuArray(linspace(xMin, xMax, Nx));
zs = gpuArray(linspace(zMin, zMax, Nz));
ts = gpuArray(linspace(tMin, tMax, Nt));
    
V = @(x) 1/2.*cos(pi.*x) + 1;
%V = @(x) exp(-x.^2) + 1;
etaM = 0;
etaV = 0;


vNow = zeros(Nx, Nz);
vNext = zeros(Nx, Nz);
UNow = zeros(Nx, Nz);
UNext = zeros(Nx, Nz);

%%

%v_x_z_t0 = @(x,z) exp(-(z-(zMax-zMin)/2).^2 - (x-(xMax-xMin)/2).^2); % v(x, z, t=0)
v_x_z_t0 = @(x,z) (1+cos(10.*pi.*(z-1/10))).*(z < 2/10) + 0.*x; % v(x, z, t=0)
U_x_z_t0 = @(x,z) (10.*pi.*V(x)*sin(10.*pi.*(z-1/10))) .*(z < 2/10); % v_t(x, z, t=0)

v_x0_z_t = @(z, t) 0; % v(x=xMin, z, t)
v_x1_z_t = @(z, t) 0; % v(x=xMax, z, t)
U_x0_z_t = @(z, t) 0; % v_t(x=xMin, z, t)
U_x1_z_t = @(z, t) 0; % v_t(x=xMax, z, t)

v_x_z0_t = @(x, t) 0; % v(x, z=zMin, t)
v_x_z1_t = @(x, t) 0; % v(x, z=zMax, t)
U_x_z0_t = @(x, t) 0; % v_t(x, z=zMin, t)
U_x_z1_t = @(x, t) 0; % v_t(x, z=zMax, t)

%%


vNow(:, :) = v_x_z_t0(xs', zs);
UNow(:,:) = U_x_z_t0(xs', zs);

r1 = max(Dt.*(V(xs(2:Nx-1))./Dz).^2)
r2 = (etaM + etaV) .* Dt./(Dx.^2)
r3 = (etaM + etaV) .* Dt./(Dz.^2)

%%
w = waitbar(0,'1','Name','Simulating',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

handle_line = surf(zs, xs,vNow(:,:),'LineWidth',1, 'LineStyle','none');
axis([0,1,0,1,-10,10]);
xlabel('z'); ylabel('x'); zlabel('v');
%handle_line = plot(zs, vNow(1,:), 'LineWidth',1);
title('Wave equation');

%%

M(Nt) = struct('cdata',[],'colormap',[]);
%v = VideoWriter("./Movies/alfven2", "MPEG-4");
%open(v);

for it=2:Nt

    if getappdata(w, 'canceling')
        break;
    end

    waitbar(it/Nt,w);
    vNext(2:Nx-1, 2:Nz-1) = vNow(2:Nx-1, 2:Nz-1) + Dt .* UNow(2:Nx-1, 2:Nz-1);
    UNext(2:Nx-1, 2:Nz-1) = UNow(2:Nx-1, 2:Nz-1) + ...
        Dt.*(V(xs(2:Nx-1)')./Dz).^2 .* (vNow(2:Nx-1, 3:Nz) - 2.*vNow(2:Nx-1, 2:Nz-1) + vNow(2:Nx-1, 1:Nz-2)) + ...
        (etaM + etaV) .* Dt./(Dx.^2) .* (UNow(3:Nx, 2:Nz-1) - 2.*UNow(2:Nx-1, 2:Nz-1) + UNow(1:Nx-2, 2:Nz-1)) + ...
        (etaM + etaV) .* Dt./(Dz.^2) .* (UNow(2:Nx-1, 3:Nz) - 2.*UNow(2:Nx-1, 2:Nz-1) + UNow(2:Nx-1, 1:Nz-2));

    vNow = vNext;
    UNow = UNext;

    % BC
    vNow(1, :) = v_x0_z_t(zs, ts(it)); % v(x=xMin, z, t)
    vNow(Nx, :) = v_x1_z_t(zs, ts(it)); % v(x=xMax, z, t)
    UNow(1, :) = U_x0_z_t(zs, ts(it)); % v_t(x=xMin, z, t)
    UNow(Nx, :) = U_x1_z_t(zs, ts(it)); % v_t(x=xMax, z, t)

    
    vNow(:, 1) = v_x_z0_t(xs, ts(it)); % v(x, z=zMin, t)
    vNow(:, Nz) = v_x_z1_t(xs, ts(it)); % v(x, z=zMax, t)
    UNow(:, 1) = U_x_z0_t(xs, ts(it)); % v_t(x, z=zMin, t)
    UNow(:, Nz) = U_x_z1_t(xs, ts(it)); % v_t(x, z=zMax, t)
    
    if mod(it, 5) == 0
        handle_line.ZData = vNow(:,:);
        %handle_line.YData = vNow(2,:);
    
        drawnow limitrate;
        %M(it) = getframe(gcf);
        %writeVideo(v, M(it));
    end
end
delete(w)
%close(v)