clear all
close all force

%%
tic
xMin = -1;
xMax = 1;
zMin = 0;
zMax = 2;
tMin = 0;
tMax = 2;

Dx = 0.02;
Dz = 0.05;
Dt = 0.0004;


Nx = round((xMax-xMin)/Dx);
Nz = round((zMax-zMin)/Dz);
Nt = round((tMax-tMin)/Dt);
xs = (linspace(xMin, xMax, Nx));
zs = (linspace(zMin, zMax, Nz));
ts = (linspace(tMin, tMax, Nt));
    
%V = @(x) 1/2.*(cos(pi.*x)) + 1;
V = @(x) 1/2 + exp(-x.^2);
eta = 5.*10.^(-4);
etaC = 1/1000;
nuV = 3/1000;
N = 10;

vNow = zeros(Nx, Nz);
vNext = zeros(Nx, Nz);
UNow = zeros(Nx, Nz);
UNext = zeros(Nx, Nz);

%%

%v_x_z_t0 = @(x,z) exp(-(z-(zMax-zMin)/2).^2 - (x-(xMax-xMin)/2).^2); % v(x, z, t=0)
%v_x_z_t0 = @(x,z) (1+cos(10.*pi.*(z-1/10))).*(z < 2/10) + (1+cos(10*pi.*(z-3/10))).*((2/10 < z) & (z < 4/10)) + 0.*x; % v(x, z, t=0)
v_x_z_t0 = @(x,z) (1+cos(N.*pi.*(z-1/N))).*(z < 2/N) + 0.*x; % v(x, z, t=0)
U_x_z_t0 = @(x,z) (N.*pi.*V(x).*sin(N.*pi.*(z-1/N))) .*(z < 2/N); % v_t(x, z, t=0)
%U_x_z_t0 = @(x,z)  V(x) .* v_x_z_t0(x, z);

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

r1 = max(Dt.*(V(xs')./Dz).^2)
r2 = (eta + nuV) .* Dt./(Dx.^2)
r3 = (etaC) .* Dt./(Dz.^2)
r4 =  (-nuV.*eta).*Dt./(Dx.^4)
r5 = (-nuV.*etaC).*Dt./(Dz.^4)
r6 = (-nuV.*eta - nuV.*etaC) .*Dt.*4./(Dx.^2 .* Dz.^2)


%%
% 
% variable = {'x'; 'z'; 't'; 'extra info'};
% Number = [Nx;Nz;Nt; 0];
% StepSize = [Dx;Dz;Dt; 0];
% Range = [xMin xMax;zMin zMax;tMin tMax; [0 0]];
% V_at_Min = {v_x0_z_t; v_x_z0_t; v_x_z_t0; 'na'};
% V_at_Max = {v_x1_z_t; v_x_z1_t; 'na'; 'na'};
% Derivative_at_Min = {U_x0_z_t; U_x_z0_t; U_x_z_t0; 'na'};
% Derivative_at_Max = {U_x1_z_t; U_x_z1_t; 'na'; 'na'};
% AlfvenSpeed = {'na'; 'na'; 'na'; V};
% eta_etaC_nuV = {'na'; 'na'; 'na'; [eta etaC nuV]};
% Constants = {'na'; 'na'; N; 'na'};
% 
% T = table(variable, Number, StepSize, Range, V_at_Min, V_at_Max, Derivative_at_Min, Derivative_at_Max, AlfvenSpeed, eta_etaC_nuV, Constants);

%%

%F = uifigure("Position",[200 200 1200 200]);
%uitable(F, "Data",T, "Position",[1 1 1200 200])

%%
w = waitbar(0,'1','Name','Simulating',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

handle_line = surf(zs, xs,vNow(:,:),'LineWidth',1, 'LineStyle','none');
axis([0,min(4, zMax),xMin,xMax,-2,2]);
xlabel('z'); ylabel('x'); zlabel('v');
%handle_line = plot(zs, vNow(1,:), 'LineWidth',1);
title('Wave equation');

%%

%M(Nt) = struct('cdata',[],'colormap',[]);
%v = VideoWriter("./Movies/alfvenHarm4", "MPEG-4");
%open(v);

for it=2:Nt

    if getappdata(w, 'canceling')
        break;
    end

    vNext(:, 3:Nz-2) = vNow(:, 3:Nz-2) + Dt .* UNow(:, 3:Nz-2);

    for ix=1:Nx
        
        UNext(ix, 3:Nz-2) = UNow(pos(ix, Nx), 3:Nz-2) + ...
        Dt.*(V(xs(pos(ix, Nx))')./Dz).^2 .* (vNow(pos(ix, Nx), 4:Nz-1) - 2.*vNow(pos(ix, Nx), 3:Nz-2) + vNow(pos(ix, Nx), 2:Nz-3) ) + ...
        (eta + nuV) .* Dt./(Dx.^2) .* (UNow(pos(ix+1, Nx), 3:Nz-2) - 2.*UNow(pos(ix, Nx), 3:Nz-2) + UNow(pos(ix-1, Nx), 3:Nz-2)) + ...
        (etaC) .* Dt./(Dz.^2) .* (UNow(pos(ix, Nx), 4:Nz-1) - 2.*UNow(pos(ix, Nx), 3:Nz-2) + UNow(pos(ix, Nx), 2:Nz-3)) + ...
        (-nuV.*eta).*Dt./(Dx.^4) .* ( vNow(pos(ix+2, Nx),  3:Nz-2) - 4.*vNow(pos(ix+1, Nx),  3:Nz-2) + 6.*vNow(pos(ix, Nx),  3:Nz-2) - 4.*vNow(pos(ix-1, Nx),  3:Nz-2) + vNow(pos(ix-2, Nx),  3:Nz-2)) + ...
        (-nuV.*etaC).*Dt./(Dz.^4) .* (vNow(pos(ix, Nx), 5:Nz) - 4.*vNow(pos(ix, Nx), 4:Nz-1) + 6.*vNow(pos(ix, Nx), 3:Nz-2) - 4.*vNow(pos(ix, Nx), 2:Nz-3) + vNow(pos(ix, Nx), 1:Nz-4)) + ...
        (-nuV.*eta - nuV.*etaC) .*Dt.*4./(Dx.^2 .* Dz.^2) .* ( ...
            vNow(pos(ix+1, Nx),4:Nz-1) - 2.*vNow(pos(ix, Nx), 4:Nz-1) + vNow(pos(ix-1, Nx), 4:Nz-1) + ...
            vNow(pos(ix+1, Nx), 2:Nz-3) -2.*vNow(pos(ix, Nx), 2:Nz-3) + vNow(pos(ix-1, Nx), 2:Nz-3) + ...
            (-2).*(vNow(pos(ix+1, Nx), 3:Nz-2) - 2.*vNow(pos(ix, Nx),3:Nz-2) + vNow(pos(ix-1, Nx), 3:Nz-2)  ));
    end

    vNow = vNext;
    UNow = UNext;

    % BC
    %vNow(1, :) = v_x0_z_t(zs, ts(it)); % v(x=xMin, z, t)
    %vNow(Nx, :) = v_x1_z_t(zs, ts(it)); % v(x=xMax, z, t)
    %UNow(1, :) = U_x0_z_t(zs, ts(it)); % v_t(x=xMin, z, t)
    %UNow(Nx, :) = U_x1_z_t(zs, ts(it)); % v_t(x=xMax, z, t)

    
    %vNow(:, 1) = v_x_z0_t(xs, ts(it)); % v(x, z=zMin, t)
    %vNow(:, Nz) = v_x_z1_t(xs, ts(it)); % v(x, z=zMax, t)
    %UNow(:, 1) = U_x_z0_t(xs, ts(it)); % v_t(x, z=zMin, t)
    %UNow(:, Nz) = U_x_z1_t(xs, ts(it)); % v_t(x, z=zMax, t)
    %m(it) = max(vNow(round(Nx/2), :));

    if mod(it, 50) == 0
        waitbar(it/Nt,w);
        handle_line.ZData = vNow(:,:);
        %handle_line.YData = vNow(2,:);
    
        %drawnow limitrate;
        %M(it) = getframe(gcf);
        %writeVideo(v, M(it));
    end
    if mod(it, 200) == 0
        drawnow
    end
end
delete(w)
%close(v)

%%
toc








