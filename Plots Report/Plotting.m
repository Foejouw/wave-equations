clear all
close all force

%%
tic


plotting = false;
save = false;
saveName = 'manyPoints';
precision = true;


prec = 10^(-10);
xMin = 0;
xMax = 1;
zMin = 0;
zMax = 10;
tMin = 0;
tMax = 5;

Dx = 0.01;
Dz = 0.01;
Dt = 0.00002;


Nx = round((xMax-xMin)/Dx);
Nz = round((zMax-zMin)/Dz);
Nt = round((tMax-tMin)/Dt);
xs = (linspace(xMin, xMax, Nx));
zs = (linspace(zMin, zMax, Nz));
ts = (linspace(tMin, tMax, Nt));
steep = round(Nx/4);
shallow = round(Nx/2);


V = @(x,z) (1/2.*(cos(2.*pi.*x))) + 1;
%V = @(x, z) 1/2 + exp(-z.^2) + 0.*z;
eta = 5*10^(-4);
%eta = 0;
etaC = 0;
nuV = 0;
N = 10;

vNow = (zeros(Nx, Nz));
vNext = (zeros(Nx, Nz));
UNow = (zeros(Nx, Nz));
UNext = (zeros(Nx, Nz));

steepAmpl = zeros(1, Nt);
shallowAmpl = zeros(1, Nt);

%%
% mu = linspace(0.5001, 0.999, 10) ; 
% 
% vA = V(1, 1) ;
% tauN = 8*10^(-5);
% tauI = 10^(-5);
% Z = 1;
% T = 8000;
% e = -1.602176634*10^(-19);
% sIn = 3.7*10^(-19);
% sEn = 10^(-19);
% mI = 1.602176634*10^(-19);
% mE = 9.1093837139*10^(-31);
% mN = 1.67492749804*10^(-27);
% kB = 1.380649*10^(-23);
% eps0 = 8.8541878188*10^(-12);
% sNN = 2.6*10^(-19) ;
% ksiN = 2-1./mu;
% ksiI = 1./mu - 1;
% c = 300000 ;
% mu0 = 1/(eps0 * c) ;
% 
% nI = 1./mu .* ((10^4 - mu * 10^4)./(2*mu - 1)) ;
% nN = nI + 10^4 ; 
% %nI = nN * ((1-mu)/(2*mu - 1)) ;
% nE = nI;
% 
% s = ksiI./ksiN ;
% nuIn = 4.*nN*sIn.*sqrt(kB*T./(mI*pi)) ;
% nuEn = nN.*sEn.*sqrt(8*T*kB./(pi*mE)) ;
% 
% lambda = 8.48*pi./sqrt(nI) .* (eps0*kB*T./e^2).^(3/2);
% 
% nuII = sqrt(mE./mI) .* nE.*e.^.4.*log(lambda) ./ (3.*mE.^2.*eps0.^2) .* (mE./(2.*pi.*kB.*T)).^(3/2) ;
% 
% nuNN = 4.*nN.*sNN.*sqrt(kB.*T./pi.*mN);
% Delta = 1 - 1./((3.*nuII./nuIn + 4).*(3*nuNN./(s.*nuIn) + 4)) ;
% gamma = 1 + s./(3*nuII./nuIn + 4) ;
% 
% 
% 
% etaA = ksiN.^2 .* vA^2 ./ (nuIn + nuEn);
% 
% eta = mE./(nE.*e^2) .*(  (mN.*nN)./(mN + mE) .* sqrt(8*kB.*T.*(mE+mN)./(pi*mE.*mN)).*sEn + ...
%     mI*3.7*10^(-6)./(mI+mE).*nI/(T.^(3/2) .* Z^2 .*log(lambda) )    ) ;
% 
% om = vA .* sqrt(e^2 .* mu0 .* nI ./ mI) ;
% nuV = nN.*kB .*T .* tauN/2 .* (Delta .* gamma + om.^2.*tauI.^2)/(Delta.^2 + om.^2.*tauI.^2) ;
% 
% 


%%

%v_x_z_t0 = @(x,z) exp(-(z-(zMax-zMin)/2).^2 - (x-(xMax-xMin)/2).^2); % v(x, z, t=0)
%v_x_z_t0 = @(x,z) (1+cos(10.*pi.*(z-1/10))).*(z < 2/10) + (1+cos(10*pi.*(z-3/10))).*((2/10 < z) & (z < 4/10)) + 0.*x; % v(x, z, t=0)
v_x_z_t0 = @(x,z) (1+cos(N.*pi.*(z-1/N))).*(z < 2/N) + 0.*x; % v(x, z, t=0)
U_x_z_t0 = @(x,z) (N.*pi.*V(x, z).*sin(N.*pi.*(z-1/N))) .*(z < 2/N); % v_t(x, z, t=0)
%v_x_z_t0 = @(x,z) 0;
%U_x_z_t0 = @(x,z) 0;

v_x0_z_t = @(z, t) 0; % v(x=xMin, z, t)
v_x1_z_t = @(z, t) 0; % v(x=xMax, z, t)
U_x0_z_t = @(z, t) 0; % v_t(x=xMin, z, t)
U_x1_z_t = @(z, t) 0; % v_t(x=xMax, z, t)

v_x_z0_t = @(x, t) 0;
%v_x_z0_t = @(x, t) 1/2.*sin(10.*t); % v(x, z=zMin, t)
v_x_z1_t = @(x, t) 0; % v(x, z=zMax, t)
%U_x_z0_t = @(x, t) 5.*cos(10.*t); % v_t(x, z=zMin, t)
U_x_z0_t = @(x, t) 0;
U_x_z1_t = @(x, t) 0; % v_t(x, z=zMax, t)

%%
for ix=1:Nx
    for iz=1:Nz
        vNow(ix, iz) = v_x_z_t0(xs(ix), zs(iz));
        UNow(ix,iz) = U_x_z_t0(xs(ix), zs(iz));
    end
end

%steepAmpl = vNow(steep, :);
%shallowAmpl = vNow(shallow, :);
steepAmpl(1) = max(vNow(steep,:));
shallowAmpl(1) = max(vNow(shallow,:));

r1 = max(Dt.*(V(xs', zs)./Dz).^2,[],  'all')
r2 = (eta + nuV) .* Dt./(Dx.^2)
r3 = (etaC) .* Dt./(Dz.^2)
r4 =  (-nuV.*eta).*Dt./(Dx.^4)
r5 = (-nuV.*etaC).*Dt./(Dz.^4)
r6 = (-nuV.*eta - nuV.*etaC) .*Dt.*4./(Dx.^2 .* Dz.^2)
%%

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
% 
% 
if plotting
 g = figure(2);
 ax = gca;
 handle_line_2 = contourf(zs, xs,vNow(:,:), 'LineStyle','-');
 f = figure(1);
 handle_line = surf(zs, xs,vNow(:,:),'LineWidth',1, 'LineStyle','none');
 axis([max(zMin, -1),min(4, zMax),xMin,xMax,-2,2]);
 xlabel('z'); ylabel('x'); zlabel('v');
% %handle_line = plot(zs, vNow(1,:), 'LineWidth',1);
 title('Wave equation');
end

%%
steepAmpl2 = zeros(1, Nt-1);
shallowAmpl2 = zeros(1, Nt-1);
%%
load("PlotAmpl2.mat")
%%

if save
    M(Nt) = struct('cdata',[],'colormap',[]);
    v = VideoWriter("./Movies/"+saveName, "Motion JPEG 2000");
    v.CompressionRatio = 80;
    open(v);
end


for it=2:Nt

    waitbar(it/Nt,w);

    if getappdata(w, 'canceling')
        break;
    end
    %vNow(:, 3) = v_x_z0_t(xs, ts(it)); % v(x, z=zMin, t)
    %UNow(:, 3) = U_x_z0_t(xs, ts(it)); % v_t(x, z=zMin, t)

    %steepAmpl = vNow(steep, :).*(steepAmpl <= vNow(steep, :)) + steepAmpl.*(steepAmpl > vNow(steep, :));
    %shallowAmpl = vNow(shallow, :).*(shallowAmpl <= vNow(shallow, :)) + shallowAmpl.*(shallowAmpl > vNow(shallow, :));

    steepAmpl3(it) = max(vNow(steep, :));
    shallowAmpl3(it) = max(vNow(shallow, :));
    for ix=1:Nx
        
        UNext(ix, 3:Nz-2) = UNow(pos(ix, Nx), 3:Nz-2) + ...
        Dt.*(V(xs(pos(ix, Nx)), zs(3:Nz-2)')'./Dz).^2 .* (vNow(pos(ix, Nx), 4:Nz-1) - 2.*vNow(pos(ix, Nx), 3:Nz-2) + vNow(pos(ix, Nx), 2:Nz-3) ) + ...
        (eta + nuV) .* Dt./(Dx.^2) .* (UNow(pos(ix+1, Nx), 3:Nz-2) - 2.*UNow(pos(ix, Nx), 3:Nz-2) + UNow(pos(ix-1, Nx), 3:Nz-2)) + ...
        (etaC) .* Dt./(Dz.^2) .* (UNow(pos(ix, Nx), 4:Nz-1) - 2.*UNow(pos(ix, Nx), 3:Nz-2) + UNow(pos(ix, Nx), 2:Nz-3)) + ...
        (-nuV.*eta).*Dt./(Dx.^4) .* ( vNow(pos(ix+2, Nx),  3:Nz-2) - 4.*vNow(pos(ix+1, Nx),  3:Nz-2) + 6.*vNow(pos(ix, Nx),  3:Nz-2) - 4.*vNow(pos(ix-1, Nx),  3:Nz-2) + vNow(pos(ix-2, Nx),  3:Nz-2)) + ...
        (-nuV.*etaC).*Dt./(Dz.^4) .* (vNow(pos(ix, Nx), 5:Nz) - 4.*vNow(pos(ix, Nx), 4:Nz-1) + 6.*vNow(pos(ix, Nx), 3:Nz-2) - 4.*vNow(pos(ix, Nx), 2:Nz-3) + vNow(pos(ix, Nx), 1:Nz-4)) + ...
        (-nuV.*eta - nuV.*etaC) .*Dt.*4./(Dx.^2 .* Dz.^2) .* ( ...
            vNow(pos(ix+1, Nx),4:Nz-1) - 2.*vNow(pos(ix, Nx), 4:Nz-1) + vNow(pos(ix-1, Nx), 4:Nz-1) + ...
            vNow(pos(ix+1, Nx), 2:Nz-3) -2.*vNow(pos(ix, Nx), 2:Nz-3) + vNow(pos(ix-1, Nx), 2:Nz-3) + ...
            (-2).*(vNow(pos(ix+1, Nx), 3:Nz-2) - 2.*vNow(pos(ix, Nx),3:Nz-2) + vNow(pos(ix-1, Nx), 3:Nz-2)  ));
    end



    vNext(:, 3:Nz-2) = vNow(:, 3:Nz-2) + Dt .* UNow(:, 3:Nz-2);

    if precision
        vNext(abs(vNext) < prec) = 0;
        UNext(abs(UNext) < prec) = 0;
    end


    vNow = vNext;
    UNow = UNext;

    

    % BC
%     vNow(1, :) = v_x0_z_t(zs, ts(it)); % v(x=xMin, z, t)
%     vNow(Nx, :) = v_x1_z_t(zs, ts(it)); % v(x=xMax, z, t)
%     UNow(1, :) = U_x0_z_t(zs, ts(it)); % v_t(x=xMin, z, t)
%     UNow(Nx, :) = U_x1_z_t(zs, ts(it)); % v_t(x=xMax, z, t)

    
    %vNow(:, 1) = v_x_z0_t(xs, ts(it)); % v(x, z=zMin, t)
    %vNow(:, Nz) = v_x_z1_t(xs, ts(it)); % v(x, z=zMax, t)
    %UNow(:, 1) = U_x_z0_t(xs, ts(it)); % v_t(x, z=zMin, t)
    %UNow(:, Nz) = U_x_z1_t(xs, ts(it)); % v_t(x, z=zMax, t)
    %m(it) = max(vNow(round(Nx/2), :));

    if plotting && mod(it, 10) == 0
        handle_line.ZData = vNow(:,:);
        %handle_line.YData = vNow(2,:);
        drawnow limitrate;
    if save && mod(it, 10) == 0
        M(it) = getframe(gcf);
        writeVideo(v, M(it));
    end
    end
    if plotting && mod(it, 10) == 0
        drawnow
        contourf(ax, zs, xs, vNow, 'LineStyle','none')

    end
end
delete(w)

if save
    close(v);
end


toc
%%
load("PlotAmpl3.mat", "shallowAmpl2", "steepAmpl2")
%%

%%
steepAmpl = cat(2, steepAmpl, steepAmpl2);
shallowAmpl = cat(2, shallowAmpl, shallowAmpl2);
%%
ts = linspace(tMin, tMax.*2, Nt.*2);
loglog(ts, steepAmpl, 'red', 'LineWidth',1.2)
hold on
loglog(ts, shallowAmpl, 'blue', 'LineWidth',1.2)
hold on
loglog(ts, ts.^(-3/2), 'LineWidth',1.2)

legend(gca, {strcat('$\max_{z} (b(',num2str(xs(steep)),', z, t))$'), strcat('$\max_{z} (b(',num2str(xs(shallow)),', z, t))$'), '$t^{-3/2}$'}, "Interpreter","latex", 'Location','southwest', 'FontSize',11)
xlim([1/2 8] )
xlabel('t')
ylabel('Maximum amplitude')

%%

%%

function ix = pos(ix, Nx)
    if ix <= 0
       %ix=1;
       ix = Nx+ix;
    elseif ix > Nx
       ix=ix-Nx;
       %ix=Nx;
    end
end

