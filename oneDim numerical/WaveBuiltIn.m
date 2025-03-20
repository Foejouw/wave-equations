


c=1;
a=0;
f=0;
m=0.1;

model = createpde;
geometryFromEdges(model, @squareg);
specifyCoefficients(model, m=m, d=0, c=c, a=a, f=f);

applyBoundaryCondition(model, "dirichlet",Edge=[2,4],u=0);
applyBoundaryCondition(model, "neumann", Edge=[1,3], g=0);
generateMesh(model);
pdemesh(model);

u0 = @(location) atan(cos(pi/2*location.x));
ut0 = @(location) 3*sin(pi*location.x).*exp(sin(pi/2*location.y));
setInitialConditions(model,u0,ut0);

n = 31;
tlist = linspace(0,5,n);

result = solvepde(model,tlist);
u = result.NodalSolution;

%%

figure
umax = max(max(u));
umin = min(min(u));
for i = 1:n
    pdeplot(model,XYData=u(:,i),ZData=u(:,i), ...
                  ZStyle="continuous",Mesh="off");
    axis([-1 1 -1 1 umin umax]); 
    xlabel x
    ylabel y
    zlabel u
    M(i) = getframe;
end
%%

figure
umax = max(max(u));
umin = min(min(u));
for i = 1:n
    pdeplot(model,XYData=u(:,1),ZData=u(:,i), ...
                  ZStyle="continuous",Mesh="off");
    axis([-1 1 -1 1 umin umax]); 
    xlabel x
    ylabel y
    zlabel u
    M(i) = getframe;
end

%%
c=1;
a=0;
f=0;
m=0.1;

model = createpde;
geometryFromEdges(model, @squareg);
specifyCoefficients(model, m=m, d=0, c=c, a=a, f=f);

applyBoundaryCondition(model, "neumann", Edge=[2,4], g=0);
applyBoundaryCondition(model, "dirichlet", Edge=[1,3], u=0);
generateMesh(model);
%pdemesh(model);

u0 = @(location) exp(-location.x.^2);
%ut0 = @(location) 3*sin(pi*location.x).*exp(sin(pi/2*location.y));
ut0 = @(location) 0;
setInitialConditions(model,u0,ut0);

n = 31;
tlist = linspace(0,5,n);

result = solvepde(model,tlist);
u = result.NodalSolution;


%%

