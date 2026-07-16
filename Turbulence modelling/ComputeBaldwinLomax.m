function turb = ComputeBaldwinLomax(mesh,prim,grad,muMole)

%% ============================================================
% COMPUTEBALDWINLOMAX
%
% Computes:
%
% 1. Wall distance
% 2. Vorticity magnitude
% 3. Friction velocity
% 4. y+
% 5. Inner layer eddy viscosity
% 6. Outer layer eddy viscosity
% 7. Final turbulent viscosity
%
%% ============================================================
kappa = 0.4;
Aplus = 26;
K = 0.0168;
Ccp = 1.6;
Cwk = 0.25;

%% ============================================================
%Allocating variables
%% ============================================================

wallDist = zeros(mesh.nxc,mesh.nrc);

%omega = zeros(mesh.nxc,mesh.nrc);

utau = zeros(mesh.nxc,1);

yplus = zeros(mesh.nxc,mesh.nrc);

D1 = zeros(mesh.nxc,mesh.nrc);

lm = zeros(mesh.nxc,mesh.nrc);

muInner = zeros(mesh.nxc,mesh.nrc);

muOuter = zeros(mesh.nxc,mesh.nrc);


for i=1:mesh.nxc

    wallRadius = mesh.wallRadius(i);

    for j=1:mesh.nrc

        wallDist(i,j) = wallRadius - mesh.YC(i,j);

    end

end

omega = abs(grad.cell.dudr - grad.cell.dvdx);

for i=1:mesh.nxc

    rhoWall = prim.rho(i,end);
    muWall  = muMole(i,end);

    tauw = muWall * abs(grad.cell.dudr(i,end));

    utau(i) = sqrt(tauw/rhoWall);

    yplus(i,:) = rhoWall .* utau(i) .* wallDist(i,:) ./ muWall;

    D1(i,:) = 1-exp(-yplus(i,:)/Aplus);

    lm(i,:) = kappa .* wallDist(i,:) .* D1(i,:);

    muInner(i,:) = ...
        prim.rho(i,:) .* lm(i,:).^2 .* omega(i,:);

    F = wallDist(i,:) .* D1(i,:) .* omega(i,:);

    Fmax = max(F);

    Fmax = max(Fmax,1e-12)

    ymax = max(wallDist(i,:));

    Umax = max(sqrt(prim.u(i,:).^2 + prim.v(i,:).^2));

    Fw = min(ymax*Fmax,Cwk*ymax*Umax^2/Fmax);

    FKIF = 1 ./ ...
        (1 + 5.5*(0.3.*wallDist(i,:)./ymax).^6);

    muOuter(i,:) = ...
        prim.rho(i,:) .* ...
        K .* Ccp .* Fw .* FKIF;
end

muT = min(muInner,muOuter);
%% ============================================================
% Store
%% ============================================================

turb.wallDist = wallDist;

turb.omega = omega;

turb.utau = utau;

turb.yplus = yplus;

turb.muInner = muInner;

turb.muOuter = muOuter;

turb.muT = muT;

end