function visc = ComputeViscousProperties(mesh,prim,grad)
%% ============================================================
% COMPUTEVISCOUSPROPERTIES
%
% Computes
%
% 1. Molecular viscosity (Sutherland)
% 2. Stress tensor
% 3. Heat flux
%
% All quantities are evaluated at East and North faces.
%% ============================================================

gamma = 1.4;

Rgas = 287.0;

Cp = gamma*Rgas/(gamma-1);

Pr = 0.72;

PrT = 0.9;

%% ============================================================
% Sutherland constants (Air) 
% Reference: COMSOL CFD Documentation- Sutherland's Law
%https://doc.comsol.com/5.5/doc/com.comsol.help.cfd/cfd_ug_fluidflow_high_mach.08.27.html
%% ============================================================

mu_ref = 1.716e-5;      % kg/m.s

T_ref = 273.15;         % K

S = 111;              % K

%% ============================================================
% Cell-centre viscosity(Molecular)
%% ============================================================

T = prim.T;

muMole= mu_ref .* (T./T_ref).^(3/2) .* (T_ref+S)./(T+S);

visc.muMole = muMole;

%% ============================================================
% Cell-centre viscosity(Turbulent)
%% ============================================================

turb = ComputeBaldwinLomax(mesh,prim,grad,muMole);

muT = turb.muT;

visc.muT = muT;

%% ============================================================
% Cell-centre viscosity(Total)
%% ============================================================

mu = muMole + muT;

visc.mu = mu;

%% ============================================================
% Face viscosity
%% ============================================================

muE = 0.5*(mu(1:end-1,:) + mu(2:end,:));

muET = 0.5*(muT(1:end-1,:) + muT(2:end,:));

muEMole = 0.5*(muMole(1:end-1,:) + muMole(2:end,:));

muN = 0.5*(mu(:,1:end-1) + mu(:,2:end));

muNT = 0.5*(muT(:,1:end-1) + muT(:,2:end));

muNMole = 0.5*(muMole(:,1:end-1) + muMole(:,2:end));

%% ============================================================
% Face velocities
%% ============================================================

uE = 0.5*(prim.u(1:end-1,:) + prim.u(2:end,:));

vE = 0.5*(prim.v(1:end-1,:) + prim.v(2:end,:));

uN = 0.5*(prim.u(:,1:end-1) + prim.u(:,2:end));

vN = 0.5*(prim.v(:,1:end-1) + prim.v(:,2:end));

%% ============================================================
% Divergence in cylinderical sections is actually ...
... du/dx + dv/dr + v/r...
% The extra v/r is due to the cylindrical geometry. 
%% ============================================================
%% ============================================================
% Divergence at EAST faces
%% ============================================================

YE = mesh.YE(1:end-1,:);

divE = grad.E.dudx + grad.E.dvdr + vE./YE;

%% ============================================================
% Divergence at NORTH faces
%% ============================================================

YN = mesh.YN(:,1:end-1);

divN = grad.N.dudx + grad.N.dvdr + vN./YN;

%% ============================================================
% EAST stresses
%% ============================================================

visc.tau_xxE = -2/3*muE.*divE + 2*muE.*grad.E.dudx;

visc.tau_rrE = -2/3*muE.*divE + 2*muE.*grad.E.dvdr;

visc.tau_xrE = muE.*( grad.E.dudr + grad.E.dvdx );

visc.tau_yyE = -2/3*muE.*divE +2*muE.*vE./YE;% That's tau_theta-theta

%% ============================================================
% NORTH stresses
%% ============================================================

visc.tau_xxN = -2/3*muN.*divN + 2*muN.*grad.N.dudx;

visc.tau_rrN = -2/3*muN.*divN + 2*muN.*grad.N.dvdr;

visc.tau_xrN = muN.*(grad.N.dudr + grad.N.dvdx );

visc.tau_yyN = -2/3*muN.*divN + 2*muN.*vN./YN;

%% ============================================================
% tau_yy_cell
%% ============================================================
tau_yy_cell = zeros(mesh.nxc,mesh.nrc);

for i = 2:mesh.nxc-1
    for j = 2:mesh.nrc-1

        tauEast  = visc.tau_yyE(i,j);

        tauWest  = visc.tau_yyE(i-1,j);

        tauNorth = visc.tau_yyN(i,j);

        tauSouth = visc.tau_yyN(i,j-1);

        tau_yy_cell(i,j) = ...
            0.25*( ...
                tauEast + ...
                tauWest + ...
                tauNorth + ...
                tauSouth );

    end
end

visc.tau_yy_cell = tau_yy_cell;

%% ============================================================
% Thermal conductivity
%% ============================================================

kE = Cp*( muEMole/Pr + muET/PrT );

kN = Cp*( muNMole/Pr + muNT/PrT );

visc.kE = kE;

visc.kN = kN;

%% ============================================================
% Heat flux
%% ============================================================

visc.qxE = -kE .* grad.E.dTdx;

visc.qrE = -kE .* grad.E.dTdr;

visc.qxN = -kN .* grad.N.dTdx;

visc.qrN = -kN .* grad.N.dTdr;

visc.face.uE = uE;

visc.face.uN = uN;

visc.face.vE = vE;

visc.face.vN = vN;

visc.divE = divE;

visc.divN = divN;
end