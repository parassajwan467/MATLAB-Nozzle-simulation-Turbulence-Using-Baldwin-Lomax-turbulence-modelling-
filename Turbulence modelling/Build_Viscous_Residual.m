function Viscous_Residual = Build_Viscous_Residual(mesh,viscFlux)

%% ============================================================
% Viscous_Residual
%
%% ------------------------------------------------------------
% Allocate Residual
%% ------------------------------------------------------------

Viscous_Residual=zeros(mesh.nxc,mesh.nrc,4);

%% ============================================================
% Loop over interior cells
%% ============================================================

for i=2:mesh.nxc-1

    for j=2:mesh.nrc-1

        V = mesh.CellArea(i,j);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VISCOUS EAST FACE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rEast = mesh.YE(i,j);

AvfEast = squeeze(viscFlux.Avf(i,j,:));

RvfEast = squeeze(viscFlux.Rvf(i,j,:));

ViscFluxE = rEast*(AvfEast*mesh.SE_x(i,j) + RvfEast*mesh.SE_r(i,j));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VISCOUS WEST FACE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rWest = mesh.YW(i,j);

AvfWest = squeeze(viscFlux.Avf(i-1,j,:));

RvfWest = squeeze(viscFlux.Rvf(i-1,j,:));

ViscFluxW = rWest*(AvfWest*mesh.SW_x(i,j) ...
    + RvfWest*mesh.SW_r(i,j));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VISCOUS NORTH FACE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rNorth = mesh.YN(i,j);

AvfNorth = squeeze(viscFlux.Avf(i,j,:));

RvfNorth = squeeze(viscFlux.Rvf(i,j,:));

ViscFluxN = rNorth*(AvfNorth*mesh.SN_x(i,j) ...
    + RvfNorth*mesh.SN_r(i,j));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VISCOUS SOUTH FACE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rSouth = mesh.YS(i,j);

AvfSouth = squeeze(viscFlux.Avf(i,j-1,:));

RvfSouth = squeeze(viscFlux.Rvf(i,j-1,:));

ViscFluxS = rSouth*(...
      AvfSouth*mesh.SS_x(i,j) ...
    + RvfSouth*mesh.SS_r(i,j));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Residual
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ViscousResidual = ...
    (ViscFluxE + ViscFluxW + ViscFluxN + ViscFluxS)/V;

Viscous_Residual(i,j,:) = reshape(ViscousResidual,1,1,4);
    end

end


end