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

        %V = mesh.CellArea(i,j);
        V = mesh.rCell(i,j)*mesh.CellArea(i,j);   % was: V = mesh.CellArea(i,j);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VISCOUS EAST FACE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rEast = mesh.YE(i,j);

AvfEast = squeeze(viscFlux.Avfe(i,j,:));

RvfEast = squeeze(viscFlux.Rvfe(i,j,:));

ViscFluxE = rEast*(AvfEast*mesh.SE_x(i,j) + RvfEast*mesh.SE_r(i,j));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VISCOUS WEST FACE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rWest = mesh.YW(i,j);

AvfWest = squeeze(viscFlux.Avfe(i-1,j,:));

RvfWest = squeeze(viscFlux.Rvfe(i-1,j,:));

ViscFluxW = rWest*(AvfWest*mesh.SW_x(i,j) ...
    + RvfWest*mesh.SW_r(i,j));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VISCOUS NORTH FACE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rNorth = mesh.YN(i,j);

AvfNorth = squeeze(viscFlux.Avfn(i,j,:));

RvfNorth = squeeze(viscFlux.Rvfn(i,j,:));

ViscFluxN = rNorth*(AvfNorth*mesh.SN_x(i,j) ...
    + RvfNorth*mesh.SN_r(i,j));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VISCOUS SOUTH FACE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rSouth = mesh.YS(i,j);

AvfSouth = squeeze(viscFlux.Avfn(i,j-1,:));

RvfSouth = squeeze(viscFlux.Rvfn(i,j-1,:));

ViscFluxS = rSouth*(...
      AvfSouth*mesh.SS_x(i,j) ...
    + RvfSouth*mesh.SS_r(i,j));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Residual
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ViscousResidual = ...
    (ViscFluxE + ViscFluxW + ViscFluxN + ViscFluxS)/V;

Viscous_Residual(i,j,:) = reshape(ViscousResidual,1,1,4);

%if i==20 && j==15

   % fprintf('\n===== CELL (20,15) =====\n');

    %disp('East Flux');
    %disp(ViscFluxE');

    %disp('West Flux');
    %disp(ViscFluxW');

    %disp('North Flux');
    %disp(ViscFluxN');

    %disp('South Flux');
    %disp(ViscFluxS');

    %disp('Residual');
    %disp(ViscousResidual');

%end
    end

end


end
