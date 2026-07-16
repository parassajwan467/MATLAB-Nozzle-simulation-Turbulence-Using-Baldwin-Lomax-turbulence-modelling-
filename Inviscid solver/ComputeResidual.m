function Residual = ComputeResidual(mesh,F,G,H,D)

%% ============================================================
% COMPUTERESIDUAL
%
% Finite Volume Residual
%
% Residual =
%
% -(Flux Balance)/CellArea
% + Source Term
% + Artificial Dissipation
%
%% ============================================================

nxc = mesh.nxc;
nrc = mesh.nrc;

Residual = zeros(nxc,nrc,4);

for i = 2:nxc-1

    for j = 2:nrc-1

      V = mesh.rCell(i,j)*mesh.CellArea(i,j);

        for k = 1:4

            %% ----------------------------------------
            % Face Fluxes
            %% ----------------------------------------

            FW = 0.5*(F(i,j,k)+F(i-1,j,k));
            FE = 0.5*(F(i,j,k)+F(i+1,j,k));

            GW = 0.5*(G(i,j,k)+G(i-1,j,k));
            GE = 0.5*(G(i,j,k)+G(i+1,j,k));

            FN = 0.5*(F(i,j,k)+F(i,j+1,k));
            FS = 0.5*(F(i,j,k)+F(i,j-1,k));

            GN = 0.5*(G(i,j,k)+G(i,j+1,k));
            GS = 0.5*(G(i,j,k)+G(i,j-1,k));
            %% ----------------------------------------
            % Net Flux Through Cell
            %% ----------------------------------------

            FluxBalance = FE*mesh.SE_x(i,j) +  GE*mesh.SE_r(i,j) + ...
                          FW*mesh.SW_x(i,j) + GW*mesh.SW_r(i,j) + ...
                          FN*mesh.SN_x(i,j) + GN*mesh.SN_r(i,j) + ...
                          FS*mesh.SS_x(i,j) + GS*mesh.SS_r(i,j);

            %% ----------------------------------------
            % Residual
            %% ----------------------------------------

            Residual(i,j,k) = -(FluxBalance/V) + H(i,j,k) + D(i,j,k);

        end

    end

end

end