function [dtGlobal,dtLocal] = ComputeTimeStep(mesh,prim,CFL)

%% ============================================================
% COMPUTETIMESTEP
%
% Computes local and global CFL time steps.
%
% OUTPUT:
%
% dtLocal(i,j) : Local cell time step
%
% dtGlobal     : Minimum time step in domain
%
%% ============================================================

nxc = mesh.nxc;
nrc = mesh.nrc;

dtLocal = zeros(nxc,nrc);

dtGlobal = 1e20;

for i = 1:nxc

    for j = 1:nrc

        %% ----------------------------------------------------
        % Local grid spacing
        %% ----------------------------------------------------

        dx_local = mesh.dx;

        dr_local = 0.5*( ...
            mesh.dr(i,j) + ...
            mesh.dr(i+1,j));

        %% ----------------------------------------------------
        % Spectral radii
        %% ----------------------------------------------------

        lambda_x = abs(prim.u(i,j)) + prim.a(i,j);

        lambda_r = abs(prim.v(i,j)) + prim.a(i,j);

        %% ----------------------------------------------------
        % Local time step
        %% ----------------------------------------------------

        dtLocal(i,j) = CFL / ...
            (lambda_x/dx_local + lambda_r/dr_local);

        %% ----------------------------------------------------
        % Global time step
        %% ----------------------------------------------------

        dtGlobal = min(dtGlobal,dtLocal(i,j));

    end

end

end