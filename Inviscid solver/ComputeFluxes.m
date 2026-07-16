function [F,G] = ComputeFluxes(mesh,prim)
%% ============================================================
% COMPUTEFLUXES
%
% Computes axisymmetric Euler flux vectors
% for every control volume.
%
% Mehta (1998) formulation:
%
% W = [ rho
%       rho*u
%       rho*v
%       rho*E ]
%
% F = x-direction flux
% G = r-direction flux
%
% Fluxes include axisymmetric radius factor r.
%% ============================================================

%gamma = 1.4;

nxc = mesh.nxc;
nrc = mesh.nrc;

F = zeros(nxc,nrc,4);
G = zeros(nxc,nrc,4);

for i = 1:nxc

    for j = 1:nrc

        %% ----------------------------------------------------
        % Primitive Variables
        %% ----------------------------------------------------

        rho = prim.rho(i,j);

        u = prim.u(i,j);

        v = prim.v(i,j);

        p = prim.p(i,j);

        E = prim.E(i,j);

        %% ----------------------------------------------------
        % Axisymmetric Radius
        %% ----------------------------------------------------

        r = mesh.rCell(i,j);

        %% ----------------------------------------------------
        % Conservative Variables
        %% ----------------------------------------------------

        rhou = rho*u;

        rhov = rho*v;

        rhoE = rho*E;

        %% ----------------------------------------------------
        % X-FLUX
        %% ----------------------------------------------------

        F(i,j,1) = r*rhou;

        F(i,j,2) = r*(rhou*u + p);

        F(i,j,3) = r*(rhou*v);

        F(i,j,4) = r*u*(rhoE + p);

        %% ----------------------------------------------------
        % R-FLUX
        %% ----------------------------------------------------

        G(i,j,1) = r*rhov;

        G(i,j,2) = r*(rhov*u);

        G(i,j,3) = r*(rhov*v + p);

        G(i,j,4) = r*v*(rhoE + p);

    end

end

end