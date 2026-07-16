function viscFlux = ComputeViscousFlux(visc)
%% ============================================================
% COMPUTEVISCOUSFLUX
%
% Computes the viscous flux vectors R and S
% at the East and North faces.
%
%       R = [0 τxx τxr uτxx+vτxr-qx]
%
%       S = [0 τxr τrr uτxr+vτrr-qr]
%% ============================================================

%% Allocate

R = zeros(size(visc.tau_xxE,1),size(visc.tau_xxE,2),4);

S = zeros(size(visc.tau_xxN,1),size(visc.tau_xxN,2),4);

%% ------------------------------------------------------------
% Axial viscous flux (R)
%% ------------------------------------------------------------

R(:,:,1) = 0;

R(:,:,2) = visc.tau_xxE;

R(:,:,3) = visc.tau_xrE;

R(:,:,4) = ...
      visc.face.uE .* visc.tau_xxE ...
    + visc.face.vE .* visc.tau_xrE ...
    - visc.qxE;

%% ------------------------------------------------------------
% Radial viscous flux (S)
%% ------------------------------------------------------------

S(:,:,1) = 0;

S(:,:,2) = visc.tau_xrN;

S(:,:,3) = visc.tau_rrN;

S(:,:,4) = ...
      visc.face.uN .* visc.tau_xrN ...
    + visc.face.vN .* visc.tau_rrN ...
    - visc.qrN;

%% ------------------------------------------------------------
% Store
%% ------------------------------------------------------------

viscFlux.Avf = R;% Axial Viscous Flux

viscFlux.Rvf= S;%Radial Viscous Flux

end