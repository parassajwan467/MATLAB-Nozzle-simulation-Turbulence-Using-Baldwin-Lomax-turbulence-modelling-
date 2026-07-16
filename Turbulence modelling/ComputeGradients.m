function grad = ComputeGradients(mesh,prim)
%% ============================================================
% COMPUTEGRADIENTS
%
% Computes:
%
% 1) Cell-centred gradients
% 2) East-face gradients
% 3) North-face gradients
%
% Returned structure:
%
% grad.cell.dudx
% grad.cell.dudr
% grad.cell.dvdx
% grad.cell.dvdr
% grad.cell.dTdx
% grad.cell.dTdr
%
% grad.E.dudx
% grad.E.dudr
% grad.E.dvdx
% grad.E.dvdr
% grad.E.dTdx
% grad.E.dTdr
%
% grad.N.dudx
% grad.N.dudr
% grad.N.dvdx
% grad.N.dvdr
% grad.N.dTdx
% grad.N.dTdr
%% ============================================================

nxc = mesh.nxc;
nrc = mesh.nrc;

X = mesh.XC;
R = mesh.YC;

u = prim.u;
v = prim.v;
T = prim.T;

%% ============================================================
% Allocate cell gradients
%% ============================================================

dudx = zeros(nxc,nrc);
dudr = zeros(nxc,nrc);

dvdx = zeros(nxc,nrc);
dvdr = zeros(nxc,nrc);

dTdx = zeros(nxc,nrc);
dTdr = zeros(nxc,nrc);

%% ============================================================
% CELL-CENTRED X GRADIENTS
%% ============================================================

for j = 1:nrc

    % Left boundary (forward difference)

    dx = X(2,j)-X(1,j);

    dudx(1,j) = (u(2,j)-u(1,j))/dx;
    dvdx(1,j) = (v(2,j)-v(1,j))/dx;
    dTdx(1,j) = (T(2,j)-T(1,j))/dx;

    % Interior (central difference)

    for i = 2:nxc-1

        dx = X(i+1,j)-X(i-1,j);

        dudx(i,j) = (u(i+1,j)-u(i-1,j))/dx;
        dvdx(i,j) = (v(i+1,j)-v(i-1,j))/dx;
        dTdx(i,j) = (T(i+1,j)-T(i-1,j))/dx;

    end

    % Right boundary (backward difference)

    dx = X(nxc,j)-X(nxc-1,j);

    dudx(nxc,j) = (u(nxc,j)-u(nxc-1,j))/dx;
    dvdx(nxc,j) = (v(nxc,j)-v(nxc-1,j))/dx;
    dTdx(nxc,j) = (T(nxc,j)-T(nxc-1,j))/dx;

end

%% ============================================================
% CELL-CENTRED RADIAL GRADIENTS
%% ============================================================

for i = 1:nxc

    % Axis

    dr = R(i,2)-R(i,1);

    dudr(i,1) = (u(i,2)-u(i,1))/dr;
    dvdr(i,1) = (v(i,2)-v(i,1))/dr;
    dTdr(i,1) = (T(i,2)-T(i,1))/dr;

    % Interior

    for j = 2:nrc-1

        dr = R(i,j+1)-R(i,j-1);

        dudr(i,j) = (u(i,j+1)-u(i,j-1))/dr;
        dvdr(i,j) = (v(i,j+1)-v(i,j-1))/dr;
        dTdr(i,j) = (T(i,j+1)-T(i,j-1))/dr;

    end

    % Wall

    dr = R(i,nrc)-R(i,nrc-1);

    dudr(i,nrc) = (u(i,nrc)-u(i,nrc-1))/dr;
    dvdr(i,nrc) = (v(i,nrc)-v(i,nrc-1))/dr;
    dTdr(i,nrc) = (T(i,nrc)-T(i,nrc-1))/dr;

end

%% ============================================================
% Allocate EAST face gradients
%% ============================================================

grad.E.dudx = zeros(nxc-1,nrc);
grad.E.dudr = zeros(nxc-1,nrc);

grad.E.dvdx = zeros(nxc-1,nrc);
grad.E.dvdr = zeros(nxc-1,nrc);

grad.E.dTdx = zeros(nxc-1,nrc);
grad.E.dTdr = zeros(nxc-1,nrc);

%% ============================================================
% EAST faces (average neighbouring cells)
%% ============================================================

for i = 1:nxc-1
    for j = 1:nrc

        grad.E.dudx(i,j) = 0.5*(dudx(i,j)+dudx(i+1,j));
        grad.E.dudr(i,j) = 0.5*(dudr(i,j)+dudr(i+1,j));

        grad.E.dvdx(i,j) = 0.5*(dvdx(i,j)+dvdx(i+1,j));
        grad.E.dvdr(i,j) = 0.5*(dvdr(i,j)+dvdr(i+1,j));

        grad.E.dTdx(i,j) = 0.5*(dTdx(i,j)+dTdx(i+1,j));
        grad.E.dTdr(i,j) = 0.5*(dTdr(i,j)+dTdr(i+1,j));

    end
end

%% ============================================================
% Allocate NORTH face gradients
%% ============================================================

grad.N.dudx = zeros(nxc,nrc-1);
grad.N.dudr = zeros(nxc,nrc-1);

grad.N.dvdx = zeros(nxc,nrc-1);
grad.N.dvdr = zeros(nxc,nrc-1);

grad.N.dTdx = zeros(nxc,nrc-1);
grad.N.dTdr = zeros(nxc,nrc-1);

%% ============================================================
% NORTH faces (average neighbouring cells)
%% ============================================================

for i = 1:nxc
    for j = 1:nrc-1

        grad.N.dudx(i,j) = 0.5*(dudx(i,j)+dudx(i,j+1));
        grad.N.dudr(i,j) = 0.5*(dudr(i,j)+dudr(i,j+1));

        grad.N.dvdx(i,j) = 0.5*(dvdx(i,j)+dvdx(i,j+1));
        grad.N.dvdr(i,j) = 0.5*(dvdr(i,j)+dvdr(i,j+1));

        grad.N.dTdx(i,j) = 0.5*(dTdx(i,j)+dTdx(i,j+1));
        grad.N.dTdr(i,j) = 0.5*(dTdr(i,j)+dTdr(i,j+1));

    end
end

%% ============================================================
% Store cell gradients
%% ============================================================

grad.cell.dudx = dudx;
grad.cell.dudr = dudr;

grad.cell.dvdx = dvdx;
grad.cell.dvdr = dvdr;

grad.cell.dTdx = dTdx;
grad.cell.dTdr = dTdr;

end