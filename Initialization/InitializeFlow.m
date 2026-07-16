function flow = InitializeFlow(mesh,iso)
%% ============================================================
% INITIALIZEFLOW
%
% Initializes the flow field using the analytical
% isentropic nozzle solution.
%
% Every axial control-volume column receives the same
% primitive variables obtained from the 1D solution.
%
% Conservative variables are then formed according
% to Mehta (1998).
%
% W(:,:,1) = rho
% W(:,:,2) = rho*u
% W(:,:,3) = rho*v
% W(:,:,4) = rho*E
%% ============================================================

fprintf('========================================\n');
fprintf('Initializing Flow Field\n');
fprintf('========================================\n');

gamma = 1.4;

%% ------------------------------------------------------------
% Allocate primitive variables
%% ------------------------------------------------------------

rho  = zeros(mesh.nxc,mesh.nrc);

u    = zeros(mesh.nxc,mesh.nrc);

v    = zeros(mesh.nxc,mesh.nrc);

p    = zeros(mesh.nxc,mesh.nrc);

T    = zeros(mesh.nxc,mesh.nrc);

Mach = zeros(mesh.nxc,mesh.nrc);

a    = zeros(mesh.nxc,mesh.nrc);

E    = zeros(mesh.nxc,mesh.nrc);

%% ------------------------------------------------------------
% Allocate conservative variables
%% ------------------------------------------------------------

W = zeros(mesh.nxc,mesh.nrc,4);

%% ------------------------------------------------------------
% Copy the analytical solution into every column
%% ------------------------------------------------------------

for i = 1:mesh.nxc

    for j = 1:mesh.nrc

        rho(i,j) = iso.rho(i);

        u(i,j) = iso.u(i);

        v(i,j) = 0.0;

        p(i,j) = iso.P(i);

        T(i,j) = iso.T(i);

        Mach(i,j) = iso.Mach(i);

        a(i,j) = iso.a(i);

    end

end

%% ------------------------------------------------------------
% Compute Total Energy
%% ------------------------------------------------------------

E = E + p ./ ((gamma-1).*rho) + 0.5*(u.^2 + v.^2);

%% ------------------------------------------------------------
% Conservative Variables (Mehta)
%% ------------------------------------------------------------

W(:,:,1) = rho;

W(:,:,2) = rho .* u;

W(:,:,3) = rho .* v;

W(:,:,4) = rho .* E;

%% ============================================================
% STORE EVERYTHING
%% ============================================================

flow.rho = rho;

flow.u = u;

flow.v = v;

flow.p = p;

flow.T = T;

flow.Mach = Mach;

flow.a = a;

flow.E = E;

flow.W = W;
%% ------------------------------------------------------------
% Print Summary
%% ------------------------------------------------------------

fprintf('\n========== INITIAL FLOW ==========\n');

fprintf('Density Range      : %.4f  - %.4f kg/m3\n',...
        min(flow.rho(:)),max(flow.rho(:)));

fprintf('Pressure Range     : %.2f  - %.2f kPa\n',...
        min(flow.p(:))/1000,...
        max(flow.p(:))/1000);

fprintf('Temperature Range  : %.2f  - %.2f K\n',...
        min(flow.T(:)),...
        max(flow.T(:)));

fprintf('Mach Range         : %.3f  - %.3f\n',...
        min(flow.Mach(:)),...
        max(flow.Mach(:)));

fprintf('Velocity Range     : %.2f  - %.2f m/s\n',...
        min(flow.u(:)),...
        max(flow.u(:)));

fprintf('==================================\n\n');

%% ------------------------------------------------------------
% Plot Initial Mach Number
%% ------------------------------------------------------------

figure(4)

contourf(mesh.XC,...
         mesh.YC,...
         flow.Mach,...
         30,...
         'LineStyle','none');

colorbar

xlabel('x (m)')

ylabel('r (m)')

title('Initial Mach Number')

axis equal;
