function PostProcess(mesh,W,iso,visc)

%% ============================================================
% POSTPROCESS
%
% Converts conservative variables to primitive variables
% and generates flow-field plots.
%
%% ============================================================

fprintf('========================================\n');
fprintf('Post Processing Results\n');
fprintf('========================================\n');

%% ------------------------------------------------------------
% Primitive variables
%% ------------------------------------------------------------

prim = PrimitiveFromW(W);

%% ------------------------------------------------------------
% Mach Contour
%% ------------------------------------------------------------

figure(10)

contourf(mesh.XC,...
         mesh.YC,...
         prim.Mach,...
         40,...
         'LineStyle','none');

colorbar

xlabel('x (m)')
ylabel('r (m)')

title('Mach Number Contour')

axis equal

%% ------------------------------------------------------------
% Pressure Contour
%% ------------------------------------------------------------

figure(11)

contourf(mesh.XC,...
         mesh.YC,...
         prim.p/1000,...
         40,...
         'LineStyle','none');

colorbar

xlabel('x (m)')
ylabel('r (m)')

title('Pressure (kPa)')

axis equal

%% ------------------------------------------------------------
% Density Contour
%% ------------------------------------------------------------

figure(12)

contourf(mesh.XC,...
         mesh.YC,...
         prim.rho,...
         40,...
         'LineStyle','none');

colorbar

xlabel('x (m)')
ylabel('r (m)')

title('Density (kg/m^3)')

axis equal

%% ------------------------------------------------------------
% Temperature Contour
%% ------------------------------------------------------------

figure(13)

contourf(mesh.XC,...
         mesh.YC,...
         prim.T,...
         40,...
         'LineStyle','none');

colorbar

xlabel('x (m)')
ylabel('r (m)')

title('Temperature (K)')

axis equal

%% ------------------------------------------------------------
% Centerline data
%% ------------------------------------------------------------

jCenter = 1;

xPlot = mesh.XC(:,jCenter);

MachCenter = prim.Mach(:,jCenter);

PressureCenter = prim.p(:,jCenter);

DensityCenter = prim.rho(:,jCenter);

TemperatureCenter = prim.T(:,jCenter);

%% ------------------------------------------------------------
% Mach Comparison
%% ------------------------------------------------------------

figure(20)

plot(xPlot,...
     MachCenter,...
     'b-',...
     'LineWidth',2)

hold on

plot(xPlot,...
     iso.Mach,...
     'ro')

grid on

xlabel('x (m)')
ylabel('Mach')

legend('Numerical','Analytical')

title('Centerline Mach Number')

%% ------------------------------------------------------------
% Pressure Comparison
%% ------------------------------------------------------------

figure(21)

plot(xPlot,...
     PressureCenter/1000,...
     'b-',...
     'LineWidth',2)

hold on

plot(xPlot,...
     iso.P/1000,...
     'ro')

grid on

xlabel('x (m)')
ylabel('Pressure (kPa)')

legend('Numerical','Analytical')

title('Centerline Pressure')

%% ------------------------------------------------------------
% Density Comparison
%% ------------------------------------------------------------

figure(22)

plot(xPlot,...
     DensityCenter,...
     'b-',...
     'LineWidth',2)

hold on

plot(xPlot,...
     iso.rho,...
     'ro')

grid on

xlabel('x (m)')
ylabel('Density (kg/m^3)')

legend('Numerical','Analytical')

title('Centerline Density')

%% ------------------------------------------------------------
% Temperature Comparison
%% ------------------------------------------------------------

figure(23)

plot(xPlot,...
     TemperatureCenter,...
     'b-',...
     'LineWidth',2)

hold on

plot(xPlot,...
     iso.T,...
     'ro')

grid on

xlabel('x (m)')
ylabel('Temperature (K)')

legend('Numerical','Analytical')

title('Centerline Temperature')

%% ------------------------------------------------------------
% Wall velocity
%% ------------------------------------------------------------

figure(24)

plot(mesh.YC(20,:),prim.u(20,:))

hold on



grid on

xlabel('x (m)')
ylabel('velocity (m/s)')

legend('velocity at wall')

title('Wall Velocity')

%% ------------------------------------------------------------
% Exit values
%% ------------------------------------------------------------

figure(25)

tau_w = visc.tau_xrN(:,end);

plot(mesh.XC(:,end),tau_w)


hold on



grid on

xlabel('x (m)')
ylabel('Wall Shear Stress')

legend('Wall Shear Stress')

title('Wall Shear Stress')
%% ------------------------------------------------------------
% Exit values
%% ------------------------------------------------------------

fprintf('\n========== EXIT CONDITIONS ==========\n');

fprintf('Exit Mach        : %.4f\n',...
        mean(prim.Mach(end,:)));

fprintf('Exit Pressure    : %.2f kPa\n',...
        mean(prim.p(end,:))/1000);

fprintf('Exit Temperature : %.2f K\n',...
        mean(prim.T(end,:)));

fprintf('Exit Density     : %.4f kg/m^3\n',...
        mean(prim.rho(end,:)));

fprintf('=====================================\n\n');


end