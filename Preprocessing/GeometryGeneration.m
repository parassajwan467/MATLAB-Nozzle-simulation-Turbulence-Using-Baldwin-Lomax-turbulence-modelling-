function [nozzle_data] = GeometryGeneration()
%% GEOMETRYGENERATION - Generates conical nozzle geometry
%  Outputs:
%    nozzle_data : struct containing nozzle geometry parameters
%                  .rinlet, .rthroat, .Linlet, .Lconv, .L_throat_total
%                  .Ldiv, .L_total, .x_exit, .P_exit, .x_wall, .r_wall
%                  .rexit, .Area, .Area_ratio

fprintf('========================================\n');
fprintf('Generating Nozzle Geometry\n');
fprintf('========================================\n');

mm_to_m = 0.001;

% Geometry parameters in meters
rinlet = 63.5 * mm_to_m;      
rthroat = 20.32 * mm_to_m;    
L_throat_total = 76.2 * mm_to_m;   
Lconv = 43.18 * mm_to_m;      
Ldiv = 38.1 * mm_to_m;        
Lstraight = 22.499 * mm_to_m; 
theta_conv = 45;              
theta_div = 15;               

nozzle_data.rinlet = rinlet;
nozzle_data.rthroat = rthroat;
nozzle_data.Linlet = Lstraight;
nozzle_data.Lconv = Lconv;
nozzle_data.L_throat_total = L_throat_total;
nozzle_data.Ldiv = Ldiv;
nozzle_data.theta_conv = theta_conv;
nozzle_data.theta_div = theta_div;
%nozzle_data.L_total = L_throat_total + Ldiv;
%nozzle_data.x_exit = nozzle_data.L_total;
nozzle_data.P_exit = 78.599e3;  % Exit pressure from paper

% Generate wall profile coordinates
n_straight = 40;
n_conv = 80;
n_div = 60;

% Inlet straight section
x1 = linspace(0, Lstraight, n_straight);
y1 = rinlet * ones(size(x1));

% Inlet fillet
x_centerf1 = 22.499 * mm_to_m;
y_centerf1 = 50.8 * mm_to_m;
radius_f1 = 12.7 * mm_to_m;
theta1 = linspace(pi/2, atan2(59.78*mm_to_m - 50.8*mm_to_m, 31.479*mm_to_m - 22.499*mm_to_m), 40);
x2 = x_centerf1 + radius_f1 * cos(theta1);
y2 = y_centerf1 + radius_f1 * sin(theta1);

% Converging section
x3 = linspace(31.479*mm_to_m, 67.22*mm_to_m, n_conv);
y3 = 59.78*mm_to_m - (x3 - 31.479*mm_to_m);

% Throat fillet
x_centerf2 = 76.2 * mm_to_m;
y_centerf2 = 33.02 * mm_to_m;
theta2 = linspace(atan2(24.04*mm_to_m - 33.02*mm_to_m, 67.22*mm_to_m - 76.2*mm_to_m), ...
                  atan2(20.748*mm_to_m - 33.02*mm_to_m, 79.458*mm_to_m - 76.2*mm_to_m), 40);
x4 = x_centerf2 + radius_f1 * cos(theta2);
y4 = y_centerf2 + radius_f1 * sin(theta2);

% Diverging section
x_fillet_end = 79.458 * mm_to_m;
y_fillet_end = 20.748 * mm_to_m;
x_exit = x_fillet_end + Ldiv * cos(deg2rad(theta_div));
x5 = linspace(x_fillet_end, x_exit, n_div);
y5 = y_fillet_end + (x5 - x_fillet_end) * tan(deg2rad(theta_div));
rexit = y_fillet_end + Ldiv * tan(deg2rad(theta_div));

% Combine all sections
xwall = [x1, x2(2:end-1), x3(2:end-1), x4(2:end-1), x5(2:end)];
ywall = [y1, y2(2:end-1), y3(2:end-1), y4(2:end-1), y5(2:end)];

if length(xwall) ~= length(ywall)
    error('xwall and ywall have different lengths');
end
radius = ywall;
nozzle_data.x_wall = xwall;
nozzle_data.r_wall = ywall;
nozzle_data.rexit = rexit;
nozzle_data.radius= radius;
nozzle_data.Area = pi * radius.^2;
nozzle_data.Area_ratio = nozzle_data.Area / (pi * rthroat^2);
nozzle_data.L_total = x_exit;
nozzle_data.x_exit = nozzle_data.L_total;

% Display geometry summary
fprintf('\n=== NOZZLE GEOMETRY SUMMARY ===\n');
fprintf('Inlet radius:           %.2f mm\n', rinlet/mm_to_m);
fprintf('Throat radius:          %.2f mm\n', rthroat/mm_to_m);
fprintf('Exit radius:            %.2f mm\n', rexit/mm_to_m);
fprintf('Exit pressure:          %.2f kPa\n', nozzle_data.P_exit/1000);
fprintf('Straight inlet length:  %.2f mm\n', Lstraight/mm_to_m);
fprintf('Converging length:      %.2f mm\n', Lconv/mm_to_m);
fprintf('Diverging length:       %.2f mm\n', Ldiv/mm_to_m);
fprintf('Total length:           %.2f mm\n', nozzle_data.L_total/mm_to_m);
fprintf('Area ratio (exit/throat): %.2f\n', nozzle_data.Area_ratio(end));
fprintf('================================\n\n');

% Plot geometry
figure(1);
plot(xwall/mm_to_m, ywall/mm_to_m, 'r-', 'LineWidth', 2);
hold on;
plot(xwall/mm_to_m, -ywall/mm_to_m, 'r-', 'LineWidth', 2);
plot(xwall/mm_to_m, zeros(size(xwall)), 'k--', 'LineWidth', 1);
xlabel('Axial Distance (mm)');
ylabel('Radius (mm)');
title('Conical Nozzle Geometry');
grid on; axis equal; hold off;

end