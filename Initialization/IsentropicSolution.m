function iso = IsentropicSolution(mesh,nozzle_data)
%% ============================================================
% ISENTROPICSOLUTION
%
% Computes the analytical 1D isentropic solution corresponding
% to the nozzle area distribution.
%
% Returns the flow variables at each axial control volume.
%
%==============================================================

fprintf('========================================\n');
fprintf('Computing Isentropic Initial Solution\n');
fprintf('========================================\n');

%% ------------------------------------------------------------
% Gas properties
%% ------------------------------------------------------------

gamma = 1.4;
R = 287;

P0 = 1.0342e6;      % Pa
T0 = 555;           % K

%% ------------------------------------------------------------
% Geometry
%% ------------------------------------------------------------

nxc = mesh.nxc;

At = pi*nozzle_data.rthroat^2;

Area = pi*mesh.wallRadiusCell.^2;

AreaRatio = Area/At;

%% ------------------------------------------------------------
% Allocate memory
%% ------------------------------------------------------------

Mach = zeros(nxc,1);

Pressure = zeros(nxc,1);

Temperature = zeros(nxc,1);

Density = zeros(nxc,1);

Velocity = zeros(nxc,1);

SoundSpeed = zeros(nxc,1);

Energy = zeros(nxc,1);

%% ------------------------------------------------------------
% Locate throat cell
%% ------------------------------------------------------------

[~,throatCell] = min(abs(mesh.XC(:,1)-nozzle_data.L_throat_total));

fprintf('Throat located at cell %d of %d\n',throatCell,nxc);

%% ------------------------------------------------------------
% Loop over every axial cell
%% ------------------------------------------------------------

for i=1:nxc

    AR = AreaRatio(i);

    if i < throatCell

        Mach(i)=AreaMachSolver(AR,gamma,'subsonic');

    elseif i==throatCell

        Mach(i)=1.0;

    else

        Mach(i)=AreaMachSolver(AR,gamma,'supersonic');

    end

    %% Isentropic relations

    Temperature(i)=T0/(1+((gamma-1)/2)*Mach(i)^2);

    Pressure(i)=P0*(Temperature(i)/T0)^(gamma/(gamma-1));

    Density(i)=Pressure(i)/(R*Temperature(i));

    SoundSpeed(i)=sqrt(gamma*R*Temperature(i));

    Velocity(i)=Mach(i)*SoundSpeed(i);

    cv = R/(gamma-1);

    Energy(i) = cv*Temperature(i) + 0.5*Velocity(i)^2;

end

%% ------------------------------------------------------------
% Store in structure
%% ------------------------------------------------------------

iso.Area = Area;

iso.AreaRatio = AreaRatio;

iso.Mach = Mach;

iso.P = Pressure;

iso.T = Temperature;

iso.rho = Density;

iso.u = Velocity;

iso.a = SoundSpeed;

iso.E = Energy;

%% ------------------------------------------------------------
% Print summary
%% ------------------------------------------------------------

fprintf('\n========== ISENTROPIC SOLUTION ==========\n');

fprintf('Inlet Mach      : %.4f\n',Mach(1));

fprintf('Throat Mach     : %.4f\n',Mach(throatCell));

fprintf('Exit Mach       : %.4f\n',Mach(end));

fprintf('Exit Pressure   : %.2f kPa\n',Pressure(end)/1000);

fprintf('Exit Temp       : %.2f K\n',Temperature(end));

fprintf('Exit Velocity   : %.2f m/s\n',Velocity(end));

fprintf('=========================================\n\n');

%% ------------------------------------------------------------
% Plot analytical solution
%% ------------------------------------------------------------

figure(3)

subplot(3,2,1)

plot(mesh.XC(:,1),Mach,'LineWidth',2)

xlabel('x (m)')

ylabel('Mach')

grid on

subplot(3,2,2)

plot(mesh.XC(:,1),Pressure/1000,'LineWidth',2)

xlabel('x (m)')

ylabel('Pressure (kPa)')

grid on

subplot(3,2,3)

plot(mesh.XC(:,1),Temperature,'LineWidth',2)

xlabel('x (m)')

ylabel('Temperature (K)')

grid on

subplot(3,2,4)

plot(mesh.XC(:,1),Density,'LineWidth',2)

xlabel('x (m)')

ylabel('Density (kg/m^3)')

grid on

subplot(3,2,5)

plot(mesh.XC(:,1),Velocity,'LineWidth',2)

xlabel('x (m)')

ylabel('Velocity (m/s)')

grid on

sgtitle('Analytical Isentropic Solution')

end