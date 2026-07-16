function sensor = ComputePressureSensor(mesh,prim)

%% ============================================================
% COMPUTEPRESSURESENSOR
%
% Jameson-Mehta Pressure Sensor
%
% Computes pressure sensor independently
% on each face of every control volume.
%
% sensor.E
% sensor.W
% sensor.N
% sensor.S
%
%% ============================================================

p = prim.p;

nxc = mesh.nxc;
nrc = mesh.nrc;

tiny = 1.0e-12;


nuX = zeros(nxc,nrc);
nuR = zeros(nxc,nrc);

for i = 2:nxc-1

    for j = 2:nrc-1

        %% ---------------------------------------
        % Axial pressure sensor
        %% ---------------------------------------

        nuX(i,j) = ...
            abs(p(i+1,j)-2*p(i,j)+p(i-1,j)) / ...
           (abs(p(i+1,j)+2*p(i,j)+p(i-1,j))+tiny);

        %% ---------------------------------------
        % Radial pressure sensor
        %% ---------------------------------------

        nuR(i,j) = ...
            abs(p(i,j+1)-2*p(i,j)+p(i,j-1)) / ...
           (abs(p(i,j+1)+2*p(i,j)+p(i,j-1))+tiny);

    end

end

sensor.axial = nuX;
sensor.radial = nuR;
end