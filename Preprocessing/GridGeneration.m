function mesh = GridGeneration(nozzle_data)

fprintf('========================================\n');
fprintf('Generating Computational Grid\n');
fprintf('========================================\n');

%% ============================================================
% USER INPUT
%% ============================================================

nx = 42;
nr = 32;

beta = 2.0;

mm_to_m = 0.001;

%% ============================================================
% AXIAL GRID
%% ============================================================

x = linspace(0,nozzle_data.L_total,nx);

dx = x(2)-x(1);

%% ============================================================
% WALL PROFILE
%% ============================================================

wallRadius = interp1( ...
    nozzle_data.x_wall,...
    nozzle_data.r_wall,...
    x,...
    'pchip');

%% ============================================================
% STRETCHED RADIAL COORDINATE
%% ============================================================

eta = linspace(0,1,nr);

stretch = 1 - ...
    (exp(beta*(1-eta))-1) ./ ...
    (exp(beta)-1);

%% ============================================================
% NODE COORDINATES
%% ============================================================

X = zeros(nx,nr);
Y = zeros(nx,nr);

for i = 1:nx

    X(i,:) = x(i);

    Y(i,:) = wallRadius(i)*stretch;

end

%% ============================================================
% REFERENCE RADIAL VECTOR
%% ============================================================

r = Y(1,:);

%% ============================================================
% CELL NUMBERS
%% ============================================================

nxc = nx-1;
nrc = nr-1;

%% ============================================================
% CELL CENTRES
%% ============================================================

XC = zeros(nxc,nrc);
YC = zeros(nxc,nrc);

for i = 1:nxc

    for j = 1:nrc

        XC(i,j) = 0.25*( X(i,j) + X(i+1,j) + X(i,j+1) + X(i+1,j+1));

        YC(i,j) = 0.25*( Y(i,j) + Y(i+1,j) + Y(i,j+1) + Y(i+1,j+1));

    end

end

%% ============================================================
% LOCAL RADIAL SPACING
%% ============================================================

dr = zeros(nx,nr-1);

for i = 1:nx

    for j = 1:nr-1

        dr(i,j) = Y(i,j+1)-Y(i,j);

    end

end

%% ============================================================
% CONTROL VOLUME AREAS
%% ============================================================

CellArea = zeros(nxc,nrc);

for i = 1:nxc

    for j = 1:nrc

        local_dr = 0.5*( ...
            dr(i,j) + ...
            dr(i+1,j));

        CellArea(i,j) = dx*local_dr;

    end

end

%% ============================================================
% WALL GEOMETRY
%% ============================================================

wallSlope = gradient(wallRadius,x);

wallAngle = atan(wallSlope);

%% ============================================================
% CELL-CENTRE WALL RADIUS
%% ============================================================

wallRadiusCell = 0.5*(wallRadius(1:end-1) + wallRadius(2:end));

%% ============================================================
% CELL-CENTRE RADIUS
%% ============================================================

rCell = YC;

%% ============================================================
% FACE CENTERS
%
% Computes the centre coordinates of the four faces
% of every control volume.
%
% East
% West
% North
% South
%% ============================================================
XE = zeros( nxc, nrc);
YE = zeros( nxc, nrc);

XW = zeros( nxc, nrc);
YW = zeros( nxc, nrc);

XN = zeros( nxc, nrc);
YN = zeros( nxc, nrc);

XS = zeros( nxc, nrc);
YS = zeros( nxc, nrc);

for i = 1: nxc
    for j = 1: nrc

        % East face
        XE(i,j) = 0.5*(X(i+1,j)+X(i+1,j+1));
        YE(i,j) = 0.5*(Y(i+1,j)+Y(i+1,j+1));

        % West face
        XW(i,j) = 0.5*(X(i,j)+X(i,j+1));
        YW(i,j) = 0.5*(Y(i,j)+Y(i,j+1));

        % North face
        XN(i,j) = 0.5*(X(i,j+1)+X(i+1,j+1));
        YN(i,j) = 0.5*(Y(i,j+1)+Y(i+1,j+1));

        % South face
        XS(i,j) = 0.5*(X(i,j)+X(i+1,j));
        YS(i,j) = 0.5*(Y(i,j)+Y(i+1,j));

    end
end

%% ============================================================
% FACE VECTORS AND LENGTHS
%
% Computes the outward face vectors
% and corresponding face lengths for
% every control volume.
%% ============================================================

SE_x = zeros(nxc,nrc);
SE_r = zeros(nxc,nrc);

SW_x = zeros(nxc,nrc);
SW_r = zeros(nxc,nrc);

SN_x = zeros(nxc,nrc);
SN_r = zeros(nxc,nrc);

SS_x = zeros(nxc,nrc);
SS_r = zeros(nxc,nrc);

LE = zeros(nxc,nrc);
LW = zeros(nxc,nrc);
LN = zeros(nxc,nrc);
LS = zeros(nxc,nrc);



for i = 1: nxc
    for j = 1: nrc
     dxEast = X(i+1,j+1)-X(i+1,j);

     drEast = Y(i+1,j+1)-Y(i+1,j);

     dxNorth = X(i+1,j+1) - X(i,j+1);
     drNorth = Y(i+1,j+1) - Y(i,j+1);

     %%-----EAST FACE-----

     SE_x(i,j)= drEast;

     SE_r(i,j)=-dxEast;

     LE(i,j)=sqrt(SE_x(i,j)^2+SE_r(i,j)^2);

     %%-----WEST FACE-----

     SW_x(i,j) = -SE_x(i,j);
     SW_r(i,j) = -SE_r(i,j);

     LW(i,j) = LE(i,j);

     %%-----NORTH FACE-----

     SN_x(i,j) = -drNorth;
     SN_r(i,j) = dxNorth;

     LN(i,j) = sqrt(SN_x(i,j)^2 + SN_r(i,j)^2);

     %%-----SOUTH FACE-----
     SS_x(i,j) = -SN_x(i,j);
     SS_r(i,j) = -SN_r(i,j);

     LS(i,j) = LN(i,j);

    end
end    
%% ============================================================
% GRID SUMMARY
%% ============================================================

fprintf('\n========== GRID SUMMARY ==========\n');

fprintf('Axial nodes             : %d\n',nx);

fprintf('Radial nodes            : %d\n',nr);

fprintf('Axial cells             : %d\n',nxc);

fprintf('Radial cells            : %d\n',nrc);

fprintf('Total control volumes   : %d\n',nxc*nrc);

fprintf('Domain length           : %.2f mm\n',...
    nozzle_data.L_total/mm_to_m);

fprintf('Average dx              : %.4f mm\n',...
    dx/mm_to_m);

fprintf('Average dr              : %.4f mm\n',...
    mean(dr(:))/mm_to_m);

fprintf('Stretching factor       : %.2f\n',beta);

fprintf('==================================\n\n');

%% ============================================================
% GRID PLOT
%% ============================================================

figure(2);
clf;
hold on;

skip_x = max(1,round(nx/30));
skip_r = max(1,round(nr/15));

for i = 1:skip_x:nx

    plot(X(i,:)/mm_to_m,...
         Y(i,:)/mm_to_m,...
         'b-');

end

for j = 1:skip_r:nr

    plot(X(:,j)/mm_to_m,...
         Y(:,j)/mm_to_m,...
         'b-');

end

plot(nozzle_data.x_wall/mm_to_m,...
     nozzle_data.r_wall/mm_to_m,...
     'r','LineWidth',2);

plot(nozzle_data.x_wall/mm_to_m,...
    -nozzle_data.r_wall/mm_to_m,...
     'r','LineWidth',2);

plot(nozzle_data.x_wall/mm_to_m,...
     zeros(size(nozzle_data.x_wall)),...
     'k--');

plot(XC/mm_to_m,...
     YC/mm_to_m,...
     'k.',...
     'MarkerSize',6);

xlabel('Axial Distance (mm)');
ylabel('Radius (mm)');

title(sprintf('Structured Grid (%d x %d Nodes)',...
      nx,nr));

axis equal;
grid on;
hold off;

%% ============================================================
% STORE EVERYTHING IN GRID STRUCTURE
%% ============================================================

mesh.nx = nx;
mesh.nr = nr;

mesh.nxc = nxc;
mesh.nrc = nrc;

mesh.beta = beta;

mesh.x = x;
mesh.r = r;

mesh.X = X;
mesh.Y = Y;

mesh.XC = XC;
mesh.YC = YC;

mesh.dx = dx;
mesh.dr = dr;

mesh.CellArea = CellArea;

mesh.wallRadius = wallRadius;
mesh.wallSlope = wallSlope;
mesh.wallAngle = wallAngle;
mesh.wallRadiusNode = wallRadius;
mesh.wallRadiusCell = wallRadiusCell;

mesh.rCell = rCell;

mesh.XE = XE;
mesh.YE = YE;

mesh.XW = XW;
mesh.YW = YW;

mesh.XN = XN;
mesh.YN = YN;

mesh.XS = XS;
mesh.YS = YS;

mesh.SE_x = SE_x;
mesh.SE_r = SE_r;

mesh.SW_x = SW_x;
mesh.SW_r = SW_r;

mesh.SN_x = SN_x;
mesh.SN_r = SN_r;

mesh.SS_x = SS_x;
mesh.SS_r = SS_r;

mesh.LE = LE;
mesh.LW = LW;
mesh.LN = LN;
mesh.LS = LS;
end