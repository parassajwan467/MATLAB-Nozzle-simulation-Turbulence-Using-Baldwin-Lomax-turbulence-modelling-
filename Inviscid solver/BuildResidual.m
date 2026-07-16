function Residual = BuildResidual(mesh,W,dtLocal)

%% ============================================================
% BUILDRESIDUAL





%% ------------------------------------------------------------
% Primitive Variables
%% ------------------------------------------------------------

prim = PrimitiveFromW(W);


sensor = ComputePressureSensor(mesh,prim);


D = ComputeArtificialDissipation(mesh,W,sensor,dtLocal);

%% ------------------------------------------------------------
%Viscous terms
%% ------------------------------------------------------------

grad = ComputeGradients(mesh,prim);

visc = ComputeViscousProperties(mesh,prim,grad);

viscFlux = ComputeViscousFlux(visc);

%% ------------------------------------------------------------
% Source Term
%% ------------------------------------------------------------

H = zeros(mesh.nxc,mesh.nrc,4);

H(:,:,3)=prim.p + visc.tau_yy_cell;

%% ------------------------------------------------------------
% Allocate Residual
%% ------------------------------------------------------------

Non_Viscous_Residual = zeros(mesh.nxc,mesh.nrc,4);

%Residual = zeros(mesh.nxc,mesh.nrc,4);

%% ============================================================
% Loop over interior cells
%% ============================================================

for i=2:mesh.nxc-1

    for j=2:mesh.nrc-1

        V = mesh.CellArea(i,j);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EAST FACE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        rho_east = 0.5*(prim.rho(i,j)+prim.rho(i+1,j));

        u_east = 0.5*(prim.u(i,j)+prim.u(i+1,j));

        v_east = 0.5*(prim.v(i,j)+prim.v(i+1,j));

        p_east = 0.5*(prim.p(i,j)+prim.p(i+1,j));

        E_east = 0.5*(prim.E(i,j)+prim.E(i+1,j));

        r_east = 0.5*(mesh.rCell(i,j)+mesh.rCell(i+1,j));

        Fx_east = [ ...
            r_east * rho_east * u_east;
            r_east * (rho_east * u_east * u_east + p_east);
            r_east * (rho_east * u_east * v_east);
            r_east * u_east * (rho_east * E_east + p_east)];

        Gx_east  = [ ...
            r_east *rho_east *v_east ;
            r_east *(rho_east *u_east *v_east );
            r_east *(rho_east *v_east *v_east +p_east );
            r_east *v_east *(rho_east *E_east +p_east )];

        FluxE = ...
            Fx_east *mesh.SE_x(i,j)+...
            Gx_east *mesh.SE_r(i,j);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% WEST FACE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        rho_west = 0.5*(prim.rho(i,j)+prim.rho(i-1,j));

        u_west = 0.5*(prim.u(i,j)+prim.u(i-1,j));

        v_west = 0.5*(prim.v(i,j)+prim.v(i-1,j));

        p_west = 0.5*(prim.p(i,j)+prim.p(i-1,j));

        E_west = 0.5*(prim.E(i,j)+prim.E(i-1,j));

        r_west = 0.5*(mesh.rCell(i,j)+mesh.rCell(i-1,j));

        Fx_west = [ ...
            r_west*rho_west*u_west;
            r_west*(rho_west*u_west*u_west+p_west);
            r_west*(rho_west*u_west*v_west);
            r_west*u_west*(rho_west*E_west+p_west)];

        Gx_west = [ ...
            r_west*rho_west*v_west;
            r_west*(rho_west*u_west*v_west);
            r_west*(rho_west*v_west*v_west+p_west);
            r_west*v_west*(rho_west*E_west+p_west)];

        FluxW = ...
            Fx_west*mesh.SW_x(i,j)+...
            Gx_west*mesh.SW_r(i,j);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NORTH FACE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        rho_north = 0.5*(prim.rho(i,j)+prim.rho(i,j+1));

        u_north = 0.5*(prim.u(i,j)+prim.u(i,j+1));

        v_north = 0.5*(prim.v(i,j)+prim.v(i,j+1));

        p_north = 0.5*(prim.p(i,j)+prim.p(i,j+1));

        E_north = 0.5*(prim.E(i,j)+prim.E(i,j+1));

        r_north = 0.5*(mesh.rCell(i,j)+mesh.rCell(i,j+1));

        Fx_north = [ ...
            r_north*rho_north*u_north;
            r_north*(rho_north*u_north*u_north+p_north);
            r_north*(rho_north*u_north*v_north);
            r_north*u_north*(rho_north*E_north+p_north)];

        Gx_north = [ ...
            r_north*rho_north*v_north;
            r_north*(rho_north*u_north*v_north);
            r_north*(rho_north*v_north*v_north+p_north);
            r_north*v_north*(rho_north*E_north+p_north)];

        FluxN = ...
            Fx_north*mesh.SN_x(i,j)+...
            Gx_north*mesh.SN_r(i,j);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SOUTH FACE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        rho_south = 0.5*(prim.rho(i,j)+prim.rho(i,j-1));

        u_south = 0.5*(prim.u(i,j)+prim.u(i,j-1));

        v_south = 0.5*(prim.v(i,j)+prim.v(i,j-1));

        p_south = 0.5*(prim.p(i,j)+prim.p(i,j-1));

        E_south = 0.5*(prim.E(i,j)+prim.E(i,j-1));

        r_south = 0.5*(mesh.rCell(i,j)+mesh.rCell(i,j-1));

        Fx_south = [ ...
            r_south*rho_south*u_south;
            r_south*(rho_south*u_south*u_south+p_south);
            r_south*(rho_south*u_south*v_south);
            r_south*u_south*(rho_south*E_south+p_south)];

        Gx_south = [ ...
            r_south*rho_south*v_south;
            r_south*(rho_south*u_south*v_south);
            r_south*(rho_south*v_south*v_south+p_south);
            r_south*v_south*(rho_south*E_south+p_south)];

        FluxS = ...
            Fx_south*mesh.SS_x(i,j)+...
            Gx_south*mesh.SS_r(i,j);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Residual
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ConvectiveResidual = -(FluxE + FluxW + FluxN + FluxS)/V;

Non_Viscous_Residual(i,j,:) = reshape(ConvectiveResidual + squeeze(H(i,j,:)) ...
    + squeeze(D(i,j,:)), 1,1,4);

        



    end

end
Viscous_Residual = Build_Viscous_Residual(mesh,viscFlux);

Residual = Non_Viscous_Residual + Viscous_Residual;


end