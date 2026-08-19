function [dtGlobal,dtLocal] = ComputeTimeStep(mesh,prim,visc,CFL)

nxc = mesh.nxc;
nrc = mesh.nrc;

dtLocal = zeros(nxc,nrc);
dtGlobal = 1e20;

for i = 1:nxc
    for j = 1:nrc

        dx_local = mesh.dx;
        dr_local = 0.5*(mesh.dr(i,j) + mesh.dr(i+1,j));

        lambda_x = abs(prim.u(i,j)) + prim.a(i,j);
        lambda_r = abs(prim.v(i,j)) + prim.a(i,j);

        % Diffusive (viscous) stability term
        nu_eff = visc.mu(i,j) / prim.rho(i,j);   % effective kinematic viscosity (molecular+turbulent)
        lambda_visc = 2*nu_eff*(1/dx_local^2 + 1/dr_local^2);

        dtLocal(i,j) = CFL / ...
            (lambda_x/dx_local + lambda_r/dr_local + lambda_visc);

        dtGlobal = min(dtGlobal,dtLocal(i,j));

    end
end

end
