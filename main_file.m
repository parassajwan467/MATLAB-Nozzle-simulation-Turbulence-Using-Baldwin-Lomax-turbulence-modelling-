clc;
clear;
close all;

%% ============================================================
% INPUTS
%% ============================================================

P0 = 1.0342e6;
T0 = 555;

gamma = 1.4;

CFL = 0.01;

MaxIter = 5000;

ResidualTolerance = 1e-8;

%% ============================================================
% GEOMETRY
%% ============================================================

nozzle = GeometryGeneration();

%% ============================================================
% GRID
%% ============================================================

mesh = GridGeneration(nozzle);

%% ============================================================
% ISENTROPIC INITIALIZATION
%% ============================================================

iso = IsentropicSolution(mesh,nozzle);

%% ============================================================
% INITIAL FLOW FIELD
%% ============================================================

flow = InitializeFlow(mesh,iso);

W = flow.W;
%fprintf('Check XC(18,31)   = %.4f mm\n', mesh.XC(18,31)*1000);
%fprintf('Check YC(18,31)   = %.4f mm\n', mesh.YC(18,31)*1000);
%fprintf('Check rho(18,31)  = %.4f kg/m3\n', flow.rho(18,31));
%fprintf('Check u(18,31)    = %.4f m/s\n', flow.u(18,31));
%fprintf('Check T(18,31)    = %.4f K\n', flow.T(18,31));
%fprintf('Check p(18,31)    = %.2f Pa\n', flow.p(18,31));

%% ============================================================
% RESIDUAL HISTORY
%% ============================================================

ResidualHistory = zeros(MaxIter,1);

fprintf('\n');
fprintf('Starting Solver...\n\n');

%% ============================================================
% MAIN ITERATION LOOP
%% ============================================================
%% ============================================================
% MAIN ITERATION LOOP
%% ============================================================
%% ============================================================
% MAIN ITERATION LOOP
%% ============================================================

for iter = 1:MaxIter

    %% --------------------------------------------------------
    % Primitive Variables
    %% --------------------------------------------------------

    prim = PrimitiveFromW(W);
    % --- viscous properties for the timestep stability limit ---
    grad = ComputeGradients(mesh,prim);
    visc = ComputeViscousProperties(mesh,prim,grad);
    %% --------------------------------------------------------
    % Time Step
    %% --------------------------------------------------------

    [dt,dtLocal] = ComputeTimeStep(mesh,prim,visc,CFL);

   % W_old = W;

    %% --------------------------------------------------------
    % RK3 Update
    %% --------------------------------------------------------

    W = RK3Step(mesh,W,dtLocal);

   % dW = norm(W(:)-W_old(:));

    %fprintf('Iteration %d   ||ΔW|| = %.12e\n',iter,dW);

    %% --------------------------------------------------------
    % Updated Primitive Variables
    %% --------------------------------------------------------

    prim = PrimitiveFromW(W);

    %% --------------------------------------------------------
    % Check for Solver Failure
    %% --------------------------------------------------------

    if min(prim.rho(:)) <= 0
        error('Negative density detected at iteration %d',iter);
    end

    if min(prim.p(:)) <= 0
        error('Negative pressure detected at iteration %d',iter);
    end

    if any(isnan(W(:))) || any(isinf(W(:)))
        error('NaN or Inf detected at iteration %d',iter);
    end

    %% --------------------------------------------------------
    % Residual Evaluation
    %% --------------------------------------------------------

    R = BuildResidual(mesh,W,dtLocal);

    ResidualHistory(iter) = max(abs(R(:)));

    %% --------------------------------------------------------
    % Monitor Solver
    %% --------------------------------------------------------

    if mod(iter,5)==0

        fprintf('\n');
        fprintf('Iter = %6d\n',iter);
        fprintf('Residual      = %.6e\n',ResidualHistory(iter));
        fprintf('Min Density   = %.6e\n',min(prim.rho(:)));
        fprintf('Min Pressure  = %.6e Pa\n',min(prim.p(:)));
        fprintf('Max Mach      = %.4f\n',max(prim.Mach(:)));

    end

    %% --------------------------------------------------------
    % Convergence Check
    %% --------------------------------------------------------

    if ResidualHistory(iter) < ResidualTolerance

        fprintf('\n');
        fprintf('Converged at iteration %d\n',iter);

        ResidualHistory = ResidualHistory(1:iter);

        break

    end

end



%% ============================================================
% FINAL SOLUTION
%% ============================================================

prim = PrimitiveFromW(W);

grad = ComputeGradients(mesh,prim);

visc = ComputeViscousProperties(mesh,prim,grad);


%% ============================================================
% VISCOSITY DIAGNOSTICS
%% ============================================================

fprintf('\n========== FINAL VISCOSITY CHECK ==========\n');

fprintf('Max molecular viscosity : %e\n', max(visc.muMole(:)));
fprintf('Max turbulent viscosity : %e\n', max(visc.muT(:)));
fprintf('Max total viscosity     : %e\n', max(visc.mu(:)));

fprintf('Max tau_xxE : %e\n', max(abs(visc.tau_xxE(:))));
fprintf('Max tau_xrE : %e\n', max(abs(visc.tau_xrE(:))));
fprintf('Max tau_rrE : %e\n', max(abs(visc.tau_rrE(:))));
fprintf('Max qxE     : %e\n', max(abs(visc.qxE(:))));
fprintf('Max qrE     : %e\n', max(abs(visc.qrE(:))));

%figure

%contourf(mesh.XC,mesh.YC,visc.muT,40,'LineStyle','none')

%colorbar

%xlabel('x (m)')
%ylabel('r (m)')

%title('Turbulent Eddy Viscosity')

%figure

%contourf(mesh.XC,mesh.YC,visc.mu,40,'LineStyle','none')

%colorbar

%xlabel('x (m)')
%ylabel('r (m)')

%title('Total Dynamic Viscosity')

%figure

%plot(mesh.XC(:,end),visc.tau_xrN(:,end),'LineWidth',2)

%grid on

%xlabel('x (m)')
%ylabel('\tau_w (Pa)')

%title('Wall Shear Stress')

%figure

%plot(prim.u(10,:),mesh.YC(10,:),'LineWidth',2)

%grid on

%xlabel('Velocity (m/s)')
%ylabel('Radius (m)')

%title('Velocity Profile near Inlet')

%figure

%plot(prim.u(20,:),mesh.YC(20,:),'LineWidth',2)

%grid on

%xlabel('Velocity (m/s)')
%ylabel('Radius (m)')

%title('Velocity Profile at Throat')

%figure

%plot(prim.u(end-5,:),mesh.YC(end-5,:),'LineWidth',2)

%grid on

%xlabel('Velocity (m/s)')
%ylabel('Radius (m)')

%title('Velocity Profile near Exit')



flowFinal.W = W;

flowFinal.rho = prim.rho;
flowFinal.u   = prim.u;
flowFinal.v   = prim.v;
flowFinal.p   = prim.p;
flowFinal.T   = prim.T;
flowFinal.a   = prim.a;
flowFinal.Mach= prim.Mach;

%% ============================================================
% CONVERGENCE HISTORY
%% ============================================================

figure;

semilogy(ResidualHistory,'LineWidth',2);

xlabel('Iteration');
ylabel('Maximum Residual');

title('Residual Convergence');

grid on;

%% ============================================================
% POST PROCESSING
%% ============================================================
%disp(size(flowFinal));
%disp(class(flowFinal));
PostProcess(mesh,flowFinal.W,iso,visc);

fprintf('\nSimulation Complete\n');
