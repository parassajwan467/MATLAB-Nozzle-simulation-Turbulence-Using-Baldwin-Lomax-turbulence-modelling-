clc;
clear;
close all;

%% ============================================================
% INPUTS
%% ============================================================

P0 = 1.0342e6;
T0 = 555;

gamma = 1.4;

CFL = 0.01;%try 0.01

MaxIter = 900;

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

    %% --------------------------------------------------------
    % Time Step
    %% --------------------------------------------------------

    [dt,dtLocal] = ComputeTimeStep(mesh,prim,CFL);

    %% --------------------------------------------------------
    % RK3 Update
    %% --------------------------------------------------------

    W = RK3Step(mesh,W,dt,dtLocal);

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

    if mod(iter,50)==0

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
disp(size(flowFinal));
disp(class(flowFinal));
PostProcess(mesh,flowFinal.W,iso,visc);

fprintf('\nSimulation Complete\n');