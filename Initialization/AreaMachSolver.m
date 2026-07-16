function Mach = AreaMachSolver(AR,gamma,branch)
%==================================================================
% Solves the Area-Mach relation using the bisection method.
%
% INPUTS
%   AR      : Area ratio (A/A*)
%   gamma   : Ratio of specific heats
%   branch  : 'subsonic' or 'supersonic'
%
% OUTPUT
%   Mach
%==================================================================

tol     = 1e-10;
maxIter = 200;

%% ---------------------------------------------------------------
% Physical checks
%% ---------------------------------------------------------------

if AR < 1

    error('Area ratio must be >= 1.');

end

if abs(AR-1) < 1e-12

    Mach = 1.0;
    return

end

%% ---------------------------------------------------------------
% Initial search interval
%% ---------------------------------------------------------------

switch lower(branch)

    case 'subsonic'

        Mlow  = 1e-8;
        Mhigh = 0.999999;

    case 'supersonic'

        Mlow  = 1.000001;
        Mhigh = 5.0;

    otherwise

        error('Unknown branch.')

end

%% ---------------------------------------------------------------
% Check if interval brackets the root
%% ---------------------------------------------------------------

fLow  = AreaFunction(Mlow,gamma)-AR;
fHigh = AreaFunction(Mhigh,gamma)-AR;

if fLow*fHigh > 0

    error('Initial interval does not bracket a root.')

end

%% ---------------------------------------------------------------
% Bisection iterations
%% ---------------------------------------------------------------

for iter = 1:maxIter

    Mmid = 0.5*(Mlow+Mhigh);

    fMid = AreaFunction(Mmid,gamma)-AR;

    if abs(fMid) < tol

        Mach = Mmid;
        return

    end

    if fLow*fMid < 0

        Mhigh = Mmid;
        fHigh = fMid;

    else

        Mlow = Mmid;
        fLow = fMid;

    end

end

%% ---------------------------------------------------------------
% Final answer
%% ---------------------------------------------------------------

Mach = 0.5*(Mlow+Mhigh);

end