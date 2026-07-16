function D = ComputeArtificialDissipation(mesh,W,sensor,dtLocal)

%% ============================================================
% COMPUTEARTIFICIALDISSIPATION
%
% Jameson-Mehta Artificial Dissipation
%
% Eq. (26)-(29) of Mehta (1998)
%
%% ============================================================

k2 = 1.0;
k4 = 0.03125;

nxc = mesh.nxc;
nrc = mesh.nrc;

D = zeros(nxc,nrc,4);



%% ============================================================
% Interior Cells
%% ============================================================

for i = 3:nxc-2

    for j = 3:nrc-2

        %% ----------------------------------------------------
        % Cell Areas
        %% ----------------------------------------------------

        Aij = mesh.CellArea(i,j);

        AE = mesh.CellArea(i+1,j);
        AW = mesh.CellArea(i-1,j);

        AN = mesh.CellArea(i,j+1);
        AS = mesh.CellArea(i,j-1);

        %% ----------------------------------------------------
        % Mehta face scaling (Eq. 27)
        %% ----------------------------------------------------

        scaleE = ...
            dtLocal(i,j)/(2*Aij) + ...
            dtLocal(i+1,j)/(2*AE);

        scaleW = ...
            dtLocal(i,j)/(2*Aij) + ...
            dtLocal(i-1,j)/(2*AW);

        scaleN = ...
            dtLocal(i,j)/(2*Aij) + ...
            dtLocal(i,j+1)/(2*AN);

        scaleS = ...
            dtLocal(i,j)/(2*Aij) + ...
            dtLocal(i,j-1)/(2*AS);

        %% ====================================================
        % Conservative Variables
        %% ====================================================

        for k = 1:4

            %%=================================================
            %% EAST FACE
            %%=================================================

            eps2 = k2 * max(sensor.axial(i,j),sensor.axial(i+1,j));

            eps4 = max(0.0,k4-eps2);

            d2E = eps2 * ...
                (W(i+1,j,k)-W(i,j,k));

            d4E = eps4 * ...
                (W(i+2,j,k) ...
               -3*W(i+1,j,k) ...
               +3*W(i,j,k) ...
               -W(i-1,j,k));

            dE = scaleE*(d2E-d4E);

            %%=================================================
            %% WEST FACE
            %%=================================================

            eps2 = k2 * max(sensor.axial(i,j),sensor.axial(i-1,j));

            eps4 = max(0.0,k4-eps2);

            d2W = eps2 * ...
                (W(i,j,k)-W(i-1,j,k));

            d4W = eps4 * ...
                (W(i+1,j,k) ...
               -3*W(i,j,k) ...
               +3*W(i-1,j,k) ...
               -W(i-2,j,k));

            dW = scaleW*(d2W-d4W);

            %%=================================================
            %% NORTH FACE
            %%=================================================

            eps2 = k2 * max(sensor.radial(i,j),sensor.radial(i,j+1));

            eps4 = max(0.0,k4-eps2);

            d2N = eps2 * ...
                (W(i,j+1,k)-W(i,j,k));

            d4N = eps4 * ...
                (W(i,j+2,k) ...
               -3*W(i,j+1,k) ...
               +3*W(i,j,k) ...
               -W(i,j-1,k));

            dN = scaleN*(d2N-d4N);

            %%=================================================
            %% SOUTH FACE
            %%=================================================

            eps2 = k2 * max(sensor.radial(i,j),sensor.radial(i,j-1));

            eps4 = max(0.0,k4-eps2);

            d2S = eps2 * ...
                (W(i,j,k)-W(i,j-1,k));

            d4S = eps4 * ...
                (W(i,j+1,k) ...
               -3*W(i,j,k) ...
               +3*W(i,j-1,k) ...
               -W(i,j-2,k));

            dS = scaleS*(d2S-d4S);

            %%=================================================
            %% Cell Dissipation (Eq. 26)
            %%=================================================

            D(i,j,k) = ...
                ( dE ...
                - dW ...
                + dN ...
                - dS ) / Aij;

        end

    end

end

end