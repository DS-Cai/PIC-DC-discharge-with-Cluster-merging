function [efield,pot]=solve_Poisson(rho,Time,DX,EPSILON0,N_G,VOLTAGE,OMEGA,INV_DX)
% solve Poisson equation (Thomas algorithm)
A = 1.0;
B = -2.0;
C = 1.0;
S = 1.0 / (2.0 * DX);
ALPHA = -DX * DX / EPSILON0;
g = zeros(1, N_G);
w = zeros(1, N_G);
f = zeros(1, N_G);
pot= zeros(1, N_G);
efield= zeros(1, N_G);

% Apply potential to the electrodes - boundary conditions
pot(1) = VOLTAGE;    % Potential at the powered electrode
pot(N_G) = 0.0;                        % Potential at the grounded electrode

% Solve Poisson equation-TDMA
f(2:N_G-1) = ALPHA * rho(2:N_G-1);
f(2) = f(2) - pot(1);
f(N_G-1) = f(N_G-1) - pot(N_G);

w(2) = C / B;
g(2) = f(2) / B;

for i = 3:N_G-1
    w(i) = C / (B - A * w(i-1));
    g(i) = (f(i) - A * g(i-1)) / (B - A * w(i-1));
end

pot(N_G-1) = g(N_G-1);

for i = N_G-2:-1:2
    pot(i) = g(i) - w(i) * pot(i+1);
end

% Compute electric field
for i = 2:N_G-1
    efield(i) = (pot(i-1) - pot(i+1)) * S;
end

efield(1) = (pot(1) - pot(2)) * INV_DX - rho(1) * DX / (2.0 * EPSILON0);       % Powered electrode
efield(N_G) = (pot(N_G-1) - pot(N_G)) * INV_DX + rho(N_G) * DX / (2.0 * EPSILON0);   % Grounded electrode