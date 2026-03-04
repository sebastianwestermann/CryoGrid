

%% Example 1: Simple A → B reaction
T = 298.15;           % Temperature in K (25°C)
c_r = 0.1;            % [A] = 0.1 M
c_p = 0.01;           % [B] = 0.01 M
G0 = -28700;          % Standard Gibbs free energy in J/mol (-28.7 kJ/mol)

dG = gibbsFreeEnergy(T, c_r, c_p, G0);

%% Example 2: N2 + 3H2 ⇌ 2NH3 with stoichiometry
T = 500;              % 500 K
c_r = [0.5; 1.5];     % [N2] = 0.5 M, [H2] = 1.5 M
c_p = 0.2;            % [NH3] = 0.2 M
G0 = -32800;          % J/mol
nu_r = [1; 3];        % 1 N2 + 3 H2
nu_p = 2;             % 2 NH3

dG = gibbsFreeEnergy(T, c_r, c_p, G0, nu_r, nu_p);

%% Example 3: Equilibrium check (Q = K)
% When Q = K, dG should be 0
K_eq = 1e5;
T = 298.15;
G0 = -8.314 * T * log(K_eq);  % dG0 = -RT ln(K)
c_r = 1/sqrt(K_eq);   % Set concentrations so Q = K
c_p = sqrt(K_eq);

dG = gibbsFreeEnergy(T, c_r, c_p, G0);



function dG = gibbsFreeEnergy(T, c_r, c_p, G0, nu_r, nu_p, units)
% GIBBSFREEENERGY Computes Gibbs free energy change for a chemical reaction
%   dG = gibbsFreeEnergy(T, c_r, c_p, G0) computes Delta G for a 1:1 reaction
%   dG = gibbsFreeEnergy(T, c_r, c_p, G0, nu_r, nu_p) specifies stoichiometry
%
%   Inputs:
%       T    - Temperature in Kelvin (K)
%       c_r  - Concentrations of reactants (cell array or vector)
%       c_p  - Concentrations of products (cell array or vector)
%       G0   - Standard Gibbs free energy change (J/mol or kJ/mol)
%       nu_r - Stoichiometric coefficients for reactants (optional, default: 1)
%       nu_p - Stoichiometric coefficients for products (optional, default: 1)
%       units- 'J' or 'kJ' for output (optional, default: same as G0)
%
%   Output:
%       dG   - Gibbs free energy change for the reaction
%
%   Formula: Delta G = Delta G^0 + RT * ln(Q)
%            where Q = (product of c_p^nu_p) / (product of c_r^nu_r)

    % Gas constant
    R = 8.314462618;  % J/(mol*K)

    % Set default stoichiometry if not provided
    if nargin < 5 || isempty(nu_r)
        nu_r = ones(size(c_r));
    end
    if nargin < 6 || isempty(nu_p)
        nu_p = ones(size(c_p));
    end
    
    % Ensure column vectors
    c_r = c_r(:);
    c_p = c_p(:);
    nu_r = nu_r(:);
    nu_p = nu_p(:);
    
    % Standard concentration (1 M)
    c0 = 1;
    
    % Calculate reaction quotient Q
    % Q = prod(c_p .^ nu_p) / prod(c_r .^ nu_r)
    Q_numerator = prod((c_p ./ c0) .^ nu_p);
    Q_denominator = prod((c_r ./ c0) .^ nu_r);
    Q = Q_numerator / Q_denominator;
    
    % Calculate Gibbs free energy: dG = G0 + RT*ln(Q)
    dG = G0 + R * T * log(Q);
    
    % Convert units if requested
    if nargin >= 7
        if strcmpi(units, 'kJ') && abs(G0) > 1000
            dG = dG / 1000;
        elseif strcmpi(units, 'J') && abs(G0) < 100
            dG = dG * 1000;
        end
    end
    
    % Display results
    fprintf('Reaction Quotient Q = %.6f\n', Q);
    fprintf('RT ln(Q) = %.3f J/mol\n', R*T*log(Q));
    fprintf('Delta G = %.3f J/mol (%.3f kJ/mol)\n', dG, dG/1000);
    
    % Determine spontaneity
    if dG < 0
        fprintf('Reaction is spontaneous under these conditions (dG < 0)\n');
    elseif dG > 0
        fprintf('Reaction is non-spontaneous under these conditions (dG > 0)\n');
    else
        fprintf('System is at equilibrium (dG = 0)\n');
    end
end
