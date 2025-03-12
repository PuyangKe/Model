function dy = Puyang_equations(t,y)
% Puyang_equations
% 11 state variables (as defined in Table 2 and the extended nitrogen cycle):
% y(1) = A: Atmospheric-ocean CO2 reservoir
% y(2) = G: Organic carbon reservoir
% y(3) = C: Carbonate reservoir
% y(4) = N: Original inorganic nitrogen reservoir (may be included in the original Puyang code; retained here but not used in the extended nitrogen cycle)
% y(5) = X_org: 13C in organic carbon
% y(6) = X_N: 15N in nitrogen
% y(7) = O: Atmospheric oxygen reservoir
% y(8) = N_atm: Atmospheric N2
% y(9) = N_am: Ammonium (NH4+)
% y(10)= N_no3: Nitrate (NO3-)
% y(11)= N_org: Organic nitrogen

global pars
dy = zeros(11,1);

%% --- Original carbon-oxygen cycle part
pCO2 = (y(1)/pars.A0)*280;
dy(1) = -1e-3*y(1);  % A
dy(2) = 1e-3*y(1) - 0.1*y(2);  % G
dy(3) = 0.1*y(2) - 0.05*y(3);  % C
dy(4) = -0.01*y(4);            % N (original, used when the extended nitrogen cycle is not applied)
dy(5) = 0.01*y(2) - 0.001*y(5);  % X_org
dy(6) = 0.005*y(4) - 0.0005*y(6);% X_N
dy(7) = 0.5*y(2) - 0.2*y(7);     % O

%% --- Extended nitrogen cycle part ---
k_fix = pars.k_nfix;       % Nitrogen fixation rate constant
k_nit = pars.k_nit;        % Maximum nitrification rate
K_nit = pars.K_nit;        % Nitrification half-saturation constant
k_assim = pars.k_assim;    % Assimilation rate constant
k_min = pars.k_min;        % Mineralization rate constant
k_denit = pars.k_denit;    % Maximum denitrification rate
K_denit = pars.K_denit;    % Denitrification half-saturation constant

% Flow calculations (units: mol/yr)
F_fix = k_fix * y(8);                          % Nitrogen fixation: Atmospheric N2 → Available inorganic nitrogen
F_nit = k_nit * y(9) / (K_nit + y(9));           % Nitrification: Ammonium → Nitrate
F_assim_am = k_assim * y(9);                     % Assimilation: Ammonium → Organic nitrogen
F_assim_no3 = k_assim * y(10);                   % Assimilation: Nitrate → Organic nitrogen
F_min = k_min * y(11);                           % Mineralization: Organic nitrogen → Ammonium
F_denit = k_denit * y(10) / (K_denit + y(10));     % Denitrification: Nitrate → Atmospheric N2

% Nitrogen cycle mass balance:
dy(8) = F_denit - F_fix;                           % Atmospheric N2 pool
dy(9) = F_fix + F_min - F_nit - F_assim_am;        % Ammonium pool
dy(10)= F_nit - F_denit - F_assim_no3;             % Nitrate pool
dy(11)= F_assim_am + F_assim_no3 - F_min;            % Organic nitrogen pool

end
