%% run_Puyang_model.m
% Main script: run the Puyang model and plot each output variable using tight_subplot
clear; clc; close all;

%% 1. Initialize the model
run_struct = Puyang_initialise(0);
pars = run_struct.pars;

%% 2. Define the time span and initial state
tspan = [pars.whenstart, pars.whenend];  % in years (negative values indicate the past)
Y0 = run_struct.startstate;

%% 3. Solve the model
options = odeset('RelTol',1e-6, 'AbsTol',1e-9, 'MaxStep',1e6);
[t_sol, Y_sol] = ode15s(@Puyang_equations, tspan, Y0, options);

%% 4. Time axis (in Ma, preserving negative values)
time_Ma = t_sol/1e6;
disp(['Time vector range: ', num2str(min(time_Ma)), ' to ', num2str(max(time_Ma))]);

%% 5. Calculate output variables (ensure the model solution returns valid data)
pCO2 = (Y_sol(:,1)/pars.A0)*280;         % ppm
O2_sim = Y_sol(:,7);                     % Atmospheric O2
T_sol = pars.T0 + pars.ESS * log2(pCO2/280);  % Temperature (°C)
organicBurial = [diff(Y_sol(:,2)); NaN]*1e-6;
carbonateBurial = [diff(Y_sol(:,3)); NaN]*1e-6;
N_flux = [diff(Y_sol(:,4)); NaN]*1e-6;
model_delta13C = (Y_sol(:,5)./Y_sol(:,2)/pars.R_std - 1)*1000;
model_delta15N = (Y_sol(:,6)./Y_sol(:,4)/pars.R_N_std - 1)*1000;

% Check if output variables are empty or contain NaN
if any(isnan(pCO2))
    warning('pCO2 contains NaN values. Please check the model solution.');
end

%% 6. External forcing signal ( using a sine function)
npts = length(t_sol);
F_long = zeros(npts,1);
F_volc = zeros(npts,1);
for i = 1:npts
    t_abs = -t_sol(i);  % in years, positive indicates the past
    if t_abs >= 232e6 && t_abs <= 234e6
        F_long(i) = pars.A_long * sin(2*pi*(t_abs-232e6)/0.405e6);
        if abs(t_abs-234e6) < pars.delta_t
            F_volc(i) = pars.A_short;
        else
            F_volc(i) = 0;
        end
    elseif t_abs >= 118e6 && t_abs <= 201e6
        F_long(i) = pars.A_long * sin(2*pi*(t_abs-118e6)/8e6);
        Fv = 0;
        volcanic_events = [201e6,200e6; 183.2e6,182.4e6; 167e6,168e6; 160e6,160e6; ...
                           155e6,155e6; 146e6,144e6; 141e6,141e6; 134.4e6,133.2e6; 130.6e6,129.4e6];
        for j = 1:size(volcanic_events,1)
            t_start = volcanic_events(j,1);
            t_end = volcanic_events(j,2);
            if t_abs >= min(t_start,t_end) && t_abs <= max(t_start,t_end)
                Fv = Fv + pars.A_short;
            end
        end
        F_volc(i) = Fv;
    else
        F_long(i) = 0;
        F_volc(i) = 0;
    end
end
F_in = pars.F0 + F_long + F_volc;

%% 7. Calculate yield variation V_prod_stage2 (for the 201–118 Ma stage)
% Ensure that Geo_data_J1.mat is in the path and contains fields 'P' and 'Ba'
s3 = load('Geo_data_J1.mat');
fields3 = fieldnames(s3);
dataJ1 = s3.(fields3{1});
[sortedAgeJ1, idxJ1] = sort(dataJ1.Age, 'ascend');
dataJ1.Age = sortedAgeJ1;
fnamesJ1 = fieldnames(dataJ1);
for i = 1:length(fnamesJ1)
    if ~strcmp(fnamesJ1{i}, 'Age')
        dataJ1.(fnamesJ1{i}) = dataJ1.(fnamesJ1{i})(idxJ1);
    end
end
validJ1 = isfinite(dataJ1.Age) & isfinite(dataJ1.P) & isfinite(dataJ1.Ba);
if sum(validJ1) < 2
    V_prod_stage2 = repmat(mean((dataJ1.P(validJ1)+dataJ1.Ba(validJ1))/2), size(time_Ma));
else
    V_prod_stage2 = interp1(dataJ1.Age(validJ1), (dataJ1.P(validJ1)+dataJ1.Ba(validJ1))/2, time_Ma, 'linear', 'extrap');
end

%% 8. Plot using tight_subplot, with each variable in its own subplot (plotted individually)
figure('Color',[1 1 1]);
ha = tight_subplot(4,2, [0.05 0.05], [0.1 0.1], [0.1 0.05]);

axes(ha(1));
plot(time_Ma, organicBurial, 'g','LineWidth',1.5);
xlabel('Time (Ma)'); ylabel('Organic Burial (mol/yr)');
title('Organic Carbon Burial');

axes(ha(2));
plot(time_Ma, carbonateBurial, 'c','LineWidth',1.5);
xlabel('Time (Ma)'); ylabel('Carbonate Burial (mol/yr)');
title('Carbonate Burial');

axes(ha(3));
plot(time_Ma, N_flux, 'm','LineWidth',1.5);
xlabel('Time (Ma)'); ylabel('Nitrogen Flux (mol/yr)');
title('Nitrogen Flux');

axes(ha(4));
plot(time_Ma, F_volc, 'r','LineWidth',1.5);
xlabel('Time (Ma)'); ylabel('Volcanic Activity (mol/yr)');
title('Volcanic Activity');

axes(ha(5));
plot(time_Ma, pCO2, 'k','LineWidth',1.5);
xlabel('Time (Ma)'); ylabel('CO₂ (ppm)');
title('Atmospheric CO₂ Concentration');

axes(ha(6));
plot(time_Ma, T_sol, 'm--','LineWidth',1.5);
xlabel('Time (Ma)'); ylabel('Temperature (°C)');
title('Global Mean Temperature');

axes(ha(7));
plot(time_Ma, O2_sim, 'g--','LineWidth',1.5);
xlabel('Time (Ma)'); ylabel('Atmospheric O₂ (mol O)');
title('Atmospheric Oxygen Concentration');

axes(ha(8));
plot(time_Ma, F_long, 'b','LineWidth',1.5);
xlabel('Time (Ma)'); ylabel('Astronomical Forcing');
title('Astronomical Forcing Signal');

% If needed, plot the yield variation in a separate figure
figure('Color',[1 1 1]);
plot(time_Ma, V_prod_stage2, 'b--','LineWidth',1.5);
xlabel('Time (Ma)'); ylabel('Yield Variation (dimensionless)');
title('Yield Variation (201–118 Ma Stage)');

fprintf('Plotting completed.\n');
