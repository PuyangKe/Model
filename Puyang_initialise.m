function run = Puyang_initialise(runcontrol)
    % Puyang_initialise
    % Initializes the Puyang model. All parameters, state variables, and external forcing data strictly follow the SCION definition.
    % Time is represented using negative numbers (e.g., -240e6 represents 240 million years ago).
    % runcontrol: If >=1, sensitivity perturbations are enabled (set to 0 in this example).

    %% Clear old variables
    clear global pars forcings workingstate state gridstate INTERPSTACK sensanal sensparams
    global pars

    %% 1. Initial reservoir states (refer to Table 2)
    pars.A0 = 3.193e18;    % Atmosphere-ocean CO2 reservoir [mol C]
    pars.G0 = 1.25e21;     % Buried organic C [mol C]
    pars.C0 = 5.0e21;      % Buried carbonate C [mol C]
    pars.S0 = 4e19;        % Ocean sulfate [mol S]
    pars.PYR0 = 1.8e20;    % Buried pyrite S [mol S]
    pars.GYP0 = 2.0e20;    % Buried gypsum S [mol S]
    pars.P0 = 3.1e15;      % Ocean phosphate [mol P]
    pars.N0 = 4.35e16;     % Ocean nitrate [mol N]
    pars.O0 = 3.7e19;      % Atmospheric O2 [mol O]
    pars.Sr0 = 1.2e17;     % Ocean strontium [mol Sr]

    %% Additional: Extended nitrogen cycle's four main reservoirs (initial values set by user, as SCION not fully define them)
    % We define:
    % N_atm: Atmospheric N2
    % N_am: Ammonium (NH4+)
    % N_no3: Nitrate (NO3-)
    % N_org: Organic nitrogen
    N_atm0 = 1e20;   %% CUSTOM: need to be verify, just examplr to try!!!!!!!!!!!!
    N_am0  = 1e16;   %% CUSTOM
    N_no30 = 5e15;   %% CUSTOM
    N_org0 = 2e16;   %% CUSTOM

    % Isotopic standard ratios (carbon and nitrogen)
    pars.R_std = 0.011;         % 13C/12C standard ratio
    pars.R_N_std = 0.003676;    % 15N/14N standard ratio

    %% 2. Parameter settings (refer to Table 6 and Puyang original code)
    % Carbon cycle parameters
    pars.k_mocb = 2.5e12;
    pars.k_locb = 2.5e12;
    pars.k_ocdeg = 1.25e12;
    pars.k_ccdeg = 1.5e13;
    pars.k_carbw = 8e12;
    pars.k_sfw = 1.75e12;
    pars.k_mccb = 2.125e13;
    % %% CUSTOM: k_oxidw calculated as per SCION (e.g., k_oxidw = k_mocb + k_locb - k_ocdeg - k_reductant_input)
    pars.k_oxidw = 7.75e12;  %% need checked again %%

    % Sulfur cycle parameters
    pars.k_mpsb = 0.7e12;
    pars.k_mgsb = 1e12;
    pars.k_pyrw = 7e11;
    pars.k_gypw = 1e12;
    pars.k_pyrdeg = 2.5e11;
    pars.k_gypdeg = 5e11;

    % Phosphorus cycle parameters
    pars.k_capb = 2e10;
    pars.k_fepb = 1e10;

    % Nitrogen cycle parameters (extended part)
    pars.k_nfix = 8.67e12;
    pars.k_denit = 4.3e12;
    pars.k_nit = 1e-3;
    pars.K_nit = 0.1;
    pars.k_assim = 1e-3;
    pars.k_min = 1e-4;
    pars.K_denit = 0.1;

    % Sr cycle, etc.
    pars.k_Sr_sedw = 17e9;
    pars.k_Sr_mantle = 7.3e9;
    pars.k_Sr_silw = 13e9;
    pars.k_Sr_granw = pars.k_Sr_silw * 0.7;  %% CUSTOM: Adjusted based on basfrac
    pars.k_Sr_basw = pars.k_Sr_silw * 0.3;   %% CUSTOM
    pars.k_Sr_sfw = [];  %% CUSTOM: As defined in SCION
    pars.k_Sr_sedb = []; %% CUSTOM
    pars.k_Sr_metam = 13e9;

    % Other parameters
    pars.k_preplant = 0.15;
    pars.k_landfrac = 0.0588;
    pars.CP_biot = 250;
    pars.CP_lam = 1000;
    pars.CN_sea = 37.5;
    pars.atfrac0 = 0.01614;
    pars.k_mr = 3.762;
    pars.k_e = 2e-4;
    pars.chi_m = 0.1;
    pars.K = 6e-5;
    pars.kw = 1e-3;
    pars.Ea = 20;
    pars.z = 10;
    pars.sigma_plus1 = 0.9;
    pars.p_minimum = 10;
    pars.p_half = 183.6;
    pars.k_fire = 3;
    pars.k_scale = 200;

    %% 3. External inputs and perturbation parameters
    pars.F0 = 5e12;
    pars.A_long = 0.2e12;
    pars.A_short = 0.1e12;
    pars.delta_t = 0.5e6;

    %% 4. Temperature and feedback parameters
    pars.ESS = 5;
    pars.T0 = 15;

    %% 5. Initial state vector (extended to 11 dimensions)
    A0 = pars.A0; 
    G0 = pars.G0; 
    C0 = pars.C0; 
    N0_orig = pars.N0; 
    O0 = pars.O0;
    X_org0 = G0 * pars.R_std;
    X_N0 = N0_orig * pars.R_N_std;

    pars.startstate = [A0; G0; C0; N0_orig; X_org0; X_N0; O0; N_atm0; N_am0; N_no30; N_org0];

    %% 6. Model run time (yr, negative numbers represent the past)
    pars.whenstart = -240e6;
    pars.whenend = -118e6;

    %% 7. Load external geochemical data (strictly following SCION definition)
    % Carnian data (234–232 Ma)
    s1 = load('Geo_data_Carnian.mat');
    fields1 = fieldnames(s1);
    dataCarnian = s1.(fields1{1});
    [sortedAge, idx] = sort(dataCarnian.Age, 'ascend');
    dataCarnian.Age = sortedAge;
    fnames = fieldnames(dataCarnian);
    for i = 1:length(fnames)
        if ~strcmp(fnames{i}, 'Age')
            dataCarnian.(fnames{i}) = dataCarnian.(fnames{i})(idx);
        end
    end

    % J2 data (201–118 Ma)
    s2 = load('Geo_data_J2.mat');
    fields2 = fieldnames(s2);
    dataJ2 = s2.(fields2{1});
    [sortedAgeJ2, idxJ2] = sort(dataJ2.Age, 'ascend');
    dataJ2.Age = sortedAgeJ2;
    fnamesJ2 = fieldnames(dataJ2);
    for i = 1:length(fnamesJ2)
        if ~strcmp(fnamesJ2{i}, 'Age')
            dataJ2.(fnamesJ2{i}) = dataJ2.(fnamesJ2{i})(idxJ2);
        end
    end

    % J1 data (for productivity and redox, 201–118 Ma)
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

    % Construct model time vector (yr)
    model_time = linspace(-240e6, -118e6, 5000);

    % Stage 1: Carnian data (234–232 Ma)
    ind_stage1 = find(model_time >= -234e6 & model_time <= -232e6);
    seg1_time = model_time(ind_stage1)';
    [uniqueAge, ia] = unique(dataCarnian.Age, 'sorted');
    uniqueAge = uniqueAge(:);
    unique_13C = dataCarnian.Carbon_isotope(ia);
    unique_15N = dataCarnian.Nitrogen_isotope(ia);
    unique_13C = unique_13C(:);
    unique_15N = unique_15N(:);
    valid = isfinite(uniqueAge) & isfinite(unique_13C) & isfinite(unique_15N);
    uniqueAge = uniqueAge(valid);
    unique_13C = unique_13C(valid);
    unique_15N = unique_15N(valid);
    if length(uniqueAge) < 2
        data_seg1.delta13C = repmat(unique_13C, size(seg1_time));
        data_seg1.delta15N = repmat(unique_15N, size(seg1_time));
    else
        data_seg1.delta13C = interp1(uniqueAge, unique_13C, seg1_time, 'linear', 'extrap');
        data_seg1.delta15N = interp1(uniqueAge, unique_15N, seg1_time, 'linear', 'extrap');
    end

    % Stage 2: J2 data (201–118 Ma)
    ind_stage2 = find(model_time >= -201e6 & model_time <= -118e6);
    seg2_time = model_time(ind_stage2)';
    validJ2 = isfinite(dataJ2.Age) & isfinite(dataJ2.C_org);
    if sum(validJ2) < 2
        data_seg2.delta13C = repmat(dataJ2.C_org(find(validJ2,1)), size(seg2_time));
    else
        data_seg2.delta13C = interp1(dataJ2.Age(validJ2), dataJ2.C_org(validJ2), seg2_time, 'linear', 'extrap');
    end
    data_seg2.delta15N = nan(length(seg2_time), 1);

    forcings_data.model_time = model_time;
    forcings_data.delta13C = nan(size(model_time));
    forcings_data.delta15N = nan(size(model_time));
    forcings_data.delta13C(ind_stage1) = data_seg1.delta13C;
    forcings_data.delta15N(ind_stage1) = data_seg1.delta15N;
    forcings_data.delta13C(ind_stage2) = data_seg2.delta13C;
    forcings_data.delta15N(ind_stage2) = data_seg2.delta15N;
    
    pars.forcings_data = forcings_data;

    %% 8. Define external forcing function
    pars.forcing.F_in = @(t) external_forcing(t, pars);

    %% 9. Initialize global variables
    global stepnumber workingstate
    stepnumber = 1;
    workingstate = [];
    
    %% Return structure
    run.pars = pars;
    run.startstate = pars.startstate;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% external_forcing.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F = external_forcing(t, pars)
    t_abs = -t;  % yr, positive numbers represent the past
    F0 = pars.F0;
    delta_t = pars.delta_t;
    A_long = pars.A_long;
    A_short = pars.A_short;
    
    if t_abs >= 232e6 && t_abs <= 234e6
        F_long = A_long * sin(2*pi*(t_abs-232e6)/0.405e6);
        if abs(t_abs-234e6) < delta_t
            F_volc = A_short;
        else
            F_volc = 0;
        end
    elseif t_abs >= 118e6 && t_abs <= 201e6
        F_long = A_long * sin(2*pi*(t_abs-118e6)/8e6);
        Fv = 0;
        volcanic_events = [201e6,200e6; 183.2e6,182.4e6; 167e6,168e6; 160e6,160e6; ...
                           155e6,155e6; 146e6,144e6; 141e6,141e6; 134.4e6,133.2e6; 130.6e6,129.4e6];
        for j = 1:size(volcanic_events,1)
            t_start = volcanic_events(j,1);
            t_end = volcanic_events(j,2);
            if t_abs >= min(t_start,t_end) && t_abs <= max(t_start,t_end)
                Fv = Fv + A_short;
            end
        end
        F_volc = Fv;
    else
        F_long = 0;
        F_volc = 0;
    end
    F = F0 + F_long + F_volc;
end