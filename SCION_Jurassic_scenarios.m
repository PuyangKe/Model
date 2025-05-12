% Run all SCION Jurassic scenarios (S1–S7)
% S1: Volcanism only
% S2: 8Ma cycle only
% S3: 4.8Ma + volcanism
% S4: 8Ma + 4.8Ma + volcanism
% S5: 4.8Ma cycle only
% S6: 8Ma + volcanism
% S7: 8Ma + 4.8Ma (no volcanism)

fprintf('\nSCION S1: Volcanism only\n');
global T_ext_data scenarioID
scenarioID = 'S1';
T_ext_data = load('T_ext_none.mat');
run = SCION_initialise(-1);
save('state_S1.mat', 'run');

fprintf('\nSCION S2: 8 Ma cycle only\n');
scenarioID = 'S2';
T_ext_data = load('8Ma_phase_signal.mat');
run = SCION_initialise(-1);
save('state_S2.mat', 'run');

fprintf('\nSCION S3: 4.8 Ma + volcanism\n');
scenarioID = 'S3';
T_ext_data = load('4p8Ma_phase_signal.mat');
run = SCION_initialise(-1);
save('state_S3.mat', 'run');

fprintf('\nSCION S4: 8 + 4.8 Ma + volcanism\n');
scenarioID = 'S4';
T_ext_data = load('T_ext_8Ma_4p8Ma.mat');
run = SCION_initialise(-1);
save('state_S4.mat', 'run');

fprintf('\nSCION S5: 4.8 Ma cycle only\n');
scenarioID = 'S5';
T_ext_data = load('4p8Ma_phase_signal.mat');
run = SCION_initialise(-1);
save('state_S5.mat', 'run');

fprintf('\nSCION S6: 8 Ma + volcanism\n');
scenarioID = 'S6';
T_ext_data = load('8Ma_phase_signal.mat');
run = SCION_initialise(-1);
save('state_S6.mat', 'run');

fprintf('\nSCION S7: 8 + 4.8 Ma only (no volcanism)\n');
scenarioID = 'S7';
T_ext_data = load('T_ext_8Ma_4p8Ma.mat');
run = SCION_initialise(-1);
save('state_S7.mat', 'run');

fprintf('\n✅ All scenarios (S1–S7) completed and saved.\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



