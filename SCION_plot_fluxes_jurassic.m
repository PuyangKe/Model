%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% SCION - Spatial Continuous Integration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Earth Evolution Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Coded by BJW Mills
%%%% b.mills@leeds.ac.uk
%%%%
%%%% plot model fluxes

%%%% output to screen
fprintf('running plotting script 2... \t')
tic
global state

%%%% make figure
figure
 
%%%% make subplot  (Carbonate Carbon)
subplot(7,2,1)
%%%% load jurassic data
J2_data = load('Geo_data_J2_updated.mat') ;
C_carb = J2_data.dataJ2.C_carb ;
Age = J2_data.dataJ2.Age ;
%%%% plot data
plot(Age, C_carb, 'x')
hold on
%%%% plot modelwo
plot(state.time_myr,state.d13c_A);
%%%% just jurassic
xlim([-215 -115])
%%%% Title
title('C_{carb}')


%%%% make subplot   (Organic Carbon)
subplot(7,2,2)
%%%% load jurassic data
J2_data = load('Geo_data_J2_updated.mat') ;
C_org = J2_data.dataJ2.C_org ;
TOC = J2_data.dataJ2.TOC ;
Age = J2_data.dataJ2.Age ;
%%%% plot data
plot(Age, C_org, 'x')
hold on
%%%% plot modelwo
plot(state.time_myr,state.delta_G);
% plot(state.time_myr,state.tempC,'k')
% plot(state.time_myr,state.SAT_equator,'r')
%%%% just jurassic
xlim([-215 -115])
%%%% Title
title('C_{org} vs Buried organic C')

%%%% make subplot   (Carbon burial)
subplot(7,2,5)
%%%% plot this model
plot(state.time_myr,state.mocb,'b');
plot(state.time_myr,state.locb,'g');
plot(state.time_myr,state.mccb,'k');
hold on
%%%% Legend
text(-590,5e12,'mocb','color','b');
text(-590,4e12,'locb','color','g');
text(-590,16e12,'mccb','color','k')
%%%% just jurassic
xlim([-215 -115])
%%%% Title
title('marine_C_{carb} Burial')



%%%% make subplot   (Palaeoproductivity)
subplot(7,2,3)
%%%% load jurassic data
J1_data = load('Geo_data_J1_updated.mat') ;
P = J1_data.forcings.P ;
Ba = J1_data.forcings.Ba ;
Age = J1_data.forcings.Age ;
%%%% plot data
% plot(Age, P, 'x')
plot(Age, Ba, 'x')
hold on
%%%% just jurassic
xlim([-215 -115])
%%%% Title
title('Palaeoproductivity')

%%%% make subplot 
subplot(7,2,6)
%%%% plot modelwo
plot(state.time_myr,state.psea);
hold on

%%%% just jurassic
xlim([-215 -115])
%%%% Title
title('P_flux')




%%%% make subplot   (Redox conditions)
subplot(7,2,4)
%%%% load jurassic data
J1_data = load('Geo_data_J1_updated.mat') ;
Mn = J1_data.forcings.Mn ;
Age = J1_data.forcings.Age ;
%%%% plot data
plot(Age, Mn, 'x')
hold on
%%%% plot modelwo
% plot(state.time_myr,state.denit)
%%%% just jurassic
xlim([-215 -115])
%%%% Title
title('Redox conditions')

%%%% make subplot   (Redox conditions)
subplot(7,2,7)
%%%% plot modelwo
plot(state.time_myr,state.denit)
%%%% just jurassic
xlim([-215 -115])
%%%% Title
title('denit')

%%%% CO2ppm;
subplot(7,2,8);
xlabel('Time (Ma)');
ylabel('Atmospheric CO_{2} (ppm)');
%%%% plot this model
plot(state.time_myr,state.RCO2.*280,'k')
%%%% just jurassic
xlim([-215 -115])
%%%% Title
title('CO2ppm')

%%%% O2 (%) 
subplot(7,2,9);
xlabel('Time (Ma)');
ylabel('Atmospheric O_{2} (%)');
%%%% plot this model
plot(state.time_myr,state.mrO2.*100,'k')
%%%% just jurassic
xlim([-215 -115])
%%%% Title
title('O2')


%%%% C input  
subplot(7,2,10);
xlabel('Time (Ma)');
ylabel('C input');
%%%% plot this model
plot(state.time_myr,state.C_volc,'k')
%%%% just jurassic
xlim([-215 -115])
%%%% Title
title('C input')

%%%%% plotting script finished
fprintf('Done: ')
endtime = toc ;
fprintf('time (s): %d \n', endtime )
% 
% %%%%% ============ Astronomical forcing and responses ============
% 
% % 如果存在天文扰动数据，则绘制
% if exist('T_ext_data', 'var') && isfield(T_ext_data, 'T_ext_data')
% 
%     % 插值 T_offset
%     T_offset = interp1(Time_Ma, T_ext_K, state.time_myr, 'linear', 0);
% 
% 
%     % 环境温度（原始 + 扰动）
%     T_env = state.tempC + T_offset;
% 
%     % === subplot 11: T_ext 与 T_offset ===
%     subplot(7,2,11)
%     plot(T_ext_data.T_ext_data.Time_Ma, T_ext_data.T_ext_data.T_ext_K, 'm', 'DisplayName', 'T_{ext}')
%     hold on
%     plot(state.time_myr, T_offset, 'r', 'DisplayName', 'T_{offset}')
%     xlabel('Time (Ma)'); ylabel('Temp (K)');
%     title('T_{ext} & T_{offset}');
%     legend show; xlim([-215 -115])
% 
%     % === subplot 12: 环境温度 T_env vs 模型原始 ===
%     subplot(7,2,12)
%     plot(state.time_myr, state.tempC, 'k--', 'DisplayName', 'T_{model}')
%     hold on
%     plot(state.time_myr, T_env, 'b', 'DisplayName', 'T_{env} = T + offset')
%     xlabel('Time (Ma)'); ylabel('Temp (°C)');
%     title('Global Temp: Model vs Env (T_{env})');
%     legend show; xlim([-215 -115])
% end
% 
% % === subplot 13: Silicate weathering (silw) ===
% figure(gcf); 
% subplot(7,2,13)
% plot(state.time_myr, state.silw, 'g', 'DisplayName', 'Silw')
% xlabel('Time (Ma)'); ylabel('Silicate weathering');
% title('Silicate Weathering Flux (silw)');
% xlim([-215 -115]); legend show
% 
% % === subplot 14: Organic C burial (mocb) ===
% subplot(7,2,14)
% plot(state.time_myr, state.mocb, 'c', 'DisplayName', 'mocb')
% xlabel('Time (Ma)'); ylabel('Org. C burial');
% title('Organic Carbon Burial (mocb)');
% xlim([-215 -115]); legend show
% 
% fprintf('Plotting astronomical forcing terms...\n');
% 
% % === 加载三个驱动项（8Ma / 4.8Ma / 组合）===
% load('8Ma_phase_signal.mat');     % 包含：time8_Ma, signal8
% load('4p8Ma_phase_signal.mat');   % 包含：time48_Ma, signal48
% load('T_ext_8Ma_4p8Ma.mat');      % 包含：Time_Ma, T_ext_K
% 
% % === 构建插值后驱动项（按叠加公式生成）===
% alpha = 3.0;
% beta = 1.5;
% signal8_interp = interp1(time8_Ma, signal8, Time_Ma, 'linear', 0);
% signal48_interp = interp1(time48_Ma, signal48, Time_Ma, 'linear', 0);
% 
% % === 绘图 ===
% figure;
% subplot(3,1,1)
% plot(time8_Ma, signal8, 'b')
% xlabel('Time (Ma ago)'); ylabel('Normalized Phase');
% title('8 Ma cycle (normalized)');
% xlim([-215 -115])
% 
% subplot(3,1,2)
% plot(time48_Ma, signal48, 'g')
% xlabel('Time (Ma ago)'); ylabel('Normalized Phase');
% title('4.8 Ma cycle (normalized)');
% xlim([-215 -115])
% 
% subplot(3,1,3)
% plot(Time_Ma, T_ext_K, 'r')
% xlabel('Time (Ma ago)'); ylabel('T_{ext} (K)');
% title('Astronomical forcing term (combined 8 + 4.8 Ma)');
% xlim([-215 -115])
% 
% fprintf('T_{ext} curves plotted.\n');


% load('state_S1.mat'); SCION_plot_fluxes_full(state, T_ext_data, 'S1');
% load('state_S2.mat'); SCION_plot_fluxes_full(state, T_ext_data, 'S2');
% load('state_S3.mat'); SCION_plot_fluxes_full(state, T_ext_data, 'S3');
% load('state_S4.mat'); SCION_plot_fluxes_full(state, T_ext_data, 'S4');

% load('state_S1.mat'); SCION_plot_fluxes_full(state, [], 'S1');
% load('state_S2.mat'); T_ext_8Ma = load('8Ma_phase_signal.mat');
% SCION_plot_fluxes_full(state, T_ext_8Ma, 'S2');
% 
% load('state_S3.mat'); T_ext_4p8Ma = load('4p8Ma_phase_signal.mat');
% SCION_plot_fluxes_full(state, T_ext_4p8Ma, 'S3');
% 
% load('state_S4.mat'); T_ext_combo = load('T_ext_8Ma_4p8Ma.mat');
% SCION_plot_fluxes_full(state, T_ext_combo, 'S4');
% 
% function SCION_plot_fluxes_full(state, T_ext_data, scenarioname)
% 
% fprintf('Running plotting for scenario %s...\n', scenarioname);
% 
% figure('Position',[100,100,1200,900])
% 
% subplot(7,2,1)
% J2_data = load('Geo_data_J2_updated.mat');
% plot(J2_data.dataJ2.Age, J2_data.dataJ2.C_carb, 'xk'); hold on
% plot(state.time_myr, state.d13c_A, 'b');
% title('\delta^{13}C (carb)'); xlim([-215 -115])
% 
% subplot(7,2,2)
% plot(J2_data.dataJ2.Age, J2_data.dataJ2.C_org, 'xk'); hold on
% plot(state.time_myr, state.delta_G, 'g');
% title('C_{org} vs \deltaG'); xlim([-215 -115])
% 
% subplot(7,2,3)
% J1_data = load('Geo_data_J1_updated.mat');
% plot(J1_data.forcings.Age, J1_data.forcings.Ba, 'xk');
% title('Ba proxy for productivity'); xlim([-215 -115])
% 
% subplot(7,2,4)
% plot(J1_data.forcings.Age, J1_data.forcings.Mn, 'xk'); hold on
% plot(state.time_myr, 1 - state.ANOX, 'm');
% title('Mn vs 1 - ANOX'); xlim([-215 -115])
% 
% subplot(7,2,5)
% plot(state.time_myr, state.mocb, 'b'); hold on
% plot(state.time_myr, state.locb, 'g');
% plot(state.time_myr, state.mccb, 'k');
% title('C Burial: mocb/locb/mccb'); xlim([-215 -115])
% 
% subplot(7,2,6)
% plot(state.time_myr, state.psea, 'c');
% title('P flux (psea)'); xlim([-215 -115])
% 
% subplot(7,2,7)
% plot(state.time_myr, state.denit, 'r');
% title('Denitrification'); xlim([-215 -115])
% 
% subplot(7,2,8)
% plot(state.time_myr, state.RCO2 * 280, 'k');
% title('CO_2 (ppm)'); xlim([-215 -115])
% 
% subplot(7,2,9)
% plot(state.time_myr, state.mrO2 * 100, 'k');
% title('O_2 (%)'); xlim([-215 -115])
% 
% subplot(7,2,10)
% plot(state.time_myr, state.C_volc, 'k');
% title('C_{volc} input'); xlim([-215 -115])
% 
% subplot(7,2,11)
% if isfield(T_ext_data, 'T_ext_K')
%     T_offset = interp1(T_ext_data.Time_Ma, T_ext_data.T_ext_K, state.time_myr, 'linear', 0);
%     plot(T_ext_data.Time_Ma, T_ext_data.T_ext_K, 'm'); hold on
%     plot(state.time_myr, T_offset, 'r');
%     title('T_{ext} & T_{offset}'); xlim([-215 -115])
% else
%     title('T_{ext} missing');
% end
% 
% subplot(7,2,12)
% if exist('T_offset', 'var')
%     T_env = state.tempC + T_offset;
%     plot(state.time_myr, state.tempC, 'k--'); hold on
%     plot(state.time_myr, T_env, 'b');
%     title('T_{env} = T + T_{offset}'); xlim([-215 -115])
% else
%     title('T_{env} unavailable');
% end
% 
% subplot(7,2,13)
% plot(state.time_myr, state.silw, 'g');
% title('Silicate weathering (silw)'); xlim([-215 -115])
% 
% subplot(7,2,14)
% plot(state.time_myr, state.mocb, 'c');
% title('Organic Carbon Burial (mocb)'); xlim([-215 -115])
% 
% saveas(gcf, ['SCION_fluxes_' scenarioname '.png']);
% 
% fprintf('Plot for scenario %s saved.\n', scenarioname);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Volcanic Forcing
% % Time (Main phase)
% time_nodes = [-201e6:1:-200e6, -190e6, -183e6, -168e6, -160e6, -155e6, -146e6, -141e6, -134.4e6, -130.6e6];
% % release ratios of CO₂ (mol/yr) Main phase
% volcanic_flux_nodes = [1.18e11, 7.79e11, 7.79e11, 0, 0, 0, 0, 0, 0, 0];
% 
% % Time (Total Duration of LIP)
% time_nodes = [-205e6:1:-191e6, -190e6, -180e6, -170e6, -160e6, -150e6, -140e6, -130e6, -120e6, -115e6];
% % release ratios of CO₂ (mol/yr) Total Duration of LIP
% volcanic_flux_nodes = [8.42e9, 4.15e10, 4.15e10, 0, 0, 0, 0, 0, 0, 0];


% %%%%%%% Volcanic CO2 injection (mol/yr)
% C_volc = 0;  % Default: no volcanic CO2
% 
% % CAMP main volcanic phase: 201–200 Ma, with estimated flux ~1.18e11 mol/yr
% if t_geol >= 200 && t_geol <= 201
%     C_volc = 1.18e11;
% end
% 
% % Karoo-Ferrar main volcanic phase: 183.2–182.4 Ma, with estimated flux ~7.79e11 mol/yr
% if t_geol >= 182.4 && t_geol <= 183.2
%     C_volc = 7.79e11;
% end
% 
% 
% %%% Carbon dioxide
% dy(3) = -locb - mocb + oxidw + ocdeg + ccdeg + carbw - mccb - sfw  + reductant_input + C_volc;
% 
% %%% delta_A * A
% dy(16) = -locb*(  delta_locb ) -mocb*( delta_mocb ) + oxidw*delta_G + ocdeg*delta_G + ccdeg*delta_C + carbw*delta_C - mccb*delta_mccb - sfw*delta_mccb + reductant_input*d13C_volc + C_volc*d13C_volc;
