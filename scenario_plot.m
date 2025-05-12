% ===============================================
% SCION Jurassic plot（S1–S7）
% include δ13C、CO₂、T、psea、Mn、Volcanic CO₂ input
% ===============================================

fprintf('Plotting SCION scenario comparison with geochemical proxies...\n');
figure;

% four different scenarios
load('state_S1.mat'); state_S1 = run.state;
load('state_S2.mat'); state_S2 = run.state;
load('state_S3.mat'); state_S3 = run.state;
load('state_S4.mat'); state_S4 = run.state;
load('state_S5.mat'); state_S5 = run.state;
load('state_S6.mat'); state_S6 = run.state;
load('state_S7.mat'); state_S7 = run.state;

% %%%% Load forcing data
% load('T_ext_8Ma_4p8Ma.mat'); 
% load('8Ma_phase_signal.mat');
% load('4p8Ma_phase_signal.mat');

%%%% Interpolate T_offset for each scenario
T_offset_S1 = zeros(size(state_S1.time_myr));
% T_offset_S2 = interp1(Time_Ma, T_ext_K, state_S2.time_myr, 'linear', 0);
% T_offset_S3 = interp1(Time_Ma, T_ext_K, state_S3.time_myr, 'linear', 0);
% T_offset_S4 = interp1(Time_Ma, T_ext_K, state_S4.time_myr, 'linear', 0);
% % S2: 8Ma
% load('8Ma_phase_signal.mat'); 
% T_offset_S2 = interp1(time8_Ma, signal8, state_S2.time_myr, 'linear', 0);
% 
% % S3: 4.8Ma
% load('4p8Ma_phase_signal.mat');
% T_offset_S3 = interp1(time48_Ma, signal48, state_S3.time_myr, 'linear', 0);
% 
% % S4: 8 + 4.8Ma
% load('T_ext_8Ma_4p8Ma.mat');
% T_offset_S4 = interp1(Time_Ma, T_ext_K, state_S4.time_myr, 'linear', 0);
load('8Ma_phase_signal.mat');
T_offset_S2 = interp1(time8_Ma, signal8, state_S2.time_myr, 'linear', 0);
T_offset_S6 = interp1(time8_Ma, signal8, state_S6.time_myr, 'linear', 0);

load('4p8Ma_phase_signal.mat');
T_offset_S3 = interp1(time48_Ma, signal48, state_S3.time_myr, 'linear', 0);
T_offset_S5 = interp1(time48_Ma, signal48, state_S5.time_myr, 'linear', 0);

load('T_ext_8Ma_4p8Ma.mat');
T_offset_S4 = interp1(Time_Ma, T_ext_K, state_S4.time_myr, 'linear', 0);
T_offset_S7 = interp1(Time_Ma, T_ext_K, state_S7.time_myr, 'linear', 0);

%%%% Plotting
figure('Position', [100, 100, 1200, 1000]);

% --- Subplot 1: δ13C ---
subplot(7,2,1)
% subplot(1,1,1)
J2 = load('Geo_data_J2_updated.mat');
plot(J2.dataJ2.Age, J2.dataJ2.C_carb, 'or', 'DisplayName', 'Obs δ^{13}C');  
plot(J2.dataJ2.Age, J2.dataJ2.C_carb, 'o', ...
    'MarkerSize', 3, ...
    'MarkerEdgeColor', [1, 0.4, 0.4], ...
    'MarkerFaceColor', [1, 0.4, 0.4], ...
    'DisplayName', 'Obs \delta^{13}C'); 
hold on; 
plot(state_S1.time_myr, state_S1.d13c_A, 'k');
plot(state_S2.time_myr, state_S2.d13c_A, 'b');
plot(state_S3.time_myr, state_S3.d13c_A, 'g');
plot(state_S4.time_myr, state_S4.d13c_A, 'r');
plot(state_S5.time_myr, state_S5.d13c_A, 'c');
plot(state_S6.time_myr, state_S6.d13c_A, 'm');
plot(state_S7.time_myr, state_S7.d13c_A, '--k');
xlim([-205 -120]); ylabel('\delta^{13}C'); title('\delta^{13}C');
legend('Location','best');

% --- Subplot 2: CO₂ ---
subplot(7,2,2)
plot(state_S1.time_myr, state_S1.RCO2*280, 'k'); hold on;
plot(state_S2.time_myr, state_S2.RCO2*280, 'b');
plot(state_S3.time_myr, state_S3.RCO2*280, 'g');
plot(state_S4.time_myr, state_S4.RCO2*280, 'r');
plot(state_S5.time_myr, state_S5.RCO2*280, 'c');
plot(state_S6.time_myr, state_S6.RCO2*280, 'm');
plot(state_S7.time_myr, state_S7.RCO2*280, 'Color', [0.5 0.5 0.5]);
xlim([-205 -120]); ylabel('ppm'); title('CO₂');

% --- Subplot 3: GAST (tempC) ---
subplot(7,2,3)
% subplot(1,1,1)
plot(state_S1.time_myr, state_S1.tempC, 'k'); 
hold on;
plot(state_S2.time_myr, state_S2.tempC, 'b');
plot(state_S3.time_myr, state_S3.tempC, 'g');
plot(state_S4.time_myr, state_S4.tempC, 'r');
plot(state_S5.time_myr, state_S5.tempC, 'c');
plot(state_S6.time_myr, state_S6.tempC, 'm');
plot(state_S7.time_myr, state_S7.tempC, 'Color', [0.5 0.5 0.5]);
xlim([-205 -120]); ylabel('°C'); title('Global Mean Temp');

% --- Subplot 4: psea ---
subplot(7,2,4)
% subplot(1,1,1)
% J1 = load('Geo_data_J1_updated.mat');
% plot(J1.forcings.Age, J1.forcings.P, 'o'); hold on;
% plot(J1.forcings.Age, J1.forcings.P, 'o', ...
%     'MarkerSize', 3, ...
%     'MarkerEdgeColor', [1, 0.4, 0.4], ...
%     'MarkerFaceColor', [1, 0.4, 0.4], ...
%     'DisplayName', 'P'); 
hold on;
plot(state_S1.time_myr, state_S1.psea, 'k'); 
plot(state_S2.time_myr, state_S2.psea, 'b');
plot(state_S3.time_myr, state_S3.psea, 'g');
plot(state_S4.time_myr, state_S4.psea, 'r');
plot(state_S5.time_myr, state_S5.psea, 'c');
plot(state_S6.time_myr, state_S6.psea, 'm');
plot(state_S7.time_myr, state_S7.psea, 'Color', [0.5 0.5 0.5]);
xlim([-205 -120]); ylabel('P flux'); title('Marine P input (psea)');

% --- Subplot 5: Mn proxy vs 1-ANOX ---
subplot(7,2,5)
J1 = load('Geo_data_J1_updated.mat');
plot(J1.forcings.Age, J1.forcings.Mn, 'xk'); hold on;
plot(J1.forcings.Age, J1.forcings.Mn, 'o', ...
    'MarkerSize', 3, ...
    'MarkerEdgeColor', [0.6, 0.3, 0.4], ...
    'MarkerFaceColor', [0.6, 0.3, 0.4], ...
    'DisplayName', 'Mn'); 
hold on;
plot(state_S1.time_myr, 1 - state_S1.ANOX, 'k');
plot(state_S2.time_myr, 1 - state_S2.ANOX, 'b');
plot(state_S3.time_myr, 1 - state_S3.ANOX, 'g');
plot(state_S4.time_myr, 1 - state_S4.ANOX, 'r');
plot(state_S5.time_myr, state_S5.ANOX, 'c');
plot(state_S6.time_myr, state_S6.ANOX, 'm');
plot(state_S7.time_myr, state_S7.ANOX, 'Color', [0.5 0.5 0.5]);
xlim([-205 -120]); ylabel('1 - ANOX'); title('Redox proxy (Mn)');

% --- Subplot 6: Volcanic CO₂ input ---
subplot(7,2,6)
plot(state_S1.time_myr, state_S1.C_volc, 'k');
xlim([-205 -120]); ylabel('mol/yr'); title('Volcanic CO₂ Input');

% --- Subplot 7: Astronomical T_ext forcing ---
subplot(7,2,7)
plot(state_S2.time_myr, T_offset_S2, 'b', 'DisplayName', 'S2 8Ma'); hold on;
plot(state_S3.time_myr, T_offset_S3, 'g', 'DisplayName', 'S3 4.8Ma');
plot(state_S4.time_myr, T_offset_S4, 'r', 'DisplayName', 'S4 8+4.8Ma');
xlim([-205 -120]); ylabel('\DeltaT (K)'); title('Astronomical T_{ext}');
legend('Location','best');

% --- Subplot 8: Astronomical effect on mocb ---
subplot(7,2,8)
plot(state_S1.time_myr, state_S1.mocb, 'k'); hold on;
plot(state_S2.time_myr, state_S2.mocb, 'b');
plot(state_S3.time_myr, state_S3.mocb, 'g');
plot(state_S4.time_myr, state_S4.mocb, 'r');
plot(state_S5.time_myr, state_S5.mocb, 'c');
plot(state_S6.time_myr, state_S6.mocb, 'm');
plot(state_S7.time_myr, state_S7.mocb, 'Color', [0.5 0.5 0.5]);
xlim([-205 -120]); ylabel('mol/yr'); title('Marine OC burial (mocb)');

% --- Subplot 9: silicate weathering (silw) ---
subplot(7,2,9)
plot(state_S1.time_myr, state_S1.silw, 'k'); hold on;
plot(state_S2.time_myr, state_S2.silw, 'b');
plot(state_S3.time_myr, state_S3.silw, 'g');
plot(state_S4.time_myr, state_S4.silw, 'r');
plot(state_S5.time_myr, state_S5.silw, 'c');
plot(state_S6.time_myr, state_S6.silw, 'm');
plot(state_S7.time_myr, state_S7.silw, 'Color', [0.5 0.5 0.5]);
xlim([-205 -120]); ylabel('mol/yr'); title('Silicate weathering');

% --- Subplot 10: NPP/VEG ---
subplot(7,2,10)
plot(state_S1.time_myr, state_S1.VEG, 'k'); hold on;
plot(state_S2.time_myr, state_S2.VEG, 'b');
plot(state_S3.time_myr, state_S3.VEG, 'g');
plot(state_S4.time_myr, state_S4.VEG, 'r');
plot(state_S5.time_myr, state_S5.VEG, 'c');
plot(state_S6.time_myr, state_S6.VEG, 'm');
plot(state_S7.time_myr, state_S7.VEG, 'Color', [0.5 0.5 0.5]);
xlim([-205 -120]); ylabel('arb'); title('VEG (NPP proxy)');

fprintf(' Scenario plotting complete.\n');





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear; clc;
% 
% load('state_S2.mat'); state_S2 = run.state;
% load('state_S3.mat'); state_S3 = run.state;
% load('state_S4.mat'); state_S4 = run.state;
% load('state_S5.mat'); state_S5 = run.state;
% load('state_S6.mat'); state_S6 = run.state;
% load('state_S7.mat'); state_S7 = run.state;
% 
% %% 2. 载入天文信号与观测 δ13C 数据
% d8   = load('8Ma_phase_signal.mat');        % 8Ma 单周期
% d4   = load('4p8Ma_phase_signal.mat');      % 4.8Ma 单周期
% d84  = load('T_ext_8Ma_4p8Ma.mat');         % 合成 8+4.8Ma
% J2   = load('Geo_data_J2_updated.mat');     % 观测 δ13C
% obs_age  = J2.dataJ2.Age;
% obs_d13C = J2.dataJ2.C_carb;
% 
% %% 3. 准备场景列表
% scenarios = {'S2','S3','S4','S5','S6','S7'};
% 
% %% 4. 循环绘图并做交叉相关
% figure('Position',[100 100 800 1800]);
% for i = 1:numel(scenarios)
%     scen = scenarios{i};
%     st   = eval(['state_' scen]);    % 对应的 state 结构体
%     t    = st.time_myr;               % 时间向量 (Ma, 负值表示过去)
% 
%     % 4.1 计算 T_offset
%     switch scen
%       case 'S2'   % only 8Ma
%         sig     = interp1(d8.time8_Ma, d8.signal8, t, 'linear', 0);
%         T_off   = (sig + 1) * 0.3;
%       case 'S3'   % 4.8Ma + volcanism
%         sig     = interp1(d4.time48_Ma, d4.signal48, t, 'linear', 0);
%         T_off   = sig;                % 与模型中一致，不缩放
%       case 'S4'   % 8Ma + 4.8Ma + volcanism
%         T_off   = interp1(d84.Time_Ma, d84.T_ext_K, t, 'linear', 0);
%       case 'S5'   % only 4.8Ma
%         sig     = interp1(d4.time48_Ma, d4.signal48, t, 'linear', 0);
%         T_off   = (sig + 1) * 0.3;
%       case 'S6'   % 8Ma + volcanism
%         sig     = interp1(d8.time8_Ma, d8.signal8, t, 'linear', 0);
%         T_off   = (sig + 1) * 0.3;
%       case 'S7'   % 8Ma + 4.8Ma (no volcanism)
%         T_off   = interp1(d84.Time_Ma, d84.T_ext_K, t, 'linear', 0);
%     end
% 
%     % 4.2 绘制
%     subplot(numel(scenarios),1,i);
%     hold on;
%     plot(t,        T_off,    '-b', 'LineWidth',1.5, 'DisplayName','T\_offset');
%     plot(t,        st.d13c_A,'-r', 'LineWidth',1.5, 'DisplayName','Model \delta^{13}C');
%     plot(obs_age,  obs_d13C, '.k','MarkerSize',8,      'DisplayName','Obs \delta^{13}C');
%     xlim([min(t) max(t)]);          % 负值表示过去，横轴无需翻转
%     ylabel('\DeltaT (K) & \delta^{13}C (‰)');
%     title(['Scenario ' scen]);
%     if i == numel(scenarios)
%         xlabel('Time (Ma, negative = past)');
%     end
%     legend('Location','best');
%     grid on;
% 
%     % 4.3 交叉相关分析：Model δ13C vs. T_offset
%     [cc, lags] = xcorr(T_off - nanmean(T_off), ...
%                       st.d13c_A - nanmean(st.d13c_A), ...
%                       50, 'coeff');
%     [~, idx] = max(abs(cc));
%     bestLag  = lags(idx);
% 
%     if cc(idx) > 0
%         corrStr = '正相关';
%     else
%         corrStr = '负相关';
%     end
% 
%     fprintf('Scenario %s: cross-corr = %.3f (%s), best lag = %d steps\n', ...
%             scen, cc(idx), corrStr, bestLag);
% end
