% Коэффициенты Тейлора из Приложения А публикации E.Y. Kutyakov, S.V. Dushin,
% A.B. Iskakov, A.N. Abramenkov, Taylor series decomposition of nonlinear 
% model of two-area four-generators power system
%% 1. Загружаем исходные данные

clear all
init_data = ID_4gen11bus();

%% 2. Вычисляем матрицу проводимостей для сети из 6-ти узлов
admit_mat = fadmit_mat(init_data.line);

%% 3. Вычисляем установившиеся значения напряжений, мощностей и токов в узлах
power_flow = fpower_flow(init_data.power_flow,admit_mat.Y_cmplx_l);

%% 4. Вычисляем матрицы коэффициентов Тейлора

% 4.1.  Матрицы генераторов A_G(i),A_G(i)_2,B_G(i),B_G(i)_2,C_G(i),C_G(i)_2,
%                           D_G(i),H_G(i),S_G(i) по формулам из Приложения А
bilin_gen = fbilin_gen (init_data.machine,power_flow);

% 4.2. Матрицы нагрузок D_L(i) по формулам из Приложения А
lin_load = flin_load (init_data.machine.n_gen,...
    init_data.power_flow.n_load,power_flow);
 
%% 5. Компоновка диагональных матриц системы (Ad, Bd, Cd, Dd)
sys_mat = fsys_mat(init_data.machine.n_gen, init_data.power_flow.n_load,...
    bilin_gen,lin_load);
%%                             ФУНКЦИИ

%============================INITIAL DATA==================================
function [init_data] = ID_4gen11bus()

% Параметры генераторов в системе машин 900MVA/20kV (ex.12.6 (p.813))
n_gen = 4; % количество генераторов
X_d(1:n_gen,1) = 1.80; X_d_ht(1:n_gen,1) = 0.3; X_d_2ht(1:n_gen,1) = 0.25;
X_q(1:n_gen,1) = 1.70; X_q_ht(1:n_gen,1) = 0.55; X_q_2ht(1:n_gen,1) = 0.25;

T_d0_ht(1:n_gen,1) = 8;    T_d0_2ht(1:n_gen,1) = 0.03;
T_q0_ht(1:n_gen,1) = 0.4;  T_q0_2ht(1:n_gen,1) = 0.05;

R_a(1:n_gen,1) = 0.0025;    X_l(1:n_gen,1) = 0.2;

A_Sat(1:n_gen,1) = 0.015; B_Sat(1:n_gen,1) = 9.6; psi_TI(1:n_gen,1) = 0.9;

H = [6.5;6.5;6.175;6.175]; %[G1;G2;G3;G4]

K_D(1:n_gen, 1) = 0;

f = 60; % частота сети в Гц
omega_0(1:n_gen, 1) = 2*pi*f; % частота сети в радианах

init_data.machine = struct('omega_0',omega_0, 'X_d',X_d, 'X_d_ht',X_d_ht,...
    'X_d_2ht',X_d_2ht, 'X_q',X_q, 'X_q_ht',X_q_ht, 'X_q_2ht',X_q_2ht,...
    'T_d0_ht',T_d0_ht,'T_d0_2ht',T_d0_2ht, 'T_q0_ht',T_q0_ht,...
    'T_q0_2ht',T_q0_2ht, 'R_a',R_a, 'X_l',X_l, 'A_Sat',A_Sat,...
    'B_Sat',B_Sat, 'psi_TI',psi_TI, 'H',H, 'K_D',K_D, 'n_gen', n_gen);

% Мощности и напряжения генераторных шин:
P = [700;700;719;700]; % активная мощность в MW, [G1;G2;G3;G4]
Q = [185;235;176;202]; % реактивная мощность в MVAr, [G1;G2;G3;G4]
absE = [1.03;1.01;1.03;1.01]; % модуль комплексного напряжения в pu системы
                              % линий 100MVA/230kV
angE = [20.2;10.5;-6.8;-17.0];% угол наколона вектора напряжений в градусах

% Мощности нагрузок:
P_load = [967;1767]; % активная мощность в MW, [bus7;bus9]
Q_load = [200-100;350-100]; % реактивная мощность в MVAr, [bus7;bus9]
n_load = 2; % количество шин с нагрузами

init_data.power_flow = struct('P',P, 'Q',Q, 'absE',absE, 'angE',angE,...
    'P_load',P_load, 'Q_load',Q_load, 'n_load',n_load);

% Параметры линий в системе линий:
r = 0.0001; % удельное активное сопротивление [p.u. / km]
x_L = 0.001; % удельное реактивное сопротивление [p.u / km]
b_C = 0.00175i; % удельная проводимость параллельного емкостного
                % сопротивления [p.u / km]

r_transf_m = 0.15i; % реактивное сопротивление трансформатора [p.u.]
                    % в системе машин 900MVA/20kV

init_data.line = struct('r',r, 'x_L',x_L, 'b_C',b_C, 'r_transf_m',r_transf_m);
end
%========================END of INITIAL DATA===============================


%=========================ADMITTANCE MATRIX================================
function [admit_mat] = fadmit_mat(data_line)

r = data_line.r;
x_L = data_line.x_L;
b_C = data_line.b_C;
r_transf_m = data_line.r_transf_m;

% Полное удельное сопротивление линии [p.u. / km] в системе линий 100MVA/230 kV
R_line = r + x_L * 1i;

% Перевод сопротивления трансформатора из системы машин 900MVA/230kV в
% систему линий 100MVA/230kV:

r_transf_line = r_transf_m * 100 / 900;

% Вычисление проводимостей элементов 11-ти узловой сети:

Y_11bus.Y11 = 1 / r_transf_line;
Y_11bus.Y15 = - Y_11bus.Y11;
Y_11bus.Y22 = Y_11bus.Y11;
Y_11bus.Y26 = Y_11bus.Y15;
Y_11bus.Y33 = Y_11bus.Y11;
Y_11bus.Y3_11 = Y_11bus.Y15;
Y_11bus.Y44 = Y_11bus.Y11;
Y_11bus.Y4_10 = Y_11bus.Y15;
Y_11bus.Y51 = Y_11bus.Y15;
Y_11bus.Y55 = 1 / r_transf_line + ...
    1 / (25 * R_line);
Y_11bus.Y56 = - 1 / (25 * R_line);
Y_11bus.Y62 = Y_11bus.Y26;
Y_11bus.Y65 = Y_11bus.Y56;
Y_11bus.Y66 = 1 / r_transf_line + ...
    1 / (25 * R_line) + 1 / (10 * R_line);
Y_11bus.Y67 = - 1 / (10 * R_line);
Y_11bus.Y76 = Y_11bus.Y67;
Y_11bus.Y77 = 1 / (10 * R_line) + ...
    2 * (1 / (110 * R_line) + b_C * 110 / 2);
Y_11bus.Y78 = - 2 / (110 * R_line);
Y_11bus.Y87 = Y_11bus.Y78;
Y_11bus.Y88 = 4 * (1 / (110 * R_line) ...
    + b_C * 110 / 2);
Y_11bus.Y89 = Y_11bus.Y78;
Y_11bus.Y98 = Y_11bus.Y89;
Y_11bus.Y99 = Y_11bus.Y77;
Y_11bus.Y9_10 = Y_11bus.Y76;
Y_11bus.Y10_4 = Y_11bus.Y4_10;
Y_11bus.Y10_9 = Y_11bus.Y9_10;
Y_11bus.Y10_10 = Y_11bus.Y66;
Y_11bus.Y10_11 = Y_11bus.Y65;
Y_11bus.Y11_3 = Y_11bus.Y3_11;
Y_11bus.Y11_10 = Y_11bus.Y10_11;
Y_11bus.Y11_11 = Y_11bus.Y55;

% Вычисление проводимостей элементов 6-ти узловой сети:

tmp_12 = Y_11bus.Y55 * Y_11bus.Y66 - Y_11bus.Y56...
    * Y_11bus.Y65;
tmp_34 = Y_11bus.Y10_10 * Y_11bus.Y11_11 ...
    - Y_11bus.Y10_11 * Y_11bus.Y11_10;

Y_6bus.Y11 = Y_11bus.Y11 - (Y_11bus.Y15...
    * Y_11bus.Y51 * Y_11bus.Y66) / tmp_12;
Y_6bus.Y12 = Y_11bus.Y15 * Y_11bus.Y62 ...
    * Y_11bus.Y56 / tmp_12;
Y_6bus.Y15 = Y_11bus.Y15 * Y_11bus.Y56 ...
    * Y_11bus.Y67 / tmp_12;
Y_6bus.Y21 = Y_6bus.Y12;
Y_6bus.Y22 = Y_11bus.Y22 - (Y_11bus.Y26 ...
    * Y_11bus.Y62 * Y_11bus.Y55) / tmp_12;
Y_6bus.Y25 = -Y_11bus.Y26 * Y_11bus.Y55 ...
    * Y_11bus.Y67 / tmp_12;
Y_6bus.Y33 = Y_11bus.Y33 - (Y_11bus.Y11_3 ...
    * Y_11bus.Y3_11 * Y_11bus.Y10_10) / tmp_34;
Y_6bus.Y34 = Y_11bus.Y10_4 * Y_11bus.Y3_11 ...
    * Y_11bus.Y11_10 / tmp_34;
Y_6bus.Y36 = Y_11bus.Y3_11 * Y_11bus.Y10_9 ...
    * Y_11bus.Y11_10 / tmp_34;
Y_6bus.Y43 = Y_6bus.Y34;
Y_6bus.Y44 = Y_11bus.Y44 - (Y_11bus.Y10_4 ...
    * Y_11bus.Y4_10 * Y_11bus.Y11_11) / tmp_34;
Y_6bus.Y46 = -Y_11bus.Y4_10 * Y_11bus.Y10_9 ...
    * Y_11bus.Y11_11 / tmp_34;
Y_6bus.Y51 = Y_6bus.Y15;
Y_6bus.Y52 = Y_6bus.Y25;
Y_6bus.Y55 = Y_11bus.Y77 - (Y_11bus.Y78...
    * Y_11bus.Y87 / Y_11bus.Y88) - (Y_11bus.Y55...
    * Y_11bus.Y67 * Y_11bus.Y76 / tmp_12);
Y_6bus.Y56 = - Y_11bus.Y78 * Y_11bus.Y89 ...
    / Y_11bus.Y88;
Y_6bus.Y63 = Y_6bus.Y36;
Y_6bus.Y64 = Y_6bus.Y46;
Y_6bus.Y65 = Y_6bus.Y56;
Y_6bus.Y66 = Y_11bus.Y99 - (Y_11bus.Y89...
    * Y_11bus.Y98 / Y_11bus.Y88) - (Y_11bus.Y10_9...
    * Y_11bus.Y9_10 * Y_11bus.Y11_11 / tmp_34);

% Матрица проводимости для сети из 6-ти узлов в комплексном виде в системе
% относительных единиц линий:
Y_cmplx_l = [Y_6bus.Y11 Y_6bus.Y12 0 0 Y_6bus.Y15 0;
        Y_6bus.Y21 Y_6bus.Y22 0 0 Y_6bus.Y25 0;
        0 0 Y_6bus.Y33 Y_6bus.Y34 0 Y_6bus.Y36;
        0 0 Y_6bus.Y43 Y_6bus.Y44 0 Y_6bus.Y46;
        Y_6bus.Y51 Y_6bus.Y52 0 0 Y_6bus.Y55...
        Y_6bus.Y56;
        0 0 Y_6bus.Y63 Y_6bus.Y64 Y_6bus.Y65...
        Y_6bus.Y66];
    
% Матрица Y_cmplx_l с действительными компонентами (в системе линий):
for i=1:6                
 for j=1:6
   Y_l(2*i-1, 2*j-1) = real(Y_cmplx_l(i,j));
   Y_l(2*i-1, 2*j) = -imag(Y_cmplx_l(i,j));
   Y_l(2*i, 2*j-1) = imag(Y_cmplx_l(i,j));
   Y_l(2*i, 2*j) = real(Y_cmplx_l(i,j));
 end
end

% Матрица проводимости с действительными компонентами в системе машин:
Y_m = Y_l * 1/9;

admit_mat.Y_cmplx_l = Y_cmplx_l;
admit_mat.Y_l = Y_l;
admit_mat.Y_m = Y_m;
end
%=======================END of ADMITTANCE MATRIX===========================


%========================POWER FLOW ANALYSIS===============================
function [power_flow] = fpower_flow(mac_dat, Y)
% Количество шагов итерационной процедуры:
val_step = 1000000;

% Номер slack bus:
N_sb = 3; % код не допускает изменения номера шины (нужно исправить)

% Стартовые данные итерационной процедуры:
P = [mac_dat.P;-mac_dat.P_load] / 100; % Активная мощность в узлах 1 - 6
                                      % в относительных единицах системы 
                                      % линий 100MVA/230kV
P(N_sb) = 0; % начальная активная мощность на slack bus
Q = zeros(6,val_step); % Реактивная мощность в узлах 1-6
Q(5,:) = mac_dat.Q_load(1) / 100; % заданная реактивная мощность нагрузки
                                  % в узле № 7 в системе линий 100MVA/230kV
Q(6,:) = mac_dat.Q_load(2) / 100; % заданная реактивная мощность нагрузки
                                  % в узле № 9 в системе линий 100MVA/230kV
V_rect_l = zeros(6,val_step);
V_rect_l (1:4,1) = mac_dat.absE; % заданные амплитуды напряжений в 
                                 % генераторных узлах 1-4 в относительных
                                 % единицах системы линий 100MVA/230kV
V_rect_l(N_sb,:) = mac_dat.absE(N_sb) * exp(mac_dat.angE(N_sb)...
    * 1i * pi / 180); % комплексное напряжение на slack bus в относительных
                      % единицах системы линий 100MVA/230kV
V_rect_l (5:6,1) = 1; % начальные значения напряжений в нагрузочных узлах 
                      % № 7,9 в pu системы линий 100MVA/230kV

% Итерационная процедура:

%acc = 1.7; % "коэфф. ускорения" (1.4 ... 1.7)

for iter = 1:val_step-1
    % Напряжение в узле № 9:
    V_rect_l(6,iter + 1) = - 1 / Y(6,6) ...
    * (V_rect_l(3,iter) * Y(6,3) + V_rect_l(4,iter) * Y(6,4) ...
    + V_rect_l(5,iter) * Y(6,5) - (P(6) - Q(6,iter) * 1i) ...
    / conj(V_rect_l(6,iter)));

%     V_rect_l(6,iter + 1) = V_rect_l(6,iter) + acc * (V_rect_l(6,iter + 1)...
%         - V_rect_l(6,iter));

    % Напряжение в узле № 7:
    V_rect_l(5,iter + 1) = - 1 / Y(5,5) ...
    * (V_rect_l(1,iter) * Y(5,1) + V_rect_l(2,iter) * Y(5,2) ...
    + V_rect_l(6,iter+1) * Y(5,6) - (P(5) - Q(5,iter) * 1i) ...
    / conj(V_rect_l(5,iter)));

%     V_rect_l(5,iter + 1) = V_rect_l(5,iter) + acc * (V_rect_l(5,iter + 1)...
%             - V_rect_l(5,iter));

    % Узел № 4:
    % Реактивная мощность:
    Q(4,iter+1) = - imag(conj(V_rect_l(4,iter)) * (V_rect_l(3,iter) * Y(4,3) ...
        + V_rect_l(4,iter) * Y(4,4) + V_rect_l(6,iter+1) * Y(4,6)));
    
    % Напряжение:
    V_rect_l(4,iter + 1) = - 1 / Y(4,4) ...
    * (V_rect_l(3,iter) * Y(4,3) + V_rect_l(6,iter+1) * Y(4,6) ...
    - (P(4) - Q(4,iter+1) * 1i) / conj(V_rect_l(4,iter)));

    % Действительная часть напряжения:
    Re_V4 = sqrt (V_rect_l(4,1)^2 - imag(V_rect_l(4,iter + 1))^2);
    % Напряжение с модулем V_rect_l(4,1):
    V_rect_l(4,iter + 1) = Re_V4 + imag(V_rect_l(4,iter + 1)) * 1i;
    
%     V_rect_l(4,iter + 1) = V_rect_l(4,iter) + acc * (V_rect_l(4,iter + 1)...
%             - V_rect_l(4,iter));
    
    % Узел № 2:
    % Реактивная мощность:
    Q(2,iter+1) = - imag(conj(V_rect_l(2,iter)) * (V_rect_l(1,iter) * Y(2,1) ...
        + V_rect_l(2,iter) * Y(2,2) + V_rect_l(5,iter+1) * Y(2,5)));
    
    % Напряжение:
    V_rect_l(2,iter + 1) = - 1 / Y(2,2) ...
    * (V_rect_l(1,iter) * Y(2,1) + V_rect_l(5,iter+1) * Y(2,5) ...
    - (P(2) - Q(2,iter+1) * 1i) / conj(V_rect_l(2,iter)));

    % Действительная часть напряжения:
    Re_V2 = sqrt (V_rect_l(2,1)^2 - imag(V_rect_l(2,iter + 1))^2);
    % Напряжение с модулем V_rect_l(2,1):
    V_rect_l(2,iter + 1) = Re_V2 + imag(V_rect_l(2,iter + 1)) * 1i;
    
%     V_rect_l(2,iter + 1) = V_rect_l(2,iter) + acc * (V_rect_l(2,iter + 1)...
%             - V_rect_l(2,iter));
    
    % Узел № 1:
    % Реактивная мощность:
    Q(1,iter+1) = - imag(conj(V_rect_l(1,iter)) * (V_rect_l(1,iter) * Y(1,1) ...
        + V_rect_l(2,iter+1) * Y(1,2) + V_rect_l(5,iter+1) * Y(1,5)));
    
    % Напряжение:
    V_rect_l(1,iter + 1) = - 1 / Y(1,1) ...
    * (V_rect_l(2,iter+1) * Y(1,2) + V_rect_l(5,iter+1) * Y(1,5) ...
    - (P(1) - Q(1,iter+1) * 1i) / conj(V_rect_l(1,iter)));

    % Действительная часть напряжения:
    Re_V1 = sqrt (V_rect_l(1,1)^2 - imag(V_rect_l(1,iter + 1))^2);
    % Напряжение с модулем V_rect_l(1,1):
    V_rect_l(1,iter + 1) = Re_V1 + imag(V_rect_l(1,iter + 1)) * 1i;
    
%     V_rect_l(1,iter + 1) = V_rect_l(1,iter) + acc * (V_rect_l(1,iter + 1)...
%             - V_rect_l(1,iter));
end

% Мощности в узле № 3 (slack bus):
P(3) = real(conj(V_rect_l(3,end)) * (Y (3,:) * V_rect_l(:,end)));
Q(3,end) = -imag(conj(V_rect_l(3,end)) * (Y (3,:) * V_rect_l(:,end)));

% Установившиеся значения мощностей, напряжений и токов (в системе линий 100MVA/230kV):
P0_l = P;
Q0_l = Q(:,end);
V0_l = V_rect_l(:,end);
I0_l = (P0_l - Q0_l * 1i) ./ conj(V0_l);

% Проверка правильности расчёта power flow:
delta = I0_l - Y * V0_l(:,end); % Разница между левой и правой частью уравнений сети
angEt = rad2deg(angle(V0_l)); %Фазы напряжений на генераторных шинах

power_flow.P0_l = P0_l;
power_flow.Q0_l = Q0_l;
power_flow.V0_l = V0_l;
power_flow.I0_l = I0_l;
power_flow.delta = delta;
power_flow.angEt = angEt;
end
%======================END of POWER FLOW ANALYSIS==========================


%====================BILINEARIZE of GEN FROM PUBLICATION===================
function [bilin_gen] = fbilin_gen(mac_dat,power_flow)

n_gen = mac_dat.n_gen;
bilin_gen = struct;

for ngen = 1:n_gen

    name_field = string(['G' num2str(ngen)]);

    % Исходные нерасчётные данные в системе машин 900MVA/20kV
    omega_0 = mac_dat.omega_0(ngen);
    X_d = mac_dat.X_d (ngen);
    X_d_ht = mac_dat.X_d_ht(ngen);
    X_d_2ht = mac_dat.X_d_2ht(ngen);
    X_q = mac_dat.X_q(ngen);
    X_q_ht = mac_dat.X_q_ht(ngen);
    X_q_2ht = mac_dat.X_q_2ht(ngen);
    T_d0_ht = mac_dat.T_d0_ht(ngen);
    T_d0_2ht = mac_dat.T_d0_2ht(ngen);
    T_q0_ht = mac_dat.T_q0_ht(ngen);
    T_q0_2ht = mac_dat.T_q0_2ht(ngen);
    R_a = mac_dat.R_a(ngen);
    X_l = mac_dat.X_l(ngen);
    A_Sat = mac_dat.A_Sat(ngen);
    B_Sat = mac_dat.B_Sat(ngen);
    psi_TI = mac_dat.psi_TI(ngen);
    H = mac_dat.H(ngen);
    K_D = mac_dat.K_D(ngen);

    % Установившиеся значения мощностей и напряжений в системе линий (почему не в системе машин 900MVA/20kV???):
    P = power_flow.P0_l(ngen);%* 100 / 900; % активная мощность
    Q = power_flow.Q0_l(ngen);%* 100 / 900; % реактивная мощность
    E_t_tild = power_flow.V0_l(ngen);%* 230 / 20; % комплексное напряжение
    E_t = abs(E_t_tild); % амплитуда терминального напряжения в системе машин

    % Определяем индуктивности (предполагается, что в системе относительных
    % едениц индуктивность равна реактивному сопротивлению,X_index = L_index):
    L_d = X_d;  L_d_ht = X_d_ht;    L_d_2ht = X_d_2ht;
    L_q = X_q;  L_q_ht = X_q_ht;    L_q_2ht = X_q_2ht;
    L_l = X_l;

    % unsaturated mutual inductances (Example 4.1, p.154):
    L_ad = L_d - L_l;
    L_aq = L_q - L_l;

    % field winding inductance (Example 4.1, p.154):
    L_fd = L_ad*(L_d_ht - L_l)/(L_ad - L_d_ht + L_l);

    % damping windings inductance (Example 4.1, p.154):
    L_1q = L_aq*(L_q_ht - L_l)/(L_aq - L_q_ht + L_l);
    L_1d = -L_ad * L_fd * (L_l - L_d_2ht) / (L_ad * L_fd + L_ad * L_l - ...
        L_ad * L_d_2ht + L_fd * L_l - L_fd * L_d_2ht);
    L_2q = -L_aq * L_1q * (L_l - L_q_2ht)/(L_1q * L_aq + L_1q * L_l -...
        L_1q * L_q_2ht + L_aq * L_l - L_aq * L_q_2ht);

    % rotor windings inductance (Example 4.1, p.154):
    R_fd = (L_ad + L_fd) / (T_d0_ht * omega_0);
    R_1d = (L_1d + (L_ad * L_fd) / (L_ad + L_fd)) / (T_d0_2ht * omega_0);
    R_1q = (L_aq + L_1q) / (T_q0_ht * omega_0);
    R_2q = (L_2q + (L_aq * L_1q) / (L_aq + L_1q)) / (T_q0_2ht * omega_0);

    % total saturation factors K_sd, K_sq (Section 3.8, p.114):
    I_t_tild = conj(P + Q * 1i) / conj(E_t_tild); %p.733
    
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
I_t_tild = I_t_tild* 1/9;%100*20/(230*900);
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    E_a_tild = E_t_tild + (R_a + X_l * 1i) * I_t_tild; %eq.3.194, p.116
    psi_at = abs(E_a_tild);
    if psi_at <= psi_TI
        psi_I = 0; %eq.3.188, p.116
    else
        psi_I = A_Sat * exp(B_Sat * (psi_at-psi_TI)); %eq.3.189, p.116
    end
    K_sd = psi_at / (psi_at + psi_I); %eq.3.187, p.116
    K_sq = K_sd;

    % Параметры рабочей точки с. 746:
    % saturated values of synchronous reactances & inductances:
    X_ds = K_sd * L_ad + L_l;
    X_qs = K_sq * L_aq + L_l;
    L_ds = X_ds;
    L_qs = X_qs;

    % terminal current & power factor angle:
    I_t = abs(I_t_tild');
    f = angle(E_t_tild) - angle(I_t_tild);% (angle(I_t_tild'));

    % internal angle:
    delta_i = rad2deg (atan((I_t * X_qs * cos((f)) -...
        I_t * R_a * sin((f))) / (E_t + I_t * R_a * cos((f))...
        + I_t * X_qs * sin((f)))));

    % rotor voltages & currents:
    e_d0 = E_t * sin (deg2rad(delta_i));
    e_q0 = E_t * cos (deg2rad(delta_i));
    i_d0 = I_t * sin (deg2rad(delta_i) + f);
    i_q0 = I_t * cos (deg2rad(delta_i) + f);

    % bus voltage:
    E_re0 = real (E_t_tild);
    E_im0 = imag (E_t_tild);

    % rotor angle:
    delta_0 = delta_i + rad2deg(atan(E_im0 / E_re0));

    % field winding current & voltage:
    i_fd0 = (e_q0 + R_a * i_q0 + L_ds * i_d0) / (K_sd * L_ad);
    E_fd0 = L_ad * i_fd0;

    % mutual flux linkage:
    psi_ad0 = K_sd * L_ad * (-i_d0 + i_fd0);
    psi_aq0 = - K_sq * L_aq * i_q0;
    psi_fd0 = i_fd0 * L_fd;
    %psi_fd0 = (K_sd * L_ad + L_fd) * i_fd0 - K_sd * L_ad * i_d0;% с.104
    %psi_1d0 = K_sd * L_ad * (i_fd0 - i_d0);% с.104
    %psi_1q0 = -K_sq * L_aq * i_q0;% с.104
    %psi_2q0 = -K_sq * L_aq * i_q0;% с.104

    % Дополнительное насыщение (с.745)
    psi_at0 = psi_at;
    % incremental saturation factor:
    K_sd_incr = 1 / (1 + B_Sat * A_Sat * exp(B_Sat * (psi_at0 - psi_TI)));
    K_sq_incr = K_sd_incr;

    % incremental saturated values of the mutual inductances:
    L_ads = K_sd_incr * L_ad;
    L_aqs = K_sq_incr * L_aq;

    % Коэффициенты a_ij

    % Пересчитываем индуктивности с учётом доп. насыщения:
    L_ads_2ht = 1 / (1 / L_ads + 1 / L_fd + 1 / L_1d);
    L_aqs_2ht = 1 / (1 / L_aqs + 1 / L_1q + 1 / L_2q);
    
    % Вычисляем параметры рабочей точки с учётом L_ads_2ht и L_aqs_2ht:
    psi_1d0 = (psi_ad0 / L_ads_2ht + i_d0 - psi_fd0/L_fd) * L_1d;
    psi_1q0 = (psi_aq0 / L_aqs_2ht + i_q0) * L_1q * L_2q / (L_1q + L_2q);
    psi_2q0 = psi_1q0;
    
    delta_0r = deg2rad(delta_0);

    % Auxiliary coefficients:
    L_ds = L_ads_2ht + L_l;
    L_qs = L_aqs_2ht + L_l;
    K_LR = L_ds * L_qs + R_a^2;
    
    K_Rd = R_a * L_ads_2ht * psi_1d0 / L_1d + R_a * L_ads_2ht * psi_fd0 / L_fd;
    K_Ld = L_qs / R_a * K_Rd;
    K_Ld_tild = L_ds / R_a * K_Rd;
    
    K_Rq = R_a * L_aqs_2ht * psi_1q0 / L_1q + R_a * L_aqs_2ht * psi_2q0 / L_2q;
    K_Lq = L_ds / R_a * K_Rq;
    K_Lq_tild = L_qs / R_a * K_Rq;
    
    alfa = (L_ads_2ht - L_aqs_2ht) / K_LR;
    betta = 1 / (4 * K_LR * H);
    K1 = 2 * R_a^2 - K_LR;
    K2 = R_a * (L_ds + L_qs);
    K3 = K1 * (E_im0^2 - E_re0^2) - 2 * E_im0 * E_re0 * K2;
    K4 = K1 * E_im0 * E_re0 + 0.5 * K2 * (E_im0^2 - E_re0^2);
    K5 = L_fd / L_1d;
    K6 = L_1q / L_2q;
    K7 = R_1d / R_fd;
    K8 = R_2q / R_1q;
    gamma_d = 2 * L_ads_2ht * betta / L_fd;
    gamma_q = 2 * L_aqs_2ht * betta / L_1q;
    
    
    % A_G(i)
    a_11 = - K_D / (2 * H);
    a_12 = 2 * alfa * betta * (K3 * (1 - cos(delta_0r)^2) + 4 * K4 * sin(delta_0r)...
        * cos(delta_0r) + ((K_Ld - K_Rq)*(L_ds * e_q0 + R_a * e_d0) + ...
        (K_Lq + K_Rd) * (L_qs * e_d0 - R_a * e_q0))) + 2 * betta * ((K_Lq_tild...
        - K_Rd) * e_d0 - (K_Ld_tild + K_Rq) * e_q0);
    a_13 = 2 * betta / L_fd * (alfa * ((L_qs * K_Lq - K1 * e_d0 - 2 * L_qs * R_a * e_q0...
        - K_Rq * R_a) * L_ads_2ht + 2 * K_Ld / (psi_1d0 / L_1d + psi_fd0 / L_fd) * K_Rd)...
        + L_ads_2ht * (e_q0 * R_a - L_ds * e_d0 - K_Lq - 2 * K_Rd + K_Lq_tild));
    a_14 = K5 * a_13;
    a_15 = 2 * betta / L_1q * (alfa * ((K1 * e_q0 - (2 * L_ds * e_d0 + K_Rd) * R_a...
        + L_ds * K_Ld) * L_aqs_2ht - 2 * K_Lq * (psi_1q0/L_1q + psi_2q0/L_2q) * K_Rq)...
        + L_aqs_2ht * (K_Ld - L_qs * e_q0 - e_d0 * R_a - 2 * K_Rq - K_Ld_tild));
    a_16 = K6 * a_15;

    a_21 = omega_0;
    a_22 = 0; a_23 = 0; a_24 = 0; a_25 = 0; a_26 = 0; b_21 = 0; b_22 = 0;

    a_31 = 0;
    a_32 = omega_0 * R_fd * L_ads_2ht * (R_a * e_q0 - L_qs * e_d0) / (L_fd * K_LR);
    a_33 = omega_0 * R_fd * ((L_ads_2ht - L_fd) * K_LR - L_ads_2ht^2 * L_qs)...
        / (L_fd^2 * K_LR);
    a_34 = omega_0 * R_fd * L_ads_2ht*(K_LR - L_ads_2ht * L_qs) / (K_LR * L_fd * L_1d);
    a_35 = omega_0 * R_fd * L_ads_2ht * R_a * L_aqs_2ht / (K_LR * L_fd * L_1q);
    a_36 = K6 * a_35;

    a_41 = 0;
    a_42 = K5 * K7 * a_32;
    a_43 = K7 * a_34;
    a_44 = omega_0 * R_1d * ((L_ads_2ht - L_1d) * K_LR - L_ads_2ht^2 * L_qs)...
        / (L_1d^2 * K_LR);
    a_45 = K5 * K7 * a_35;
    a_46 = K6 * a_45;

    a_51 = 0;
    a_52 = -omega_0 * R_1q * L_aqs_2ht * (R_a * e_d0 + L_ds * e_q0) / (L_1q * K_LR);
    a_53 = - R_1q / R_fd * a_35;
    a_54 = - R_1q / R_1d * a_45;
    a_55 = omega_0 * R_1q * ((L_aqs_2ht - L_1q) * K_LR - L_aqs_2ht^2 * L_ds)...
        / (L_1q^2 * K_LR);
    a_56 = omega_0 * R_1q * L_aqs_2ht * (K_LR - L_aqs_2ht * L_ds) / (L_1q * L_2q * K_LR);

    a_61 = 0;
    a_62 = K6 * K8 * a_52;
    a_63 = -R_2q / R_fd * a_36;
    a_64 = -R_2q / R_1d * a_46;
    a_65 = K8 * a_56;
    a_66 = omega_0 * R_2q * ((L_aqs_2ht - L_2q) * K_LR - L_aqs_2ht^2 * L_ds)...
        / (L_2q^2 * K_LR);

    % B_G(i)
    b_11 = 2 * betta * (alfa * ((K1 * E_im0 - 2 * L_qs * R_a * E_re0) * ...
        (1 - 2 * cos(delta_0r)^2) + 2 * (K1 * E_re0 + K2 * E_im0) * sin(delta_0r)...
        * cos(delta_0r) + (L_ds * sin(delta_0r) - R_a * cos(delta_0r)) * ...
        (K_Ld - K_Rq) - (R_a * sin(delta_0r) + L_qs * cos(delta_0r)) * ...
        (K_Lq + K_Rd)) + (K_Rd - K_Lq_tild) * cos(delta_0r) - (K_Ld_tild + K_Rq) * sin(delta_0r));

    b_12 = 2 * betta * (alfa * ((K1 * E_re0 + 2 * L_ds * R_a * E_im0) * ...
        (1 - 2 * cos(delta_0r)^2) - 2 * (K1 * E_im0 - K2 * E_re0) * sin(delta_0r)...
        * cos(delta_0r) + (R_a * cos(delta_0r) - L_qs * sin(delta_0r)) * ...
        (K_Lq + K_Rd) - (R_a * sin(delta_0r) + L_ds * cos(delta_0r)) * ...
        (K_Ld - K_Rq)) + (K_Ld_tild + K_Rq) * cos(delta_0r) + (K_Rd - K_Lq_tild)...
        * sin(delta_0r));

    b_31 = omega_0 * R_fd * L_ads_2ht * (L_qs * cos(delta_0r) + R_a * sin(delta_0r))...
        / (L_fd * K_LR);
    b_32 = omega_0 * R_fd * L_ads_2ht * (L_qs * sin(delta_0r) - R_a * cos(delta_0r))...
        / (L_fd * K_LR);
    
    b_41 = K5 * K7 * b_31;
    b_42 = K5 * K7 * b_32;
        
    b_51 = omega_0 * R_1q * L_aqs_2ht * (R_a * cos(delta_0r) - L_ds * sin(delta_0r))...
        / (L_1q * K_LR);
    b_52 = omega_0 * R_1q * L_aqs_2ht * (L_ds * cos(delta_0r) + R_a * sin(delta_0r))...
        / (L_1q * K_LR);
    
    b_61 = K6 * K8 * b_51;
    b_62 = K6 * K8 * b_52;
    
    
    % C_G(i)
    c_11 = 0;
    c_12 = -1/K_LR * ((L_qs - L_ds) * ((2 * cos(delta_0r)^2 - 1) * E_re0 + ...
        2 * E_im0 * cos(delta_0r) * sin(delta_0r)) + (K_Lq + K_Rd) * sin(delta_0r)...
        + (K_Rq - K_Ld) * cos(delta_0r));
    c_13 = L_ads_2ht * (L_qs * sin(delta_0r) + cos(delta_0r) * R_a) / (L_fd * K_LR);
    c_14 = K5 * c_13;
    c_15 = L_aqs_2ht*(L_ds * cos(delta_0r) - sin(delta_0r) * R_a) / (L_1q * K_LR);
    c_16 = K6 * c_15;
    
    c_21 = 0;
    c_22 = 1/K_LR * ((L_qs - L_ds) * ((2 * cos(delta_0r)^2 - 1) * E_im0 - 2 * E_re0...
        * cos(delta_0r) * sin(delta_0r)) + (K_Lq + K_Rd) * cos(delta_0r) + ...
        (K_Ld - K_Rq) * sin(delta_0r));
    c_23 = L_ads_2ht * (R_a * sin(delta_0r) - L_qs * cos(delta_0r)) / (L_fd * K_LR);
    c_24 = K5 * c_23;
    c_25 = L_aqs_2ht * (L_ds * sin(delta_0r) + cos(delta_0r) * R_a) / (L_1q * K_LR);
    c_26 = K6 * c_25;
    
    
    % D_G(i)
    d_11 =  1 / K_LR * ((L_ds - L_qs) * sin(delta_0r) * cos(delta_0r) - R_a);
    d_12 = ((L_qs - L_ds) * cos(delta_0r)^2 - L_qs) / K_LR;
    
    d_21 = d_12 + (L_qs + L_ds) / K_LR;
    d_22 = -(d_11 + 2 * R_a / K_LR);
    
    
    % A2_G(i)
    a_18 = alfa * betta * (4 * K4 * (2 * cos(delta_0r)^2 - 1) + 4 * K3 * ...
        sin(delta_0r) * cos(delta_0r) + ((K_Lq + K_Rd) * (L_qs * e_q0 + ...
        R_a * e_d0) + (K_Ld - K_Rq) * (R_a * e_q0 - L_ds * e_d0))) + betta * ...
        ((K_Ld_tild + K_Rq) * e_d0 + (K_Lq_tild - K_Rd) * e_q0);
    a_19 = 2 * L_ads_2ht * betta * ((2 * alfa * L_qs - 1) * R_a * e_d0 - ...
        (alfa * K1 + L_ds) * e_q0) / L_fd;
    a_1_10 = K5 * a_19;
    a_1_11 = 2 * L_aqs_2ht * betta * ((L_qs - alfa * K1) * e_d0 - R_a * e_q0...
        * (2 * L_ds * alfa + 1))/ L_1q;
    a_1_12 = K6 * a_1_11;
    a_1_15 = 2 * betta * R_a * L_ads_2ht^2 * (alfa * L_qs - 1) / L_fd^2;
    a_1_16 = 2 * K5 * a_1_15;
    a_1_17 = 2 * L_aqs_2ht * L_ads_2ht * betta * (L_qs - L_ds - alfa * K1)...
        / (L_fd * L_1q);
    a_1_18 = K6 * a_1_17;
    a_1_22 = K5^2 * a_1_15;
    a_1_23 = K5 * a_1_17;
    a_1_24 = K5 * K6 * a_1_17;
    a_1_29 = - 2 * betta * R_a * L_aqs_2ht^2 * (alfa * L_ds + 1) / L_1q^2;
    a_1_30 = 2 * K6 * a_1_29;
    a_1_36 = K6^2 * a_1_29;
    
    a_3_8 = - omega_0 * R_fd * L_ads_2ht * (L_qs * e_q0 + e_d0 * R_a) /...
        (2 * L_fd * K_LR);
    
    a_4_8 = K5 * K7 * a_3_8;
    
    a_5_8 = omega_0 * R_1q * L_aqs_2ht * (L_ds * e_d0 - R_a * e_q0) / ...
        (2 * L_1q * K_LR);
    
    a_6_8 = K6 * K8 * a_5_8;
    
    
    % B2_G(i)
    b_11_2 = - 2 * alfa * betta * (L_qs * cos(delta_0r) + R_a * sin(delta_0r))...
        * (L_ds * sin(delta_0r) - R_a * cos(delta_0r));
    b_12_2 = 4 * alfa * betta * (K1 * (0.5 - cos(delta_0r)^2) + K2 * cos(delta_0r)...
        * sin(delta_0r));
    b_13_2 = b_12_2;
    b_14_2 = 2 * alfa * betta * (L_qs * sin(delta_0r) - R_a * cos(delta_0r))...
        * (L_ds * cos(delta_0r) + R_a * sin(delta_0r));
    
    % C2_G(i)
    c_18 = - 1/K_LR * ((L_qs - L_ds) * ((2 * cos(delta_0r)^2 - 1) * E_im0...
        - 2 * E_re0 * cos(delta_0r) * sin(delta_0r)) + 0.5 * ((K_Lq + K_Rd)...
        * cos(delta_0r) + sin(delta_0r) * (K_Ld - K_Rq)));
    c_19 = L_ads_2ht * (L_qs * cos(delta_0r) - R_a * sin(delta_0r)) / (L_fd * K_LR);
    c_1_10 = K5 * c_19;
    c_1_11 = - L_aqs_2ht * (L_ds * sin(delta_0r) + R_a * cos(delta_0r)) / (L_1q * K_LR);
    c_1_12 = K6 * c_1_11;
    
    c_28 = - 1/K_LR * ((L_qs - L_ds) * ((2 * cos(delta_0r)^2 - 1) * E_re0...
        + 2 * E_im0 * cos(delta_0r) * sin(delta_0r)) + 0.5 * ((K_Lq + K_Rd)...
        * sin(delta_0r) + cos(delta_0r) * (K_Rq - K_Ld)));
    c_29 = L_ads_2ht * (L_qs * sin(delta_0r) + R_a * cos(delta_0r)) / (L_fd * K_LR);
    c_2_10 = K5 * c_29;
    c_2_11 = L_aqs_2ht * (L_ds * cos(delta_0r) - R_a * sin(delta_0r)) / (L_1q * K_LR);
    c_2_12 = K6 * c_2_11;
    
   
    % H_G(i)
    h_12 = 2 * betta * (2 * alfa * (K1 * E_re0 + K2 * E_im0) * (2 * cos(delta_0r)^2 ...
        - 1) + 4 * alfa * (K1 * E_im0 - K2 * E_re0) * sin(delta_0r) * cos(delta_0r) +...
        alfa * ((K_Lq + K_Rd) * (L_qs * sin(delta_0r) - R_a * cos(delta_0r)) +...
        (K_Ld - K_Rq) * (R_a * sin(delta_0r) + L_ds * cos(delta_0r))) -...
        (K_Ld_tild + K_Rq) * cos(delta_0r) + (K_Lq_tild - K_Rd) * sin(delta_0r));
    h_13 = gamma_d * ((1 - 2 * alfa * L_qs) * R_a * cos(delta_0r) - (alfa *...
        K1 + L_ds) * sin(delta_0r));
    h_14 = K5 * h_13;
    h_15 = gamma_q * ((alfa * K1 - L_qs) * cos(delta_0r) - R_a * (2 * alfa ...
        * L_ds + 1) * sin(delta_0r));
    h_16 = K6 * h_15;
    h_18 = 2 * betta * (2 * alfa * (K1 * E_im0 - K2 * E_re0) * (1 - 2 * ...
        cos(delta_0r)^2) + 4* alfa * (K1 * E_re0 + K2 * E_im0) * sin(delta_0r)...
        * cos(delta_0r) + alfa * ((K_Ld - K_Rq) * (L_ds * sin(delta_0r) - R_a ...
        * cos(delta_0r)) - (K_Lq + K_Rd) * (R_a * sin(delta_0r) + L_qs *...
        cos(delta_0r))) + (K_Rd - K_Lq_tild) * cos(delta_0r) - (K_Ld_tild ...
        + K_Rq) * sin(delta_0r));
    h_19 = gamma_d * ((alfa * K1 + L_ds) * cos(delta_0r) - R_a * (2 * alfa *...
        L_qs - 1) * sin(delta_0r));
    h_1_10 = K5 * h_19;
    h_1_11 = gamma_q * ((2 * alfa * L_ds + 1) * cos(delta_0r) * R_a + (alfa ...
        * K1 - L_qs) * sin(delta_0r));
    h_1_12 = K6 * h_1_11;

    h_32 = - b_32;
    h_38 = b_31;
    
    h_42 = -b_42;
    h_48 = b_41;
    
    h_52 = -b_52;
    h_58 = b_51;
    
    h_62 = -b_62;
    h_68 = b_61;
    
    % S_G(i)
    s_12 = ((2 * cos(delta_0r)^2 - 1) * (L_ds - L_qs)) / K_LR;
    s_18 = 2 / K_LR * ((L_ds - L_qs) * cos(delta_0r) * sin(delta_0r));
    
    s_22 = s_18;
    s_28 = -s_12;
    
    
    
    % Матрица Аi Bi:

    A = [   a_11 a_12 a_13 a_14 a_15 a_16;
            a_21 a_22 a_23 a_24 a_25 a_26;
            a_31 a_32 a_33 a_34 a_35 a_36;
            a_41 a_42 a_43 a_44 a_45 a_46;
            a_51 a_52 a_53 a_54 a_55 a_56;
            a_61 a_62 a_63 a_64 a_65 a_66];
    B = [   b_11 b_12;
            b_21 b_22;
            b_31 b_32;
            b_41 b_42;
            b_51 b_52;
            b_61 b_62];
    C = [   c_11 c_12 c_13 c_14 c_15 c_16;
            c_21 c_22 c_23 c_24 c_25 c_26];
    D = [   d_11 d_12;
            d_21 d_22];

        
% Матрицы квадратичной части:
    A_2 = zeros(6,36);
    A_2(1,8:12) = [a_18 a_19 a_1_10 a_1_11 a_1_12];
    A_2(1,14:18) = [0 a_1_15 a_1_16 a_1_17 a_1_18];
    A_2(1,20:24) = [0 0 a_1_22 a_1_23 a_1_24];
    A_2(1,26:30) = [0 0 0 a_1_29 a_1_30];
    A_2(1,32:36) = [0 0 0 0 a_1_36];
    A_2(3:6,8) = [a_3_8; a_4_8; a_5_8; a_6_8];
    
    B_2 = zeros(6,4);
    B_2 (1,:) = [b_11_2 b_12_2 b_13_2 b_14_2];
    
    H_ = zeros(6,12);
    H_(1,2:6) = [h_12 h_13 h_14 h_15 h_16];
    H_(1,8:12) = [h_18 h_19 h_1_10 h_1_11 h_1_12];
    H_(3:6,2) = [h_32; h_42; h_52; h_62];
    H_(3:6,8) = [h_38; h_48; h_58; h_68];
    
    C_2 = zeros(2,36);
    C_2(:,8:12) = [ c_18 c_19 c_1_10 c_1_11 c_1_12;
                    c_28 c_29 c_2_10 c_2_11 c_2_12];
                
    S = zeros(2,12);
    S(:,2) = [ s_12; s_22];
    S(:,8) = [ s_18; s_28];
    
    % Формирует структуру результатов линеаризации генератора № ngen:
    
    bilin_gen = setfield(bilin_gen,name_field,'E_re0',E_re0);
    bilin_gen = setfield(bilin_gen,name_field,'E_im0',E_im0);
    bilin_gen = setfield(bilin_gen,name_field,'e_d0',e_d0);
    bilin_gen = setfield(bilin_gen,name_field,'e_q0',e_q0);
    bilin_gen = setfield(bilin_gen,name_field,'i_d0',i_d0);
    bilin_gen = setfield(bilin_gen,name_field,'i_q0',i_q0);
    bilin_gen = setfield(bilin_gen,name_field,'i_fd0',i_fd0);
    bilin_gen = setfield(bilin_gen,name_field,'psi_fd0',psi_fd0);
    bilin_gen = setfield(bilin_gen,name_field,'psi_1d0',psi_1d0);
    bilin_gen = setfield(bilin_gen,name_field,'psi_1q0',psi_1q0);
    bilin_gen = setfield(bilin_gen,name_field,'psi_2q0',psi_2q0);
    bilin_gen = setfield(bilin_gen,name_field,'delta_0',delta_0);
    bilin_gen = setfield(bilin_gen,name_field,'A',A);
    bilin_gen = setfield(bilin_gen,name_field,'B',B);
    bilin_gen = setfield(bilin_gen,name_field,'C',C);
    bilin_gen = setfield(bilin_gen,name_field,'D',D);
    bilin_gen = setfield(bilin_gen,name_field,'A2',A_2);
    bilin_gen = setfield(bilin_gen,name_field,'B2',B_2);
    bilin_gen = setfield(bilin_gen,name_field,'C2',C_2);
    bilin_gen = setfield(bilin_gen,name_field,'H',H_);
    bilin_gen = setfield(bilin_gen,name_field,'S',S);
    
    %lin_gen = setfield(G,name,'A',A);
    %lin_gen = setfield(G,name,'B',B);
end


end

%=================END of BILINEARIZE of GEN FROM PUBLICATION===============

%===================LINEARIZE of LOAD FROM PUBLICATION=====================

function [lin_load] = flin_load(n_gen,n_load,power_flow)

E0 = power_flow.V0_l;
P0 = power_flow.P0_l;
Q0 = power_flow.Q0_l;
lin_load = struct;
for nl = 1:n_load
    name_field = string(['D_L' num2str(nl)]);
    
    E_re = real(E0(nl+n_gen));
    E_im = imag(E0(nl+n_gen));
    Et = abs(E0(nl+n_gen));

    Di(1,1) =   P0(nl+n_gen) / (E_re^2 + E_im^2) * (1 - E_re^2 / (E_re^2 + E_im^2));

    Di(1,2) =   Q0(nl+n_gen) / (E_re^2 + E_im^2) - P0(nl+n_gen) ...
        * E_re * E_im / (E_re^2 + E_im^2)^2;

    Di(2,1) = - Q0(nl+n_gen) / (E_re^2 + E_im^2) - P0(nl+n_gen) ...
        * E_re * E_im / (E_re^2 + E_im^2)^2;

    Di(2,2) =   P0(nl+n_gen) / (E_re^2 + E_im^2) * (1 - E_im^2 / (E_re^2 + E_im^2));
    
    lin_load = setfield(lin_load,name_field,Di);
end
end

%==========================================================================


%===============CONSTRUCTION of DIAGONAL MATRICES==========================
function sys_mat = fsys_mat(n_gen, n_load,bilin_gen,lin_load)

ns = size(bilin_gen.G1.A,1); % количество переменных состояния в линейной моделе генератора
lin2mac = 1/9; % коэфф перевода из системы линий в систему машин

Ad = zeros(ns * n_gen);
Bd = zeros(ns * n_gen, 2 * (n_gen+n_load));
Cd = zeros(2 * (n_gen+n_load), ns * n_gen);
Dd = zeros(2 * (n_gen+n_load));
Ad_2 = zeros(ns * n_gen, ns * ns * n_gen);
Bd_2 = zeros(ns * n_gen, 4 * (n_gen+n_load));
Cd_2 = zeros(2 * (n_gen+n_load), ns * ns * n_gen);
Hd = zeros(ns * n_gen, 2 * ns * n_gen);
Sd = zeros(2 * (n_gen+n_load), 2 * ns * n_gen);

for k = 1:n_gen               
  name_field = string(['G' num2str(k)]);
  Ad((ns*k-ns+1):ns*k,(ns*k-ns+1):ns*k) = getfield(bilin_gen,name_field,'A');  
  Bd((ns*k-ns+1):ns*k,(2*k-1):2*k) = getfield(bilin_gen,name_field,'B');
  Cd((2*k-1):2*k,(ns*k-ns+1):ns*k) = getfield(bilin_gen,name_field,'C');
  Dd((2*k-1):2*k,(2*k-1):2*k) = getfield(bilin_gen,name_field,'D');
  Ad_2((ns*k-ns+1):ns*k,(ns^2*k-ns^2+1):ns^2*k) = getfield(bilin_gen,name_field,'A2');
  Bd_2((ns*k-ns+1):ns*k,(4*k-3):4*k) = getfield(bilin_gen,name_field,'B2');
  Cd_2((2*k-1):2*k,(ns^2*k-ns^2+1):ns^2*k) = getfield(bilin_gen,name_field,'C2');
  Hd((ns*k-ns+1):ns*k, (2*ns*k-2*ns+1):2*ns*k) = getfield(bilin_gen,name_field,'H');
  Sd((2*k-1):2*k, (2*ns*k-2*ns+1):2*ns*k) = getfield(bilin_gen,name_field,'S');
end

n = 0;
for k = n_gen+1 : n_gen+n_load
    n = n + 1;
    name_field = string(['D_L' num2str(n)]);
    Di = getfield(lin_load,name_field) * lin2mac;
    Dd((2*k-1):2*k,(2*k-1):2*k) = Di(1:2,1:2);
end
 
 sys_mat.Ad = Ad;
 sys_mat.Bd = Bd;
 sys_mat.Cd = Cd;
 sys_mat.Dd = Dd;
 sys_mat.Ad_2 = Ad_2;
 sys_mat.Bd_2 = Bd_2;
 sys_mat.Cd_2 = Cd_2;
 sys_mat.Hd = Hd;
 sys_mat.Sd = Sd;

end
%========================END of CONSTRUCTION===============================


%=============================CHECK of A===================================
function verific = fverific(A)

% Вычисляем собственные значения матрицы динамики ЭЭС:
eigA = eig(A);

% Сортируем собсвтенные значения по типу моды:
count_o = 0;
count_a = 0;
for n = 1 : length(eigA)
    if imag(eigA(n)) ~= 0
        count_o = count_o + 1;
        oscillatory(count_o,1) = eigA(n);
    else
        count_a = count_a + 1;
        aperiodic(count_a,1) = eigA(n);
    end
end

% Сортируем собственные значения по убывания действительной части:
eigAsort = sort(eigA,'descend','ComparisonMethod','real');
oscillatory_sort = sort(oscillatory,'descend','ComparisonMethod','real');
aperiodic_sort = sort(aperiodic,'descend');

% Значения колебательных и апериодических мод из книги Кундура (с.815):
eig_Kundur_osc = [  -0.76e-3+0.22e-2i;
                    -0.76e-3-0.22e-2i;
                    -0.111+3.43i;
                    -0.111-3.43i;
                    -0.492+6.82i;
                    -0.492-6.82i;
                    -0.506+7.02i;
                    -0.506-7.02i;
                    -37.89+0.142i;
                    -37.89-0.142i;
                    -38.01+0.38e-1i;
                    -38.01-0.38e-1i];

eig_Kundur_ap = [  -0.96e-1;
                    -0.117;
                    -0.265;
                    -0.276;
                    -3.428;
                    -4.139;
                    -5.287;
                    -5.303;
                    -31.03;
                    -32.45;
                    -34.07;
                    -35.53];

% Разница между собственными значениями:
delta_oscill = oscillatory_sort - eig_Kundur_osc
delta_aperiod = aperiodic_sort - eig_Kundur_ap
            
% Положение собственных значений на комплексной плоскости:            
figure()
plot(real(eigA),imag(eigA),'o', real([eig_Kundur_osc;eig_Kundur_ap]),...
    imag([eig_Kundur_osc;eig_Kundur_ap]),'x')
legend('from code','from Kundur')

verific.eigAsort = eigAsort;
verific.oscillatory_sort = oscillatory_sort;
verific.aperiodic_sort = aperiodic_sort;
verific.eig_Kundur = [eig_Kundur_osc;eig_Kundur_ap];
verific.delta_oscill = delta_oscill;
verific.delta_aperiod = delta_aperiod;
end
%============================END of CHECK==================================

