clear, clc
fs = 15; 

%% Input Data
data = xlsread('CleanData.xlsx'); 
time = data(:, 1); %time arary from 0 hours to 23 hours
L = data(:, 2); %hourly load demand from SBS Paper
E_grid = data(:, 4); %available amount of energy to import from grid
Z = 0.8; %outrage scenarios

%% Solar Parameters
I = data(:, 5); % hourly solar irradiance at time t 
M = 0.215;  %Module Efficiency
lf = 0.862; %(1-%losses) in distribution
n_solar = 1; %number of solar panels [-]; 50 previously 
c_s = 0.1255; %cost of solar per kWh generated
W_solar = 1.63; %footprint of each panel [m^2] 
g_solar = n_solar*W_solar*M*lf*I/1000; %hourly generation [kW]
S_0 = 300/1000; %******rated power of one solar panel [kW]*********
g_solar_cost = S_0*c_s; %capital cost of PV Panel [$/kW]*[kW]

%% Wind Parameters
rho = 1.2; %air density [kg/m^3]
V = data(:, 6);  %wind speed at time t [m/s]
A = 254.5; %swept area of the turbine [m^2]
P_wind = 0.5*rho*A*(V.^3)/1000; %power from wind [kW]
C_p = 0.593; %Betz limit
P_turbine = P_wind*C_p;
eta = 0.96; %gearbox transmission efficiency
P_gen = P_turbine*eta; %power generated from one turbine [kW]
n_turbines = 1; %**********number of turbine units; 50 previously***************
g_wind = n_turbines*P_gen; %total renewable energy generated from wind at time t
%W_wind = 2*2.8; %m^2 footprint of one turbine, double by 2 for spacing
W_wind = 1000;
%A_wind = W_wind*n_turbines; %total area occupied by turbines
c_w = 0.0426; %turbine cost per kWh generated [$/kWh]
W_0 = 50; %*******rated power of one turbine [kW]*********
g_wind_cost = W_0*c_w; %capital cost of turbine [$/kW]*[kW]

%% Diesel Parameters
%cap the number of diesel units at 21 
n_diesel = 1  ; %current # of diesel generator units at Puerto Rico Airport
P_diesel = 1200; %diesel generation rate [kW]
g_diesel = n_diesel*P_diesel*ones(24, 1); %power generated from diesel at time t
c_d = 0.239; %cost of diesel per kW generated
beta = 0.077818182; %diesel consumption rate [gal/kWh]
fc = 1.9; %fuel cost per gallon [$/gallon]
D_0 = P_diesel; 
g_diesel_cost = c_d; 
W_diesel = 14.7;  %[m^2] 

%% Battery Parameters
U_0 = 3.3; %Coefficient in open-circuit-voltage model [V]
C = 51782; %Equivalent Cell capacitance [F]
R = 0.01; %Cell resistance [ohm]
W_battery =  0.009; %footprint of each cell [m^2] 
B_0 = 1000; %**********rated battery capacity [kW]*********
g_battery_cost = 0.3; % LCOE of battery [$/kWh]

%% LCA Values
CO2_b = 456 *10^-6; %[ton CO2/ kWh]
CO2_s = 64 *10^-6; %[ton CO2/ kWh]
CO2_w = 31 *10^-6;%[ton CO2/ kWh]
CO2_d = 331.2 *10^-6; %[ton CO2/ kWh]
CO2_G = 844.86*10^-6 ; %[ton CO2/ kWh] Grid Lifecycle Emission Factor

% CO2_b = 456 * B_0 *10^-6; %[ton CO2/ kWh] * kW * 1h
% CO2_s = 64 * S_0 *10^-6 ; %[ton CO2/ kWh] * kW * 1h
% CO2_w = 31 * W_0 *10^-6;%[ton CO2/ kWh] * kW * 1h
% CO2_d = 331.2 * D_0 *10^-6; %[ton CO2/ kWh]
% CO2_G = 844.86*10^-6 ; %[ton CO2/ kWh] Grid Lifecycle Emission Factor

carbon_cost = 50; %social cost of carbon [$/metric ton of CO2]

%% Grid Parameters
c_grid = 0.2; %[$/kWh]

%% Scaling Parameters
%Scaling parameters were found by dividing the maximum area of each source
%by the footprints 
SOC_min = B_0*0.9; %minimum allowable cell energy levels [kWh] 
SOC_max = B_0*0.15; %maximum allowable cell energy levels [kWh] 
P_max = 270; %maximum discharge of battery [W] 
b_min = 0; %minimum battery scaling
b_max = 70*10^6; %scaling area constraint 
w_min = 0; %minimum wind scaling
w_max = 112000; %maximum wind scaling 
s_min = 0; %minimum wind scaling
s_max = 348900; %maximum solar scaling
d_min = 0; %minimum diesel scaling
d_max = 42000; %maximum diesel scaling 

A_max = 627465.2; %total area constraint [m^2]
A_used = 0.1 * A_max;

%% Optimization Program 
cvx_begin 
    variables b s w SOC(24) E(24) B_c(24) B_d(24) D(24) 
    
    minimize(sum(B_c)*g_battery_cost + sum(D)*c_d + ...
        sum(g_solar)*c_s*s + sum(g_wind)*c_w*w + ...
        c_grid*sum(E) + CO2_b*carbon_cost*b + CO2_s*sum(g_solar)*carbon_cost*s ...
        + CO2_w*sum(g_wind)*carbon_cost*w + CO2_G*carbon_cost*sum(E) + ...
        CO2_d*sum(D)*carbon_cost)
    
    subject to 
        W_battery*b + W_solar*s + W_wind*w...
            <= A_used; 
        
    for i = 1:length(time) - 1
        SOC(i+1) == SOC(i) + 0.9*B_c(i) - 1/0.9*B_d(i);
    end
    I = length(time);
    SOC(1) == SOC(I) + 0.9*B_c(I) - 1/0.9*B_d(I);
    SOC(8) == 0;
    for i = 1:length(time)
        SOC(i) <= b*SOC_max; 
        SOC(i) >= 0; 
        B_c(i) <= b*B_0;
        B_d(i) <= b*B_0;
        B_c(i) >= 0;
        B_d(i) >= 0;
        D(i) >= 0;
        D(i) <= P_diesel;
        E(i) <= E_grid*Z; 
        E(i) >= 0; 
        s*g_solar(i) + w*g_wind(i) + D(i) + B_d(i) - B_c(i) + E(i) == L(i); 
    end
    
    s <= s_max;
    s >= s_min;
    b <= b_max;
    b >= b_min;
    w <= w_max;
    w >= w_min;

cvx_end 

%% Results
fprintf(1,'------------------- Results --------------------\n');
fprintf(1,'--------------------------------------------------\n');
fprintf(1,'Minimum microgrid Cost : %4.2f USD\n',cvx_optval);
fprintf(1,'\n');
fprintf(1,'Solar %f | Wind %f | Battery %f',s,w,b);
%%
g_battery_cost*b + g_solar_cost*s + g_wind_cost*w + c_grid*sum(E)
