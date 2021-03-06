clear, clc
fs = 15; 

%% Input Data
data = xlsread('CleanData.xlsx'); 
time = data(:, 1); %time arary from 0 hours to 23 hours
L = data(:, 2); %hourly load demand from SBS Paper
E_grid = data(:, 4); %available amount of energy to import from grid
Z = 1;%outage scenarios !!!!parameter that can be adjusted!!!

%% Solar Parameters
I = data(:, 5); % hourly solar irradiance at time t 
T = data(:, 7); 
M = 0.215;  %Module Efficiency
lf = 0.862; %(1-%losses) in distribution
n_solar = 144; %number of solar panels [-]; 50 previously 
kt = -0.0029; %temperature coefficient 
Tref = 25; %�C
Iref = 1000; %W/m^2
c_s = 0.1255; %LCOE of solar per kWh generated
W_solar = 1.63; %footprint of each panel [m^2] 
%g_solar = n_solar*W_solar*M*lf*I/1000; %hourly generation [kW]
S_0 = 345/1000; %rated power of one solar panel [kW]
g_solar = NaN(24, 1); 
for h = 1:length(time)
    g_solar(h, 1) = n_solar*S_0*(I(h)/Iref)*(1 + kt*(T(h)+(0.0256*I(h))- Tref));
end
g_solar_cost = S_0*c_s; %capital cost of PV Panel [$/kW]*[kW]

%% Wind Parameters
rho = 1.2; %air density [kg/m^3]
V = data(:, 6);  %wind speed at time t [m/s]
A = 254.5; %swept area of the turbine [m^2]
% P_wind = 0.5*rho*A*(V.^3)/1000; %power from wind [kW]
% C_p = 0.593; %Betz limit
% P_turbine = P_wind*C_p;
% eta = 0.96; %gearbox transmission efficiency
% P_gen = P_turbine*eta; %power generated from one turbine [kW]
v_ci = 1.8; %[m/s] cut in wind speed
v_co = 28; %[m/s] cut out wind speed 
v_rt = 10; %[m/s] rated wind speed 
W_0 = 50; %rated power of one turbine [kW]
P_gen = NaN(24,1); 
%Source: https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0211642&type=printable
for h = 1:length(time) 
    if V(h) <= v_ci || V(h) >= v_co
        P_gen(h, 1) = 0 ;
    elseif V(h) >= v_ci && V(h) <= v_rt 
        P_gen(h, 1) = W_0*((V(h)-v_ci)/(v_rt-v_ci)); 
    elseif V(h) >= v_rt && V(h) <= v_co 
        P_gen(h, 1) = W_0;
    end
end
n_turbines = 1; %**********number of turbine units; 50 previously***************
g_wind = n_turbines*P_gen; %total renewable energy generated from wind at time t
%W_wind = 2*2.8; %m^2 footprint of one turbine, double by 2 for spacing
W_wind = A*2;
%A_wind = W_wind*n_turbines; %total area occupied by turbines
c_w = 0.038; %turbine cost per kWh generated [$/kWh]
g_wind_cost = W_0*c_w; %capital cost of turbine [$/kWh]*[kW]*[1hr]

%% Diesel Parameters
%cap the number of diesel units at 21 
n_diesel = 21  ; %current # of diesel generator units at Puerto Rico Airport
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
W_battery =  14.86 * 0.01; %footprint of each cell [m^2] !!Assumption!!
B_0 = 70; %rated battery capacity [kW]
%B_0 = 730; %rated battery capacity [kW]
g_battery_cost = 0.7; % LCOE of battery [$/kWh]
gamma = 0.9; 
i_min = -35; %[A]
i_max = 70; %[A]

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
SOC_min = B_0*0.15*1000; %minimum allowable cell energy levels [Wh] 
SOC_max = B_0*0.9*1000; %maximum allowable cell energy levels [Wh] 
P_max = 270; %maximum discharge of battery [W] 
%b_min = 0; %minimum battery scaling
%b_max = 100000; %scaling area constraint 
w_min = 0; %minimum wind scaling
w_max = 112000; %maximum wind scaling 
s_min = 0; %minimum wind scaling
s_max = 2656; %maximum solar scaling
d_min = 0; %minimum diesel scaling
d_max = 42000; %maximum diesel scaling 

A_max = 627465.2; %total area constraint [m^2]
A_used = A_max * 1; %!!!!parameter that can be adjusted!!!

%% Optimization Program 
cvx_precision best
cvx_begin  
    variables s w SOC(24) E(24) P(24) P_c(24) D(24) n_battery
    minimize(norm(P)/1000*g_battery_cost + sum(D)*c_d + n_battery*0.01 +...
        sum(g_solar)*c_s*s + sum(g_wind)*c_w*w + ...
        c_grid*sum(E) + CO2_b*carbon_cost*norm(P/1000) + CO2_s*sum(g_solar)*carbon_cost*s ...
        + CO2_w*sum(g_wind)*carbon_cost*w + CO2_G*carbon_cost*sum(E) + ...
        CO2_d*sum(D)*carbon_cost)
    
    subject to 
        W_battery*n_battery + W_solar*s + W_wind*w...
            <= A_used; 
    SOC(11) == 0; % SOC=0 at 11AM gives us the lowest cost

     I = length(time);

    SOC(1) == SOC(I) - 1*P_c(I); 
     %Battery energy level is the same at hour 0 and hour 24
    
    P(I) - P_c(I) + R*C*quad_over_lin(P_c(I), 2*SOC(I) +...
            C*n_battery*U_0^2) <= 0; 
     %Battery energy level at hr 24 needs to comply the rule
   
    
    for i = 1:length(time) - 1
        SOC(i+1) == SOC(i) - 1*P_c(i); 
        P(i) - P_c(i) + R*C*quad_over_lin(P_c(i), 2*SOC(i) +...
            C*n_battery*U_0^2) <= 0; 
    end
    
   M_c = [15:18]; % battery charging time
   M_d = [1:14,19:24]; % battery discharging time 
    
   for i = M_c
        P(i) >= -n_battery*640; 
        % [W] assume power limit of each battery is 640 W 
        P(i) <= 0;
   end
   
   for i = M_d
       P_c(i) >= i_min*geo_mean([n_battery, (2/C*SOC(i)+U_0^2*n_battery)]);
       P_c(i) <= i_max*geo_mean([n_battery, (2/C*SOC(i)+U_0^2*n_battery)]); 
   end
         
    for i = 1:length(time)
        SOC(i) <= n_battery*SOC_max; %[Wh] 
        SOC(i) >= 0; %[Wh]
        D(i) >= 0;
        D(i) <= (1-Z)*n_diesel*P_diesel;
        E(i) >= 0; 
        s*g_solar(i) + w*g_wind(i) + D(i) + + E(i) + P(i)/1000 == L(i); 
    end
    

    
    n_battery >= 0;
    sum(E) <= sum(E_grid)*Z; 
    s <= s_max;
    s >= s_min;
    w <= w_max;
    w >= w_min;


cvx_end 

CO2 = CO2_b*n_battery + CO2_s*sum(g_solar)*s+ CO2_w*sum(g_wind)*w + CO2_G*sum(E) + CO2_d*sum(D);

%% Results
fprintf(1,'------------------- Results --------------------\n');
fprintf(1,'--------------------------------------------------\n');
fprintf(1,'Minimum microgrid Cost : %4.2f USD\n',cvx_optval);
fprintf(1,'\n');
fprintf(1,'Solar %f | Wind %f | Battery %f',s,w,n_battery);
fprintf(1,'\n');
fprintf(1,'Total daily emissions %f metric tons CO2eq',CO2, '\n');

%Renewable Generation
total_wind = w.*g_wind;
total_solar = s.*g_solar;
figure(1)
plot(time/24, E, time/24, L, time/24, total_solar, time/24, total_wind, 'cyan',...
    time/24, D, 'black', 'LineWidth', 5)
title(['Outage Scenario (Z = ', num2str(Z), ')'])
legend('Hourly Energy Imported from Grid', 'Hourly Demand', ...
    'Solar Generation', 'Wind Generation', 'Diesel Generation','location','southeast') 
%xticks([0:1:23])
%xlim([0 24])
set(gca,'FontSize',20)
datetick('x', 'HHPM')
xlabel('Hour')
ylabel('Power [kW]')

%Battery Charging & Discharging'
figure(2) 
plot(time/24, P/1000, 'green', time/24, P_c/1000, 'red', 'LineWidth', 5) 
title(['Battery Dynamics (Z = ', num2str(Z), ')'])
legend('Battery Power ', 'Electrochemical Battery Power')
set(gca, 'FontSize', 20)
datetick('x', 'HHPM')
xlabel('Hour')
ylabel('Power [kW]')

%Battery SOC
figure(3)
plot(time/24, SOC/1000, 'LineWidth', 5)
title(['Battery SOC (Z = ', num2str(Z), ')'])
set(gca, 'FontSize', 20)
datetick('x', 'HHPM')
xlabel('Hour')
ylabel('Energy [kWh]')


