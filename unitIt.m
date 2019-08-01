%%%%%%%% unitIt to convert ppm data to unitary and area data  %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Data_ppm Data_unit Data_area Data_mol...
          Data_N2_O2_ratio X_without_time]...
          = unitIt(Experimental_settings, data, no_blanks)

% For clarity, actual values from user input array should follow
Irradiance = Experimental_settings(1); % mw/cm2
Temp       = Experimental_settings(2); % deg celcius
Argon_flow = Experimental_settings(3); % ml/min
CO2_flow   = Experimental_settings(4); % ml/min
Amount_H2O = Experimental_settings(5); % ppm volume
Pressure   = Experimental_settings(6); % kPa
Mass_cat   = Experimental_settings(7); % mg 
Area_coated = Experimental_settings(8); % cm2

Total_flow = Argon_flow + CO2_flow;

CO2_Ar_flow_after_sat = Total_flow * (1-(Amount_H2O/1E6));
H2O_flow              = Total_flow - CO2_Ar_flow_after_sat;

% Moles/min from flow rates
Mol_H2O = (Pressure * H2O_flow)/(8314 * (Temp + 273.15));
Mol_CO2 = (((CO2_Ar_flow_after_sat) * (CO2_flow/(Argon_flow + CO2_flow))...
           * 101) / (8314 * (Temp + 273.15)));
Mol_Ar  = (((CO2_Ar_flow_after_sat) * (Argon_flow/(Argon_flow + CO2_flow))...
           * 101) / (8314 * (Temp + 273.15)));
Total_flow = Mol_H2O + Mol_CO2 + Mol_Ar;

% Mole fraction for partial pressure calc
X_H2O = Mol_H2O / Total_flow;
X_CO2 = Mol_CO2 / Total_flow;
X_Ar  = Mol_Ar / Total_flow;

% Partial pressures
PH2O = X_H2O * Pressure;
PCO2 = X_CO2 * Pressure;
PAr  = X_Ar  * Pressure;

% Remove average values in blanks from CO, CH4 and H2 exp data
Data_sub_blanks        = data((no_blanks + 1):end,[2:3 5]) - ...
                             mean(data(1:no_blanks,[2:3 5]));

%Data_sub_blanks(:,1:2) = abs(Data_sub_blanks(:,1:2));
Data_sub_blanks = max(0, Data_sub_blanks);

Data_ppm               = [data((no_blanks + 1):end,1) Data_sub_blanks(:,1:2)...
                          data((no_blanks +1):end, 4) Data_sub_blanks(:,3)...
                          data((no_blanks + 1):end,6) data((no_blanks + 1):...
                                                          end,7)];

Data_unit       = Data_ppm.*((60 * Total_flow) / (Mass_cat / 1000));
Data_area       = Data_ppm.*((60 * Total_flow) / (Area_coated));
Data_mol        = Data_ppm.*((60 * Total_flow));
Data_N2_O2_ratio = Data_ppm(:,7)./Data_ppm(:,6); 

X_without_time  = ones(length(Data_unit),3);
X_without_time  = X_without_time .* [Irradiance PH2O PCO2];

end