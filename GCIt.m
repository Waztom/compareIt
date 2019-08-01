%%%%%%%%%%%%%%%%%Import data from Agilent GC  %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%        GCIt function         %%%%%%%%%%%%%%%%%%%%%%%%%

function [y , X, calib_models] = GCIt(calib_models)

if length(calib_models) == 0
    % Import FID and TCD calibration data. Save calibration files in 
    % ascending order with the lowest number having the lowest calibration. 
    % Need to import twO separate times as different calibration values 
    % likely used due to difference in sensitivity of FID vs TCD detector
    
    disp("Please select FID calibration GC folder")
    FID_calibration_data = ImportAgilentFID();
    disp("Please select TCD calibration GC folder")
    TCD_calibration_data = ImportAgilentFID();

    % Input ppm values and number of calibrations completed. %Enter as [], 
    % so if 10 and 100 ppm used, enter '[100 10]'
    FID_user_prompt_ppm_used    = 'What ppms were used for FID calibration? '; 
    FID_calibration_ppms_used   = input(FID_user_prompt_ppm_used)';
    TCD_user_prompt_ppm_used    = 'What were the ppms used for TCD calibration? '; 
    TCD_calibration_ppms_used   = input(TCD_user_prompt_ppm_used)';

    % FID and TCD calibration area arrays
    FID_CO2_Area = [];
    FID_CH4_Area = [];
    FID_CO_Area  = [];

    TCD_H2_Area  = [];
    TCD_CO2_Area = [];
    TCD_O2_Area  = [];
    TCD_N2_Area  = [];

    % Manually select peaks for FID calibration data
    findpeaks(FID_calibration_data(1).Signal,FID_calibration_data(1).Time,...
              'MinPeakheight',10,'MinPeakDistance',0.1,'MinPeakWidth',0.01);

    % Select FID peak ranges from plot
    disp("Set range of FID peaks");
    disp("Press a key once zoomed in on the  CO2 peak");
    pause;
    [FID_CO2_range,~] = ginput(2);
    disp("Press a key once zoomed in on the CH4 peak");
    pause;
    [FID_CH4_range,~] = ginput(2);
    disp("Press a key once zoomed in on the CO peak");
    pause;
    [FID_CO_range,~] = ginput(2);

    % Manually select peaks for TCD calibration data
    findpeaks(TCD_calibration_data(2).Signal,TCD_calibration_data(2).Time,...
              'MinPeakheight',10,'MinPeakDistance',0.1,'MinPeakWidth',0.01);

    % Select TCD peak ranges from plot
    disp("Set range of TCD peaks");
    disp("Press a key once zoomed in on the  H2 peak");
    pause;
    [TCD_H2_range,~] = ginput(2);
    disp("Press a key once zoomed in on the CO2 peak");
    pause;
    [TCD_CO2_range,~] = ginput(2);
    disp("Press a key once zoomed in on the O2 peak");
    pause;
    [TCD_O2_range,~] = ginput(2);
    disp("Press a key once zoomed in on the N2 peak");
    pause;
    [TCD_N2_range,~] = ginput(2);

    % Get areas and peak locations using calibration data
    for i = 1:numel(FID_calibration_data)

        if FID_calibration_data(i).instrument == "GC/FID"

            FID_x_time_data_points      = FID_calibration_data(i).Time;
            FID_y_data_response_points  = FID_calibration_data(i).Signal;

            FID_CO2_Area = [FID_CO2_Area; GuassFitCalibration(FID_CO2_range(1),...
                            FID_CO2_range(2),FID_x_time_data_points,...
                            FID_y_data_response_points)];
            FID_CH4_Area = [FID_CH4_Area; GuassFitCalibration(FID_CH4_range(1),...
                            FID_CH4_range(2),FID_x_time_data_points,...
                            FID_y_data_response_points)];
            FID_CO_Area  = [FID_CO_Area; GuassFitCalibration(FID_CO_range(1),...
                            FID_CO_range(2),FID_x_time_data_points,...
                            FID_y_data_response_points)];
         end
    end

    for i = 1:numel(TCD_calibration_data)
         if TCD_calibration_data(i).instrument == "GC/TCD"

            TCD_x_time_data_points      = TCD_calibration_data(i).Time;
            TCD_y_data_response_points  = TCD_calibration_data(i).Signal;

            TCD_H2_Area  = [TCD_H2_Area; GuassFitCalibration(TCD_H2_range(1),...
                            TCD_H2_range(2), TCD_x_time_data_points,...
                            TCD_y_data_response_points)];
            TCD_CO2_Area = [TCD_CO2_Area; GuassFitCalibration(TCD_CO2_range(1),...
                            TCD_CO2_range(2), TCD_x_time_data_points,...
                            TCD_y_data_response_points)];
            TCD_O2_Area  = [TCD_O2_Area; GuassFitCalibration(TCD_O2_range(1),...
                            TCD_O2_range(2), TCD_x_time_data_points,...
                            TCD_y_data_response_points)];
            TCD_N2_Area  = [TCD_N2_Area; GuassFitCalibration(TCD_N2_range(1),...
                            TCD_N2_range(2), TCD_x_time_data_points,...
                            TCD_y_data_response_points)];            
         end
    end

    % Fit linear functions to FID and TCD calibration data
    Linear_model = 'y ~ x1 - 1';

    FID_CO2_Linear_Fit = fitlm(FID_CO2_Area, FID_calibration_ppms_used,...
                               Linear_model);
    FID_CH4_Linear_Fit = fitlm(FID_CH4_Area, FID_calibration_ppms_used,...
                               Linear_model);
    FID_CO_Linear_Fit  = fitlm(FID_CO_Area, FID_calibration_ppms_used,...
                               Linear_model);

    TCD_H2_Linear_Fit  = fitlm(TCD_H2_Area, TCD_calibration_ppms_used,...
                               Linear_model);
    TCD_CO2_Linear_Fit = fitlm(TCD_CO2_Area, TCD_calibration_ppms_used,...
                               Linear_model);
    TCD_O2_Linear_Fit  = fitlm(TCD_O2_Area, TCD_calibration_ppms_used,...
                               Linear_model);
    TCD_N2_Linear_Fit  = fitlm(TCD_N2_Area, TCD_calibration_ppms_used,...
                               Linear_model);

    % Plot the calibration curves
    figure;
    subplot(3,3,1);plot(FID_CO2_Linear_Fit);...
                        title('FID CO_{2} calibration curve',...
                        'Interpreter','tex');
    subplot(3,3,2);plot(FID_CH4_Linear_Fit);...
                        title('FID CH_{4} calibration curve',...
                        'Interpreter','tex');
    subplot(3,3,3);plot(FID_CH4_Linear_Fit);...
                        title('FID CO calibration curve',...
                        'Interpreter','tex');
    subplot(3,3,4);plot(TCD_H2_Linear_Fit);...
                        title('TCD H_{2} calibration curve',...
                        'Interpreter','tex');
    subplot(3,3,5);plot(TCD_CO2_Linear_Fit);...
                        title('TCD CO_{2} calibration curve',...
                        'Interpreter','tex');
    subplot(3,3,6);plot(TCD_O2_Linear_Fit);...
                        title('TCD O_{2} calibration curve',...
                        'Interpreter','tex');
    subplot(3,3,7);plot(TCD_N2_Linear_Fit);...
                        title('TCD N_{2} calibration curve',...
                        'Interpreter','tex');
    
    calib_models = {FID_CO2_Linear_Fit FID_CH4_Linear_Fit FID_CO_Linear_Fit...
                    TCD_H2_Linear_Fit TCD_CO2_Linear_Fit TCD_O2_Linear_Fit...
                    TCD_N2_Linear_Fit};
end

if ~isempty(calib_models) 
    
    % Import experimental GC data
    disp("Please select experimental GC folder")
    Experimental_data = ImportAgilentFID();

    % Manually select peaks for FID calibration data
    % Note run 15 is chosen as usually production is seen at this point, although 
    % it will change if the flow rate is changed. NB data imported in order of 
    % FID then TCD data. Odd numbers == FID data and even numbers == TCD data
    figure; findpeaks(Experimental_data(25).Signal,Experimental_data(25).Time,...
                     'MinPeakheight',10,'MinPeakDistance',0.1,...
                     'MinPeakWidth',0.01);

    % Select FID peak ranges from plot
    disp("Set range of FID peaks");
    disp("Press a key once zoomed in on the CO2 peak");
    pause;
    [FID_CO2_range,~] = ginput(2);
    disp("Press a key once zoomed in on the CH4 peak");
    pause;
    [FID_CH4_range,~] = ginput(2);
    disp("Press a key once zoomed in on the CO peak");
    pause;
    [FID_CO_range,~] = ginput(2);

    % Manually select peaks for TCD calibration data
    findpeaks(Experimental_data(16).Signal,Experimental_data(16).Time,...
              'MinPeakheight',10,'MinPeakDistance',0.1,'MinPeakWidth',0.01);

    % Select TCD peak ranges from plot
    disp("Set range of TCD peaks");
    disp("Press a key once zoomed in on the H2 peak");
    pause;
    [TCD_H2_range,~] = ginput(2);
    disp("Press a key once zoomed in on the CO2 peak");
    pause;
    [TCD_CO2_range,~] = ginput(2);
    disp("Press a key once zoomed in on the O2 peak");
    pause;
    [TCD_O2_range,~] = ginput(2);
    disp("Press a key once zoomed in on the N2 peak");
    pause;
    [TCD_N2_range,~] = ginput(2);

    % Use calibration data and area to calculate ppm for peak
    % You may need to run this again to set appropriate range for peaks 
    FID_CO2_ppm_data  = ppmIt(Experimental_data,FID_CO2_range(1),...
                        FID_CO2_range(2),200,50,50,calib_models{1},...
                        "FID","CO_{2}");
    FID_CH4_ppm_data  = ppmIt(Experimental_data,FID_CH4_range(1),...
                        FID_CH4_range(2),200,50,50,calib_models{2},...
                        "FID","CH_{4}");
    FID_CO_ppm_data   = ppmIt(Experimental_data,FID_CO_range(1),...
                        FID_CO_range(2),200,50,50,calib_models{3},...
                        "FID","CO");

    TCD_H2_ppm_data   = ppmIt(Experimental_data,TCD_H2_range(1),...
                        TCD_H2_range(2),20,10,10,calib_models{4},...
                        "TCD","H_{2}");
    TCD_CO2_ppm_data  = ppmIt(Experimental_data,TCD_CO2_range(1),...
                        TCD_CO2_range(2),20,1,10,calib_models{5},...
                        "TCD","CO_{2}"); 
    TCD_O2_ppm_data   = ppmIt(Experimental_data,TCD_O2_range(1),...
                        TCD_O2_range(2),20,10,10,calib_models{6},...
                        "TCD","O_{2}");
    TCD_N2_ppm_data   = ppmIt(Experimental_data,TCD_N2_range(1),...
                        TCD_N2_range(2),20,10,10,calib_models{7},...
                        "TCD","N_{2}");

    % Put data into final table 
    Data_ppm = [FID_CO2_ppm_data FID_CH4_ppm_data FID_CO_ppm_data...
                TCD_CO2_ppm_data TCD_H2_ppm_data TCD_O2_ppm_data TCD_N2_ppm_data];

    % Get experimental settings
    % Enter as [Irradiance(mw/cm2); Temp; Ar_flow_rate(ml/min); 
    % CO2_flow_rate(ml/min; Amount_H2O(ppm); Pressure(kPa); 
    % Mass_catalyst(mg); Area_coated(cm2)))]
    % NB for area coated, need to check mesh guage and area open - see 
    % 'https://www.meshdirect.co.uk/woven-stainless-wire-cloth-80-mesh-0.18-
    % mm-aperture.html'
    Experimental_ui       = 'Enter experimental settings '; 
    Experimental_settings = input(Experimental_ui)';

    % Convert ppm data into unitary production
    No_Blanks_ui = 'Enter number of blanks run '; 
    No_Blanks    = input(No_Blanks_ui)';

    [Data_ppm, Data_unit, Data_area,...
     Data_mol,Data_N2_O2_ratio,X_without_time] = unitIt(Experimental_settings,...
                                                        Data_ppm, No_Blanks);

    % Get date, time and method info
    [FID_date_time, FID_sample_info] = infoIt(Experimental_data,"FID");

    FID_date_time   = FID_date_time((No_Blanks + 1):end,:);
    FID_sample_info = FID_sample_info((No_Blanks +1):end,:);

    % Convert time to minutes
    Time_vec = datevec(FID_date_time-FID_date_time(1));
    Time     = Time_vec(:,4:6);
    Time_min = Time(:,1).*60 + Time(:,2) + Time(:,3)./60;
    Time_hr  = Time_min./60;

    % Construct table including sample name, time, ppm and unitary production
    y = table(FID_sample_info, FID_date_time, Data_ppm, Data_unit,...
              Data_area, Data_mol, Data_N2_O2_ratio);

    % Construct variable predictors array with equal number of data rows
    % X = repelem(Experimental_settings,size(Data_unit,1),1);
    X = [X_without_time Time_hr Time_min];
end
 
end