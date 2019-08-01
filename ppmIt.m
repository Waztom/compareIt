%%%%%%%%%%%%%%%%%%%%%%% ppmIt to get ppm data %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find peaks, fit gaussian to peaks (if any) and integrate 
% to calcluate area under the curve

function ppm_data = ppmIt(data,x_min,x_max,MinPeakHeight,...
                          MinPeakDistance, MinPeakWidth,... 
                          data_fit,detector,compound)
                      
if detector == "FID"    
    ppm_data =[];
    for i = 1:numel(data)
        if data(i).instrument == "GC/FID" &&...
           data(i).MethodName ~= "Warren Short Gas_break.M"
            
            x = data(i).Time;
            y = data(i).Signal;
            
            range = x >= x_min & x <= x_max;
            x = x(range);
            y = y(range);
            
            y_normalised = Y_normalised(x,y);
            plot(x,y_normalised);
            title(strcat(detector," ",compound));
            
            % Find peaks and integrate to find ppm values
            [pks, loc]   = findpeaks(y_normalised,'MinPeakHeight',...
                                     MinPeakHeight,...
                                    'MinPeakDistance', MinPeakDistance,...
                                    'MinPeakWidth',MinPeakWidth);           
            if numel(loc) == 0
                ppm_data = [ppm_data; 0];
            elseif numel(loc) > 1;
                ppm_data = [ppm_data; 0];
            elseif numel(loc) == 1
                area = GuassFit(x,y_normalised);
                ppm_data = [ppm_data; predict(data_fit, area)];
            end
        end
    end
    
  elseif detector == "TCD"    
    ppm_data =[];
    for i = 1:numel(data)
        if data(i).instrument == "GC/TCD" &&...
           data(i).MethodName ~= "Warren Short Gas_break.M"
            
            x = data(i).Time;
            y = data(i).Signal;
            
            range = x >= x_min & x <= x_max;
            x = x(range);
            y = y(range);
            
            y_normalised = Y_normalised(x,y);
            plot(x,y_normalised)            
            title(strcat(detector," ",compound));
             
            %Find peaks and integrate to find ppm values
            [pks, loc]   = findpeaks(y_normalised,'MinPeakHeight',...
                                     MinPeakHeight,...
                                    'MinPeakDistance', MinPeakDistance,...
                                    'MinPeakWidth',MinPeakWidth);           
            if numel(loc) == 0
                ppm_data = [ppm_data; 0];
            elseif numel(loc) > 1;
                ppm_data = [ppm_data; 0];
            elseif numel(loc) == 1
                x_new = GuassFit(x,y_normalised);
                ppm_data = [ppm_data; predict(data_fit, x_new)];
            end
        end    
    end 
    
end