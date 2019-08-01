%%%%%%%%%%%%%%% GuassFitCalibration to fit Guassian peaks %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Fit Gaussian to peaks identified (loc) 
% and calculate area under fitted Gaussian. loc = peak location from 
% FindPeaks function. x and y are data points, r is value added or 
% subtracted to loc to find range for x and data for Gausssian fit 

function area = GuassFitCalibration(peak_start,peak_end,x,y)

    range = x >= peak_start & x <= peak_end;
    x = x(range);
    y = y(range);
    
    % Normalise y values in case they are on a slope
    y_normalised = Y_normalised(x,y);
    
    %Fit Gaussian to x and y values
    f = fit(x,y_normalised,'gauss2');
    plot(f)
    y_guass = f(x);
    area = trapz(x,abs(y_guass));
end