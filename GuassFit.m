%%%%%%%%%%%%%%%%%%% GuassFit to fit Guassian peaks %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fit gaussian to peaks identified (loc) and
% calculate area under fitted Gaussian. loc = peak location from FindPeaks 
% function. x and y are data points, r is value added or subtracted to 
% loc to find range for x and data for Gausssian fit 

function area = GuassFit(x,y)

    f = fit(x,y,'gauss2');
    y_guass = f(x);
    area = trapz(x,abs(y_guass));
end