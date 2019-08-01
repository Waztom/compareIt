%%%%%%%%%%%%% Y_normalised to account for peaks on slopes %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Returns normalised y-data points after subtracting y values
% from line under curve 

function y_normalised = Y_normalised(x,y)
    
    % Fit straight line between two points to normalise signals on
    % slopes
    
    y_1 = y(1);
    y_2 = y(end);

    x_1 = x(1);
    x_2 = x(end);

    m = (y_1 - y_2) / (x_1 - x_2);
    c = y_2 - m * x_2;

    y_normalised     = smoothdata(y - (m * x + c));
    min_y_normalised = min(y_normalised);

    if min_y_normalised < 0
        y_normalised = y_normalised + abs(min_y_normalised);
    end   

end