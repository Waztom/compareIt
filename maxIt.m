%%%%%%%%%%%%%%%%%%%% numIT to estimate kinetic coefficients %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  [values] = numIt(modelFun, no_coefficients,...
                           no_starting_points, std_dev_iterations,...
                           no_iterations, lb, ub)
    
% Set up array for random array of upper and lower bounds
random_bounds    = 0*ones(1,no_coefficients);
                                                
% Set up array sizes for storing coefficient estimates
coefficient_values       = zeros(no_iterations, no_coefficients);
coefficient_values_total = zeros(no_iterations*std_dev_iterations, no_coefficients);
median_values            = zeros(std_dev_iterations, no_coefficients);

for f = 1:std_dev_iterations
     
    for i = 1:no_iterations
 
        % Randomly select intial values from lower bound to upper bound
        % ranges for coeficient estimates
        for index = 1:no_coefficients
            random_bounds(index) = (ub(index)-lb(index)).*rand(1,1) + lb(index);
        end
        
        problem = createOptimProblem('fmincon','x0',random_bounds,'objective', modelFun,...
                                    'lb',lb,'ub',ub);
      
        ms = MultiStart('UseParallel',true,'Display','off');
        [xmulti] = run(ms,problem,no_starting_points);
         
        for value = 1:no_coefficients
            coefficient_values(i,value) = xmulti(1,value);       
        end
        
    end
    
        if f == 1
            coefficient_values_total  = [coefficient_values];
        elseif f ~= 1
            coefficient_values_total  = [coefficient_values_total;coefficient_values];
        end

     for median_value = 1:no_coefficients
            median_values(f,median_value) = median(coefficient_values(:,median_value));       
     end
end

values = {coefficient_values_total, median_values};

end