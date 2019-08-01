%%%%%%%%%%%%%%%%%%%% histIt to plot histograms %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function histogram = histIt(coefficient_values,no_coefficients,hist_title)

histogram = figure;
sgtitle(hist_title);

Coefficients = {'k', 'k_d', '\alpha', 'K_{H_2O}', 'K_{CO_2}'};

% Set median line height  
max_height=8;  

if no_coefficients == 4
    plot_dimension_rows = 2;
    plot_dimension_cols = 2;
else
    plot_dimension_rows = 2;
    plot_dimension_cols = 3;
end
    
for i = 1:no_coefficients
    subplot(plot_dimension_rows,plot_dimension_cols,i); 
    histfit(coefficient_values(:,i),30);
    title(Coefficients{i});
    ylabel('Frequency');
    xlabel('Coefficient estimate');
    hold on
    line([median(coefficient_values(:,i)),median(coefficient_values(:,i))],...
         [0,max_height],'Color','r','LineStyle','--');
end