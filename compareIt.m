%% Make list of kinetic models %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ModelCH4 = @(b,X) b(1).* (((X(:,1).^b(3))./(b(2))).*((X(:,4)./(b(2))).^((X(:,1).^b(3))-1)).*(exp(-1.*((X(:,4)./(b(2))).^(X(:,1).^b(3))))))...
           .* ((((b(4).*X(:,2)).^4).*(b(5).*X(:,3)))./((1 + b(4).*X(:,2) + b(5).*X(:,3)).^(5)));
       
ModelCO = @(b,X) b(1).* (((X(:,1).^b(3))./(b(2))).*((X(:,4)./(b(2))).^((X(:,1).^b(3))-1)).*(exp(-1.*((X(:,4)./(b(2))).^(X(:,1).^b(3))))))...
           .* ((b(4).*X(:,2).*b(5).*X(:,3))./((1 + b(4).*X(:,2) + b(5).*X(:,3)).^2));
       
ModelH2 = @(b,X) b(1).* (((X(:,1).^b(3))./(b(2))).*((X(:,4)./(b(2))).^((X(:,1).^b(3))-1)).*(exp(-1.*((X(:,4)./(b(2))).^(X(:,1).^b(3))))))...
           .* ((((b(4).*X(:,2)).^2))./((1 + b(4).*X(:,2) + b(5).*X(:,3)).^2));

%% Group model handles in cell and set upper and lower bound limts for coefficient estimates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myData    = {All_CH4_Model,            All_CO_Model,        All_H2_Model};

myModels  = {ModelCH4,                 ModelCO,             ModelH2};

lb        = {[0,0.1,0.01,0,0],   [0,0.1,0.01,0,0],    [0,0.1,0.01,0,0]}; 
ub        = {[10,100,0.5,100,100], [10,100,0.5,100,100], [10,100,0.5,100,100]};

compounds = {'CH_4', 'CO', 'H_2'};
%% Test various models and collect the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parpool;

for K = 1:length(myModels) 
tic

model_coeffients{K} = numIt(myModels{K},All_Model_X,myData{K},length(lb{K}),...
                                10,10,50,lb{K},ub{K});
                           
toc
end   

delete(gcp);
%% Get statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Each cell describes standard deviation/mean of median value/s obtained for
% each model 

sdev_multi_runs = {};
mean_multi_runs = {};

for K = 1:length(model_coeffients)
    
    sdev_multi_runs{K} = std(model_coeffients{K}{2});
    mean_multi_runs{K} = mean(model_coeffients{K}{2});
    
end

%% Plot histograms and box plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for K = 1:length(model_coeffients)
    
   histIt(model_coeffients{K}{1},size(model_coeffients{K}{1},2),compounds{K})
   
end    


%% Plot model vs data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y_units  = '(\mumol.cm^{-2}.h^{-1})';
x_units  = 'Time (h)';

title_plots= {'(a)', '(b)', '(c)'};
color_plots = {[0, 0.4470, 0.7410],[0, 0.5, 0],[0.4940, 0.1840, 0.5560]};

% Standard error
Standard_error =[];

% Data for change in CO2 
% data_range_all = {(73:141), (142:212), (1:72)};
% legend_plot={'25.72 kPa CO_2 experimental values (PCO_{2}-L)',...
%              '25.72 kPa CO_2 Weibull LH kinetic model',...
%              '48.44 kPa CO_2 experimental values (PCO_{2}-I)',...
%              '48.44 kPa CO_2 Weibull LH kinetic model',...
%              '98.38 kPa CO_2 experimental values (PCO_{2}-H)',...
%              '98.38 kPa CO_2 Weibull LH kinetic model'};

% Data for change in H2O
data_range_all = {(73:141), (213:281), (282:350)};
legend_plot={'2.66 kPa H_2O experimental values (PH_{2}O-L)',...
             '2.66 kPa H_2O Weibull LH kinetic model',...
             '3.47 kPa H_2O experimental values (PH_{2}O-I)',...
             '3.47 kPa H_2O Weibull LH kinetic model',...
             '6.64 kPa H_2O experimental values (PH_{2}O-H)',...
             '6.64 kPa H_2O Weibull LH kinetic model'};

for K = 1:length(model_coeffients)
   
   figure;
  
   for data = 1:length(data_range_all)     
            
        b = mean_multi_runs{K};
        data_range = data_range_all{data};
        
        y_pred = myModels{K}(b,All_Model_X(data_range,:));
        time = All_Model_X(data_range,4);
        plot(time,myData{K}(data_range),'LineStyle','none','MarkerSize',10,'MarkerEdgeColor',color_plots{data},...
             'MarkerFaceColor',color_plots{data},'Marker','d');
        hold on;
        plot(time,y_pred, 'Linewidth',3,'Color',color_plots{data});
        
        std_err = sqrt((sum((myData{K}(data_range)-y_pred).^2))/(numel(myData{K})));
        Standard_error = [Standard_error std_err];
     

   end
   
   legend(legend_plot)
    
   xlabel(x_units);
   ylabel([compounds{K} ' production ' y_units]);
   title(title_plots{K});

   set(findall(gcf,'-property','FontSize','-property', 'FontWeight'),'FontSize',...
              25, 'Fontweight','bold');

 

end

%% Get product distribution plots
y_units  = {'CH_4 product distribution (%)', 'CO product distribution (%)', 'H_2 product distribution (%)'};
y_units_S  = {'CH_4/CO selectivity', 'CO/H_2 selectivity', 'CH_{4}/H_2 selectivity'};
x_units  = 'Time (h)';

title_plots= {'(a)', '(b)', '(c)'};
color_plots = {[0, 0.4470, 0.7410],[0, 0.5, 0],[0.4940, 0.1840, 0.5560]};

% Data for change in CO2 
% data_range_all = {(73:141), (142:212), (1:72)};
% legend_plot_S = {'25.72 kPa CO_2 selectivity (PCO_{2}-L)',...
%                  '48.44 kPa CO_2 selectivity (PCO_{2}-I)',...
%                  '98.38 kPa CO_2 selectivity (PCO_{2}-H)'};
%          
% legend_plot = {'25.72 kPa CO_2 product distribution (PCO_{2}-L)',...
%                '48.44 kPa CO_2 product distribution (PCO_{2}-I)',...
%                '98.38 kPa CO_2 product distribution (PCO_{2}-H)'};

% Data for change in H2O
data_range_all = {(73:141), (213:281), (282:350)};
legend_plot_S = {'2.66 kPa H_2O selectivity (PH_{2}O-L)',...
                 '3.47 kPa H_2O selectivity (PH_{2}O-I)',...
                 '6.64 kPa H_2O selectivity (PH_{2}O-H)'};

legend_plot = {'2.66 kPa H_2O product distribution (PH_{2}O-L)',...
               '3.47 kPa H_2O product distribution (PH_{2}O-I)',...
               '6.64 kPa H_2O product distribution (PH_{2}O-H)'};
         
CH4_product_distribution = myData{1}./(myData{1} + myData{2} + myData{3}).*100;
CO_product_distribution  = myData{2}./(myData{1} + myData{2} + myData{3}).*100;
H2_product_distribution  = myData{3}./(myData{1} + myData{2} + myData{3}).*100;

product_distribution = {CH4_product_distribution, CO_product_distribution, H2_product_distribution};

selectivity_CH4_CO = myData{1}./myData{2};
selectivity_CO_H2  = myData{2}./myData{3};
selectivity_CH4_H2 = myData{1}./myData{3};

selectivity_all = {selectivity_CH4_CO, selectivity_CO_H2, selectivity_CH4_H2};

for S = 1:length(selectivity_all)
   
    figure;
       
    for data = 1:length(data_range_all)     
        
        select = selectivity_all{S};
        data_range = data_range_all{data};
        
        time = All_Model_X(data_range,4);
        
        hold on;
        plot(time,select(data_range),'LineStyle','none','MarkerSize',10,'MarkerEdgeColor',color_plots{data},...
             'MarkerFaceColor',color_plots{data},'Marker','o');
        
    end

   legend(legend_plot_S)
    
   xlabel(x_units);
   ylabel(y_units_S{S});
   title(title_plots{S});

   set(findall(gcf,'-property','FontSize','-property', 'FontWeight'),'FontSize',...
              25, 'Fontweight','bold');

end
  
for P = 1:length(product_distribution)
    figure;
   
    
    for data = 1:length(data_range_all)     
        
        product = product_distribution{P};
        data_range = data_range_all{data};
        
        time = All_Model_X(data_range,4);
        
        hold on;
        plot(time,product(data_range),'LineStyle','none','MarkerSize',10,'MarkerEdgeColor',color_plots{data},...
             'MarkerFaceColor',color_plots{data},'Marker','o');
        
    end

   legend(legend_plot)
    
   xlabel(x_units);
   ylabel(y_units{P});
   title(title_plots{P});

   set(findall(gcf,'-property','FontSize','-property', 'FontWeight'),'FontSize',...
              25, 'Fontweight','bold');

end

%% Global optimisation 
bCH4 = mean_multi_runs{1};
bCO  = mean_multi_runs{2};
bH2  = mean_multi_runs{3};
% 
% CH4_rate = @(x) bCH4(1).* (((400.^bCH4(3))./(bCH4(2))).*((x(1)./(bCH4(2))).^((400.^bCH4(3))-1)).*(exp(-1.*((x(1)./(bCH4(2))).^(400.^bCH4(3))))))...
%            .* ((((bCH4(4).*x(2)).^4).*(bCH4(5).*x(3)))./((1 + bCH4(4).*x(2) + bCH4(5).*x(3)).^(5)));
% 
% CO_rate = @(x) bCO(1).* (((400.^bCO(3))./(bCO(2))).*((x(1)./(bCO(2))).^((400.^bCO(3))-1)).*(exp(-1.*((x(1)./(bCO(2))).^(400.^bCO(3))))))...
%            .* ((bCO(4).*x(2).*bCO(5).*x(3))./((1 + bCO(4).*x(2) + bCO(5).*x(3)).^2));
%        
% H2_rate = @(x) bH2(1).* (((400.^bH2(3))./(bH2(2))).*((x(1)./(bH2(2))).^((400.^bH2(3))-1)).*(exp(-1.*((x(1)./(bH2(2))).^(400.^bH2(3))))))...
%            .* ((((bH2(4).*x(2)).^2))./((1 + bH2(4).*x(2) + bH2(5).*x(3)).^2));

CH4_rate = @(x) bCH4(1).* (((400.^bCH4(3))./(bCH4(2))).*((x(:,1)./(bCH4(2))).^((400.^bCH4(3))-1)).*(exp(-1.*((x(:,1)./(bCH4(2))).^(400.^bCH4(3))))))...
           .* ((((bCH4(4).*x(:,2)).^4).*(bCH4(5).*x(:,3)))./((1 + bCH4(4).*x(:,2) + bCH4(5).*x(:,3)).^(5)));

CO_rate = @(x) bCO(1).* (((400.^bCO(3))./(bCO(2))).*((x(:,1)./(bCO(2))).^((400.^bCO(3))-1)).*(exp(-1.*((x(:,1)./(bCO(2))).^(400.^bCO(3))))))...
           .* ((bCO(4).*x(:,2).*bCO(5).*x(:,3))./((1 + bCO(4).*x(:,2) + bCO(5).*x(:,3)).^2));
       
H2_rate = @(x) bH2(1).* (((400.^bH2(3))./(bH2(2))).*((x(:,1)./(bH2(2))).^((400.^bH2(3))-1)).*(exp(-1.*((x(:,1)./(bH2(2))).^(400.^bH2(3))))))...
           .* ((((bH2(4).*x(:,2)).^2))./((1 + bH2(4).*x(:,2) + bH2(5).*x(:,3)).^2));

       
CH4_best = @(x) CH4_rate(x) ./ (CO_rate(x) + H2_rate(x));
CH4_best_opt = @(x) -CH4_best(x);

CO_best  = @(x) CO_rate(x) ./ (CH4_rate(x) + H2_rate(x));
CO_best_opt  = @(x) -CO_best(x);

H2_best  = @(x) H2_rate(x) ./ (CH4_rate(x) + CO_rate(x));
H2_best_opt  = @(x) -H2_best(x);

best_all = {CH4_best_opt, CO_best_opt, H2_best_opt};

lb = [0.5 2.66 25.72];
ub = [5 6.64 98.38];

parpool;

for model = 1:length(best_all) 
tic

max_values{model} = maxIt(best_all{model},length(lb),...
                                10,1,1,lb,ub);
                           
toc
end   

delete(gcp);
