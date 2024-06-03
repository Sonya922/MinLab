%% IMPORT Data
clear
load('LiveHeLa_0429_20x_MCR_Struct.mat')
close all
s = 0 ;

%need to update
Obj = [];
Obj = LiveHeLa_0429_20x;
Water_Obj = Water_0429_20x;  
Dry_cell = DryHeLa_AD_Summary{4,1}.Ave_Nuc';  %0329 AD Ave 
Power_factor = 1;

FOV = 2;
Raw_Str = 'LiveHeLa 0429 20x ';
Raw_Str_short = 'LiveHeLa 0429 ';
SpecStr = [Raw_Str, num2str(FOV), ''];
Ref_Str = 'Pure Water Ave4to8 50mW ';
Dry_Str = 'DryMass 0329 7days AD ';
Raw = Obj.bgcorr(FOV,:)';
Ref = Water_Obj.Ave';   
Power = mode(Obj.Power);    % in Power
RunTime = mode(Obj.RunTime)*5;  %in seconds

for i =1:1
figure
hold on
Index_3400 = find(abs(Wavenumber-3400)<1);
scale_ref = Raw(Index_3400)./Ref(Index_3400);

Ref_corr = Ref*scale_ref + 0;
Raw_corr = Raw + 0;

handle = plot(Wavenumber,Ref_corr,'LineWidth',1.0);
handle = plot(Wavenumber,Raw_corr,'LineWidth',1.0);
scale_dry = max(Raw_corr(1065:1075))/max(Dry_cell(1065:1075));   %2930cm-1
handle = plot(Wavenumber,Dry_cell * scale_dry,'LineWidth',1.0);
handle.LineWidth = 1.0;

plot_yline(0);
ylim([-0.2e05 max(Raw_corr)*1.1])
xlabel('Raman shift (cm^{-1})')
ylabel('Intensity (counts)')
legend([Ref_Str,'* ',num2str(scale_ref)],Raw_Str,Dry_Str,'Location','northwest')
title([SpecStr,' ',num2str(Power),' mW for ',num2str(RunTime),' s'])

if s == 1
    cd ..
    cd('Raw')
    saveas(gca,['Result MinusDryLast ',SpecStr,'-01.png'])
    saveas(gca,['Result MinusDryLast ',SpecStr,'-01.fig'])
end

end
flag = true;
SC_Spec = [];

for n =1 
    
peak_2930 = 2936;
peak_3050 = 3061;
peak1_index = 3275;

Index_2930 = find(abs(Wavenumber-peak_2930)<1);
Index_3050 = find(abs(Wavenumber-peak_3050)<1);
Index_3860 = find(abs(Wavenumber-3864)<1);  %end

index_val = [2200,2800,3000,3030,3040,3080,3100,3150,3200,3275,3300,3400,3450,3800]';
convertWNtoIndex(Wavenumber, index_val)

end   %Index Inside

while flag 
factor = double(input('Enter factor: '));

if factor == 1
    break
end

s = 0 ;

for i=1:length(factor)
    
SC_Spec =  Raw_corr(:,1) - Ref_corr(:,1).*factor(i);
scaleq = 0.2;  

%smooth
span = 5;
SC_Spec(1:Index_3030) = smooth(SC_Spec(1:Index_3030), 1, 'moving');
SC_Spec(Index_3030:end) = smooth(SC_Spec(Index_3030:end), span, 'moving');

Dry_cell_sm= Dry_cell;
Dry_cell_sm(1:Index_3030) = smooth(Dry_cell(1:Index_3030), 1, 'moving');
Dry_cell_sm(Index_3030:end) = smooth(Dry_cell(Index_3030:end), span, 'moving');

dry_index = Index_2200: Index_3000;
scale_dry = max(SC_Spec(dry_index))/max(Dry_cell_sm(dry_index));   %2930cm-1
upperlimit = max(max(SC_Spec(dry_index)),max(Raw_corr*scaleq));

Hydr_shell = [];
Hydr_shell = (SC_Spec) - (Dry_cell_sm).* scale_dry;
Hydr_shell = smooth(Hydr_shell, span, 'moving');

Cell_water = Raw_corr - Dry_cell_sm.* scale_dry;
scale_ref = Cell_water(Index_3400)./Ref(Index_3400);
Ref_corr = Ref.* scale_ref;

figure
set(0,'DefaultLineLineWidth',1.5)
hold on
handle = plot(Wavenumber,Ref_corr*scaleq,'--','LineWidth',1.0);
handle = plot(Wavenumber,Raw_corr*scaleq); % Ref scale
handle = plot(Wavenumber,SC_Spec);
handle = plot(Wavenumber,(Dry_cell_sm) * scale_dry);
index= [801:1600];
handle = plot(Wavenumber(index),(Hydr_shell(index)),'-');
handle = plot(Wavenumber(index),(Cell_water(index))*scaleq,'-','LineWidth',1.0);
box on 

h = legend(['Water * ',num2str(scaleq),'* ',num2str(round(scale_ref,2))],...
    [Raw_Str_short, '*',num2str(scaleq)],...
    [Raw_Str_short, 'SC'],...
    'Scaled Dry mass', 'Hydration Shell',...
    ['Intracellular Water * ',num2str(scaleq)],...
    'Location','northwest');
    
h.FontSize = 7;
xlabel('Raman shift (cm^{-1})')
ylabel ('Intenisty (a.u.)')

line_x = [3150 3250 3400 3450 3650 peak_3050,peak_2930,peak1_index];
plot_xline(line_x);
plot_yline(0);

index_min = Index_3040:Index_3080;
Hydr_shell_sm = (Hydr_shell);
Ratio_sm = Hydr_shell_sm(Index_3450)/Hydr_shell_sm(Index_3275);

Ratio_sm = round(Ratio_sm,2);
Min01 = min(Hydr_shell_sm(index_min));
Min02 = Hydr_shell_sm(Index_3150);

Hydr_Ratio = Hydr_shell_sm(Index_3450)/Cell_water(Index_3400,1);
Hydr_Ratio = round(Hydr_Ratio,4);

% index_area_Cell = Index_2930 : Index_3860;
% index_area_Hydr = Index_3100 : Index_3860;
% index_area_Water = Index_2800 : Index_2930;

Bulk_water_CHArea = sum(Ref_corr(Index_2800 : Index_3050));
Bulk_water_OHArea = sum(Ref(Index_2800 : Index_3860))./Power_factor;    %Same Power Ref

%Paint Hydr Area
Hydr_shell_sm(Hydr_shell_sm <0) = 0; %non-negative
Hydr_Area_OHArea = sum(Hydr_shell_sm(Index_3150 : Index_3860));
OHarea = area(Wavenumber(Index_3150:Index_3860),...
    Hydr_shell_sm(Index_3150:Index_3860),...
    'FaceColor','g','FaceAlpha',0.05,'EdgeAlpha',0);
OHarea.HandleVisibility = "Off";

% Cell water
Cell_water_OHArea = sum(Cell_water(Index_3050 : Index_3860))...
    + Bulk_water_CHArea;

OHarea_water1 = area(Wavenumber(Index_2800:Index_3050),...
    Ref_corr(Index_2800:Index_3050)*scaleq,...
    'FaceColor','y','FaceAlpha',0.1,'EdgeAlpha',0);
OHarea_water1.HandleVisibility = "Off";    %Yellow

OHarea_water2 = area(Wavenumber(Index_3050:Index_3860),...
    Cell_water(Index_3050:Index_3860)*scaleq,...
    'FaceColor','r','FaceAlpha',0.02,'EdgeAlpha',0);
OHarea_water2.HandleVisibility = "Off";  %Pink

Hydr_Ratio_Area = Hydr_Area_OHArea / Cell_water_OHArea;
Hydr_Ratio_Area = round(Hydr_Ratio_Area,4);

Cell_water_conc = Cell_water_OHArea / Bulk_water_OHArea * 55.5;

lowerlimit = min(Hydr_shell(Index_2800: Index_3860));
ylim([lowerlimit*1.4 upperlimit *1.1])
xlim([2400 3900])
title(['#',num2str(FOV),' ',num2str(factor(i)*100),' %   Ratio ',num2str(Ratio_sm),...
    '  Height Ratio ',num2str(Hydr_Ratio*100),' %',...
    '  Area Ratio ',num2str(Hydr_Ratio_Area *100),' %'])
end

end

factor_array(FOV,1) = factor;
saveas(gca,['Result MinusDryLast ', SpecStr,'-Contri ',num2str(factor(i)*100),'% ','.png'])

%% Save the SC
if Output == []
    count = 1;
else
    count = size(Output.Final_factor,1);
    count = count + 1;
end

Output.Final_factor(count,1) = factor;
Output.Ratio(count,1) = Ratio_sm;
Output.Height_ratio(count,1) = Hydr_Ratio *100;
Output.Area_ratio(count,1) = Hydr_Ratio_Area *100;

Output.Cell_water(count,:) = Cell_water;
Output.SC_Spec(count,:) = SC_Spec;
Output.Hydr_shell(count,:) = Hydr_shell;
Output.Cell_water_conc(count,1) = Cell_water_conc;

%%
Output_to_excel(:,1)= Output.Final_factor;
Output_to_excel(:,2)= Output.Ratio;
Output_to_excel(:,3)= Output.Height_ratio;
Output_to_excel(:,4)= Output.Area_ratio;
Output_to_excel(:,5)= Output.Cell_water_conc;


%% Gaussian fitting
figure
clear yFitted yFitted1 yFitted2 tbl 

fontSize = 10;
X = Wavenumber(index_gau);
X_origin = Wavenumber(index_gau);
Y = smooth(Residual_SC(index_gau));
tbl = table(X', Y);
% Define the model as Y = a + b*x + c*exp(-(x-d)^2/e) + d * exp(-(x-f)^2/g)
% Note how this "x" of modelfun is related to big X and big Y.
% x(:, 1) is actually X and x(:, 2) is actually Y - the first and second columns of the table.
modelfun = @(b,x) b(1) + b(2) * x(:, 1) + b(3) * exp(-(x(:, 1) - b(4)).^2/b(5)) + b(6) * exp(-(x(:, 1) - b(7)).^2/b(8));  
beta0 = [6, 100, 3275, 3300, 70, 2500, 3450, 50]; % Guess values to start with.  Just make your best guess.
% Now the next line is where the actual model computation is done.
mdl = fitnlm(tbl, modelfun, beta0);
% Now the model creation is done and the coefficients have been determined.
% YAY!!!!

% Extract the coefficient values from the the model object.
% The actual coefficients are in the "Estimate" column of the "Coefficients" table that's part of the mode.
coefficients = mdl.Coefficients{:, 'Estimate'};

% Let's do a fit, but let's get more points on the fit, beyond just the widely spaced training points,
% so that we'll get a much smoother curve.
X = linspace(min(X), max(X), 500); % Let's use 1920 points, which will fit across an HDTV screen about one sample per pixel.
peak1 = round(coefficients(4),0);
peak2 = round(coefficients(7),0);

% Create smoothed/regressed data using the model:
yFitted = coefficients(1) + coefficients(2) * X + coefficients(3) * exp(-(X - coefficients(4)).^2 / coefficients(5)) + ...
	coefficients(6) * exp(-(X - coefficients(7)).^2 / coefficients(8));
% Now we're done and we can plot the smooth model as a red line going through the noisy blue markers.
yFitted1 =  coefficients(3) * exp(-(X - coefficients(4)).^2 / coefficients(5));
yFitted2 =  coefficients(6) * exp(-(X - coefficients(7)).^2 / coefficients(8));

peak1_sum  = sum(yFitted1);
peak2_sum  = sum(yFitted2);
ratio_height = coefficients(6)/coefficients(3);
ratio_height = round(ratio_height,2);
ratio_int = sum(yFitted2)./sum(yFitted1);
ratio_int = round(ratio_int, 2);

hold on;
plot(X_origin,Y,'o','LineWidth', 0.5)
plot(X,yFitted1);
plot(X,yFitted2);
plot(X, yFitted, 'r-', 'LineWidth', 1.5);
grid on;
title(['#',num2str(FOV),' Peak1 = ', num2str(peak1),'  Peak2 = ', num2str(peak2),'  HeightRt = ', num2str(ratio_height),'   IntgRt =', num2str(ratio_int)],'FontSize', fontSize);
xlabel('X', 'FontSize', fontSize);
ylabel('Y', 'FontSize', fontSize);
% legendHandle = legend('Noisy Y', 'Fitted Y', 'Location', 'northeast');
% legendHandle.FontSize = 8;
ylim([-2000 12000])
% Set up figure properties:
% Enlarge figure to full screen.
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% Get rid of tool bar and pulldown menus that are along top of figure.
% set(gcf, 'Toolbar', 'none', 'Menu', 'none');
% Give a name to the title bar.
% set(gcf, 'Name', 'Demo by ImageAnalyst', 'NumberTitle', 'Off') 
box on

saveas(gca,['Result MinusDryLast ', SpecStr,'-fiting ','.png'])
saveas(gca,['Result MinusDryLast ', SpecStr,'-fiting ','.fig'])

