
%% IMPORT Data
% clear
% load('LiveHeLa_0429_20x_MCR_Struct.mat')
close all
s = 0 ;
fg = 1;              %loop count
Power_factor = 1;   %Power measure water and sample
% HeLa0302, Water measured at 50mW; 

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

%Obj_array
Obj = [];
Obj = LiveHeLa_0429_20x;

%update name
tag_num = 20;
tag = [1: tag_num]';
Raw_Str = [];
Raw_Str = "HeLa 0429 20x FOV " + string(tag) + " ";
Raw_Str = convertStringsToChars(Raw_Str);

Raw_Str_short = "HeLa Live 50mW " + string(tag)+ ' ';
Raw_Str_short = convertStringsToChars(Raw_Str_short);

%building Ref_array
Water_Obj = Water_0429_20x;
Ref_array = mean(Water_Obj.bgcorr(4:8,:),1);
Ref_Str = 'Pure Water 50mW ';
factor = factor_Ref4to8_Mix_includeC;  %from workspace
Dry_Str = 'DryMass 7days AD ';

Output = [];

for FOV = 1: tag_num
    
    SpecStr = Raw_Str(FOV,:);
    Power_factor(FOV) = (50 * 80)./(Obj.Power(FOV)*Obj.RunTime(FOV));
    
%  %Dry Mass array
    if isequal(factor(FOV,2),2)
        Dry_cell = DryHeLa_AD_Summary{4,1}.Ave_Nuc';
        DM_Symbol = ' N ';
    end
    if isequal(factor(FOV,2),3)
        Dry_cell = DryHeLa_AD_Summary{4,1}.Ave_Cyto';
        DM_Symbol = ' C ';
    end
    if isequal(factor(FOV,2),1)
        Dry_cell = DryHeLa_AD_Summary{4,1}.Ave_All';
        DM_Symbol = ' A ';
    end    
    
    Raw=[];
    Ref=[];
    Raw = Obj.bgcorr(FOV,:)';
    Ref = Ref_array';
    
    m = 0; n= 0;
    scale_ref = Raw(Index_3400)./Ref(Index_3400);
    Ref_corr = Ref*scale_ref + m;
    Raw_corr = Raw + n;

for rept = 1 : 2

SC_Spec = [];
SC_Spec =  Raw_corr(:,1) - Ref_corr(:,1).*factor(FOV);
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

Cell_water = Raw_corr - (Dry_cell_sm).* scale_dry;
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
    [Raw_Str_short{FOV}, '*',num2str(scaleq)],...
    [Raw_Str_short{FOV}, 'SC'],...
    'Scaled Dry mass', 'Hydration Shell',...
    ['Intracellular Water * ',num2str(scaleq)],...
    'Location','northwest');
    
h.FontSize = 7;
xlabel('Raman shift (cm^{-1})')
ylabel ('Intenisty (a.u.)')

line_x = [3150 3250 3400 3450 3650 peak_3050,peak_2930,peak1_index];
plot_xline(line_x );
plot_xline([2851,2973])
plot_yline(0);

index_min = Index_3040:Index_3080;
Hydr_shell_sm = (Hydr_shell);
Ratio_sm = Hydr_shell_sm(Index_3450)/Hydr_shell_sm(Index_3275);

Ratio_sm = round(Ratio_sm,2);
Min01 = min(Hydr_shell_sm(index_min));
Min02 = Hydr_shell_sm(Index_3150);

Hydr_Ratio = Hydr_shell_sm(Index_3450)/Cell_water(Index_3400,1);
Hydr_Ratio = round(Hydr_Ratio,4);

Bulk_water_CHArea = sum(Ref_corr(Index_2800 : Index_3050));
Bulk_water_OHArea = sum(Ref(Index_2800 : Index_3860))./Power_factor(FOV);    %Same Power Ref

%Paint Hydr Area (Green)
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
OHarea_water1.HandleVisibility = "Off";  %Yellow

OHarea_water2 = area(Wavenumber(Index_3050:Index_3860),...
    Cell_water(Index_3050:Index_3860)*scaleq,...
    'FaceColor','r','FaceAlpha',0.02,'EdgeAlpha',0);
OHarea_water2.HandleVisibility = "Off"; %Pink

Hydr_Ratio_Area = Hydr_Area_OHArea / Cell_water_OHArea;
Hydr_Ratio_Area = round(Hydr_Ratio_Area,4);

Cell_water_conc = Cell_water_OHArea / Bulk_water_OHArea * 55.5;

lowerlimit = min(Hydr_shell(Index_2800: Index_3860));
ylim([lowerlimit*1.5 upperlimit *1.1])
xlim([2400 3900])
title(['#',num2str(FOV),DM_Symbol,num2str(factor(FOV)*100),' %   Ratio ',num2str(Ratio_sm),...
    '  Height Ratio ',num2str(Hydr_Ratio*100),' %',...
    '  Area Ratio ',num2str(Hydr_Ratio_Area *100),' %'])

if mod(fg,2) == 1
    close
else
    
if ismember(FOV,[6 7 9 15])
saveas(gca,['Batch MCR ', SpecStr{1},'.png'])
saveas(gca,['Batch MCR ', SpecStr{1},'.fig'])
close 
end

end

fg = fg +1;

end

%Save the SC
Output.Final_factor(FOV,1) = factor(FOV);
Output.DM_Symbol(FOV,1) = {DM_Symbol};
Output.Ratio(FOV,1) = Ratio_sm;
Output.Height_ratio(FOV,1) = Hydr_Ratio *100;
Output.Area_ratio(FOV,1) = Hydr_Ratio_Area *100;

Output.Cell_water(FOV,:) = Cell_water;
Output.SC_Spec(FOV,:) = SC_Spec;
Output.Hydr_shell(FOV,:) = Hydr_shell;
Output.Cell_water_conc(FOV,1) = Cell_water_conc;
Output.SpecStr(FOV,1) = SpecStr;
Output.OHtoCH3_raw(FOV,1) = Raw_corr(Index_2930);
Output.OHtoCH3_raw(FOV,2) = Raw_corr(Index_3400);
Output.OHtoCH3_raw(FOV,3) = Raw_corr(Index_3400)./Raw_corr(Index_2930);
Output.OHtoCH3_raw(FOV,4) = Hydr_Area_OHArea;
Output.OHtoCH3_raw(FOV,5) = Cell_water_OHArea;


end

%%
Output_to_excel(:,1)= Output.Final_factor;
Output_to_excel(:,2)= Output.Ratio;
Output_to_excel(:,3)= Output.Height_ratio;
Output_to_excel(:,4)= Output.Area_ratio;
Output_to_excel(:,5)= Output.Cell_water_conc;

%% Overlaid Hydration Shell
figure
set(0,'DefaultLineLineWidth',1.0)
for i =1 : size(Output.Final_factor,1)
handle = plot(Wavenumber, Output.Hydr_shell(i,:),'--');
hold on
end

xlabel('Raman shift (cm^{-1})')
ylabel ('Intenisty (a.u.)')
xlim([2600 3850])
ylim([-1.8e04 3.5e04])

line_x = [3150 3250 3275 3400 3450 3650 peak_2930];
plot_xline(line_x)
line_y = [0];
plot_yline(line_y)

%% Calculate the averaged frequency
s= 0;
for FOV = 1:3
SpecStr = Raw_Str(FOV,:);

Index = Index_3200:Index_3800;
Pure_Water = Ref;
Final_SC = Output.Hydr_shell(FOV,:);
Final_SC = Final_SC -Final_SC(Index_3800);
Final_SC = smooth(Final_SC);
Final_SC (Final_SC <0) = 0;

Ave_WN_water = plot_AveFreq(Wavenumber, Index, Pure_Water);
Ave_WN_Hydrshell = plot_AveFreq(Wavenumber, Index, Final_SC);

figure(1)
% ylim([-1.5e04 5.5e05])
xlabel('Raman shift (cm^{-1})')
ylabel ('Intenisty (a.u.)')
legend('Pure water',SpecStr,'Location','northwest')

figure(2)
% ylim([-1.4e-3 7.5e-3])
xlabel('Raman shift (cm^{-1})')
ylabel ('Normalized Intenisty (a.u.)')

delta_WN(FOV,1) = Ave_WN_Hydrshell - Ave_WN_water;
delta_WN(FOV,1)= round(delta_WN(FOV,1), 1);

if s ==1
    figure(1)
    saveas(gca,[SpecStr,'Calculated_Q_',num2str(start_WN),'-1.png'])
    saveas(gca,[SpecStr,'Calculated_Q_',num2str(start_WN),'-1.fig'])
end

end




