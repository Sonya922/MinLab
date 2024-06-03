
close all
% clear

s= specr;
s.importdata;
Wavenum = s.data.wavenum{1,1};
Wavenumber = Wavenum(1,1:1600);
Intensity= s.data.spc{1,1};
Intensity = Intensity(:,1:1600);
% WN_Water(:,4)= Wavenumber';
% Water(:,4)= Intensity';


Run_num = size(Intensity,1);

i = 1;

while i<= Run_num
    temp_sum = Intensity(i,:);
    Diff = diff(temp_sum);
    %     sort_Diff = sort(Diff,'descend')
    figure(i)
    handle = plot(Wavenumber,temp_sum);
    handle.LineWidth = 1.5;
    hold on
    handle = plot(Wavenumber(2:end),Diff);
    
    line_y = [temp_sum(end),temp_sum(1070),temp_sum(1327)];  %2930 3400 peak
    plot_yline(line_y)
%   handle = yline (temp_sum(847),'--',num2str(temp_sum(847)));  %D2O peak
    answer = input('  Do you want to delete burst? ','s');

    if strcmpi(answer,'b')
        burst_size = 0;
        i = i -1;
        close
    end
    
    if isempty(answer)
        burst_size = 0;
        i= i +1;
        close
    else
        burst_size = str2num(answer);
        for jj=1:burst_size
            
            dcm_obj = datacursormode(gcf);
            set(dcm_obj,'DisplayStyle','datatip','SnapToDataVertex','on','Enable','on')
            disp(['Click on pixel of interest, then press any key',num2str(jj)])
            % Wait while the user does this.
            pause
            
            peak_info = getCursorInfo(dcm_obj);
            burst_peaks(jj)=[peak_info.Position(1)];
            burst_locs(jj)=find(Wavenumber ==burst_peaks(jj));
            
        end
    end
    
    if burst_size > 0
        
        for burst_i = 1: burst_size
            n= 5;
            lower_bound = max(burst_locs(burst_i)-n, 1);
            upper_bound = min(burst_locs(burst_i)+n, 1599);
            
            index = lower_bound:upper_bound+1;
            temp_sum(index) = NaN;
            
%             for index = lower_bound:upper_bound
% %                 if abs(Diff(index))> std(Diff)
%                     temp_sum(index) = NaN;
%                     temp_sum(index+1) = NaN;
% %                 end
%             end
        end
        Intensity(i,:) = fillgaps(temp_sum);
        figure(i)
        handle = plot(Wavenumber,Intensity(i,:));
        handle.LineWidth = 1.5;
        i = i + 1 ;
    end
    
end

%% Plot the trend
Obj = HeLa_DI_20x_0820;
for FOV = 1:1
SpecStr = 'HeLa_DI_20x_0820_';
% Intensity = LiveHeLa_0311_PowerDep_20x.IntData(:,:,FOV);
Intensity = Obj.IntData(:,:,FOV);

index_val = [923,2490,3490,2930,3300,3060,3400, 3864]';
convertWNtoIndex(Wavenumber, index_val)

% range = [1:20,60:80];
range = 1:30;
MinusBG = false;
plotArray = [Index_3400, Index_923, Index_2930, Index_3864];
% plotArray = [Index_2490, Index_2930, Index_3060, Index_3400];  %in D2O PBS
titleArray = {'3400 Signal','Coverslip Signal','2930 Signal','BG Signal'};

figure
for i =1:4
    subplot(2,2,i)
    handle = plot (Intensity(range,plotArray(i)),'LineWidth',1.5);  %3300 peak
    if MinusBG
        handle = plot (Intensity(range,plotArray(i))-Intensity(range,end),'LineWidth',1.5);  %3300 peak
    end
    title(titleArray{i});
end

saveas(gca, ['Intenisty_trend_',SpecStr,num2str(FOV),'.png'])
end

%%
% Cutoff_Index = size(Intensity,1);
Cutoff_Index = 80;
start = 1;
Sum = sum(Intensity(start:Cutoff_Index,:),1);
% 
% Cutoff_Index_range = [1:52,54:60];
% Sum = sum(Intensity(Cutoff_Index_range,:),1);

figure
handle = plot(Wavenumber,Sum);
handle.LineWidth = 1.5;
handle = yline (Sum(end),'--','0');
handle.HandleVisibility = "Off";

s=0;

if s==1
    figure
    handle = plot(Wavenumber,Sum);
    handle.LineWidth = 1.5;
    hold on
    for jj = 1:1
        dcm_obj = datacursormode(gcf);
        set(dcm_obj,'DisplayStyle','datatip','SnapToDataVertex','on','Enable','on')
        disp(['Click on pixel of interest, then press any key',num2str(jj)])
        % Wait while the user does this.
        pause
        peak_info = getCursorInfo(dcm_obj);
        burst_peaks(jj)=[peak_info.Position(1)];
        burst_locs(jj)=find(Wavenumber ==burst_peaks(jj));
        n=5;
        index = burst_locs(jj)-n:burst_locs(jj)+n;
        Sum(index) = NaN;
        Sum = fillgaps(Sum);
    end
    handle = plot(Wavenumber,Sum);
    handle.LineWidth = 1.5;

end

index = 7;
if index == 1
    import_Sum_Data = [];
    import_Intensity_Data = [];
    Cutoff_array = [];
end

import_Sum_Data(index,:) = Sum;
tempIntensity = Intensity(start:Cutoff_Index,:);
row = size(tempIntensity,1);
import_Intensity_Data(1:row,:,index) = tempIntensity;
import_Intensity_Data(row+1:end,:,index)=NaN;
Cutoff_array(index) = Cutoff_Index;


%% after import all Index to the object, run the following
Obj = [];
index_val = [2930]';
convertWNtoIndex(Wavenumber, index_val)
Obj.SumData = import_Sum_Data;
Obj.IntData = import_Intensity_Data;
Obj.bgcorr = import_Sum_Data - import_Sum_Data(:,end);
Obj.MaxNorm = Obj.bgcorr./Obj.bgcorr(:,Index_2930);
Obj.Power = ["50mW"];
Obj.RunTime = ["400s 80Run"];
Water_DI_20x_0820 = Obj;

