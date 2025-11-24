close all
clear
clc
%% add basic_functions to path
filePath = matlab.desktop.editor.getActiveFilename;
pos=strfind(filePath,filesep);
pathloc=filePath(1:pos(end-1));
addpath([pathloc,'basic_functions'])
makePretty


PathName='D:\data1\martinf\results\Yeranuhi\';
%PathName= 'D:\KLT_illustration\';

Corr_file =  readtable([PathName 'Comparaison_BF_Ca\Correspondance_BF_Ca.xlsx'], 'Sheet', 'MFM');
filenames = Corr_file.File_contrac;
F=length(filenames);
%%
list_param = {'File_Name', 'Strain' 'Num_Particles','Contractile_Particles', 'Beating_Period', 'Contraction_Amplitude','Contraction_time', 'Relaxation_time', 'Contraction_slope', 'Relaxation_slope','Relaxation_slope_00_50', 'Relaxation_slope_50_100', 'Calcium_Period', 'Calcium_Amplitude','Ascend_time', 'Ascend_slope', 'Decay_time', 'Decay_slope'};
Param_contrac_Ca = array2table(zeros(F,length(list_param)));
Param_contrac_Ca.Properties.VariableNames = list_param;

Num_Particles = NaN*ones(F,1);
Contractile_Particles = NaN*ones(F,1);

Beating_Period = NaN*ones(F,1);
Contraction_Amplitude = NaN*ones(F,1);
Contraction_time = NaN*ones(F,1);
Relaxation_time = NaN*ones(F,1);
Contraction_slope = NaN*ones(F,1);
Relaxation_slope = NaN*ones(F,1); 
Relaxation_slope_00_50 = NaN*ones(F,1);
Relaxation_slope_50_100 = NaN*ones(F,1);

Calcium_Period = NaN*ones(F,1);
Calcium_Amplitude = NaN*ones(F,1);
Ascend_time = NaN*ones(F,1);
Ascend_slope = NaN*ones(F,1);
Decay_time = NaN*ones(F,1);
Decay_slope = NaN*ones(F,1);
%%
for f=1:F
    disp([num2str(f) ' | ' num2str(F)]);
    Contrac_FileName = Corr_file.File_contrac{f};
    files = dir([PathName, 'BF\Analysis\MFM\' Contrac_FileName '*PS.mat']);
    load([PathName, 'BF\Analysis\MFM\', files(1).name])

    pos=strfind(Contrac_FileName,'-');

    St = Contrac_FileName(pos(3)+1:pos(4)-1);
    if length(strfind(St, ' '))>0
        pp = strfind(St, ' ');
        Strain{f} = St(1:pp(1)-1);
    else
        Strain{f} = St;
    end
    
    Num_Particles(f) = sum(PS.signalarea);
    Contractile_Particles(f) = sum(PS.signalarea)/PS.NumberOfParticles;


    if sum(PS.signalarea)/PS.NumberOfParticles>0

        Beating_Period(f) = PS.Signal.Period_median;
        Contraction_Amplitude(f) = PS.Signal.Amp_asc_median;
        Contraction_time(f) = PS.Signal.Asc_time_median;
        Relaxation_time(f) = PS.Signal.Decay_time_median;
        Contraction_slope(f) = PS.Signal.Asc_slope_median;
        Relaxation_slope(f) = PS.Signal.Decay_slope_median;
        Relaxation_slope_00_50(f) = PS.Signal.Decay_slope_0_50_median;
        Relaxation_slope_50_100(f) = PS.Signal.Decay_slope_50_100_median;
    end

    Calcium_FileName = Corr_file.File_Ca{f};
    files = dir([PathName, 'Ca\Analysis\MFM\' Calcium_FileName '*PSCa.mat']);
    load([PathName, 'Ca\Analysis\MFM\', files(1).name])
    if length(PSCa.Signal)>0
        Calcium_Period(f) = PSCa.Signal.Period_median;
        Calcium_Amplitude(f) = PSCa.Signal.Amp_asc_median;
        Ascend_time(f) = PSCa.Signal.Asc_time_median;
        Decay_time(f) = PSCa.Signal.Decay_time_median;
        Ascend_slope(f) = PSCa.Signal.Asc_slope_median;
        Decay_slope(f) = PSCa.Signal.Decay_slope_median;
    end

end

Param_contrac_Ca.File_Name = filenames;
Param_contrac_Ca.Strain = Strain';
Param_contrac_Ca.Num_Particles = Num_Particles;
Param_contrac_Ca.Contractile_Particles = Contractile_Particles;
Param_contrac_Ca.Beating_Period = Beating_Period;
Param_contrac_Ca.Contraction_Amplitude = Contraction_Amplitude;
Param_contrac_Ca.Contraction_time = Contraction_time;
Param_contrac_Ca.Relaxation_time = Relaxation_time;
Param_contrac_Ca.Contraction_slope = Contraction_slope;
Param_contrac_Ca.Relaxation_slope = Relaxation_slope;
Param_contrac_Ca.Relaxation_slope_00_50 = Relaxation_slope_00_50;
Param_contrac_Ca.Relaxation_slope_50_100 = Relaxation_slope_50_100;
Param_contrac_Ca.Calcium_Period = Calcium_Period;
Param_contrac_Ca.Calcium_Amplitude = Calcium_Amplitude;
Param_contrac_Ca.Ascend_time = Ascend_time;
Param_contrac_Ca.Decay_time = Decay_time;
Param_contrac_Ca.Ascend_slope = Ascend_slope;
Param_contrac_Ca.Decay_slope = Decay_slope;

writetable(Param_contrac_Ca,[PathName 'Comparaison_BF_Ca\Param_contrac_Ca.xlsx'], 'Sheet','MFM')  ;

%%
Data_CTL = readtable([PathName 'Comparaison_BF_Ca\Param_contrac_Ca.xlsx'], 'Sheet', 'CTL');
Data_MFM = readtable([PathName 'Comparaison_BF_Ca\Param_contrac_Ca.xlsx'], 'Sheet', 'MFM');
Data(1:height(Data_CTL),:) = Data_CTL;
Data(height(Data_CTL)+1:height(Data_CTL)+height(Data_MFM),:) = Data_MFM;
%%
close all

Amplitude = figure;
a_a = axes(Amplitude);
xlabel('Contrac Amplitude')%('$\sigma_{Amp}$')
ylabel('Ca Amplitude')
legend('Location', 'northeast')
hold on


Contraction_time = figure;
a_ct = axes(Contraction_time);
xlabel('Contraction time (s)')
ylabel('Ascend time (s)')
legend('Location', 'northeast')
hold on

Relaxation_time = figure;
a_rt = axes(Relaxation_time);
xlabel('Relaxation time (s)')
ylabel('Decay time (s)')
legend('Location', 'northeast')
hold on

Contraction_slope = figure;
a_cs = axes(Contraction_slope);
xlabel('Contraction speed (um/s)')
ylabel('Ascend speed (s$^{-1}$)')
legend('Location', 'northeast')
hold on

Relaxation_slope = figure;
a_rs = axes(Relaxation_slope);
xlabel('Relaxation speed (um/s)')
ylabel('Decay speed (s$^{-1}$)')
legend('Location', 'northeast')
xticks('auto')
hold on

ind = (Data.Contractile_Particles >0.4);

[rho,pval] = corr(Data.Contraction_Amplitude(ind),Data.Calcium_Amplitude(ind));
scatter(a_a, Data.Contraction_Amplitude(ind), Data.Calcium_Amplitude(ind), 'k','filled','DisplayName','$\rho$ = '+string(rho)+', pval = '+string(pval));

[rho,pval] = corr(Data.Contraction_time(ind),Data.Ascend_time(ind));
scatter(a_ct, Data.Contraction_time(ind), Data.Ascend_time(ind), 'k', 'filled','DisplayName','$\rho$ = '+string(rho)+', pval = '+string(pval));

[rho,pval] = corr(Data.Relaxation_time(ind),Data.Decay_time(ind));
scatter(a_rt, Data.Relaxation_time(ind), Data.Decay_time(ind), 'k', 'filled','DisplayName','$\rho$ = '+string(rho)+', pval = '+string(pval));

[rho,pval] = corr(Data.Contraction_slope(ind),Data.Ascend_slope(ind));
scatter(a_cs, Data.Contraction_slope(ind), Data.Ascend_slope(ind), 'k','filled','DisplayName','$\rho$ = '+string(rho)+', pval = '+string(pval));

[rho,pval] = corr(-Data.Contraction_slope(ind),-Data.Decay_slope(ind));
scatter(a_rs, -Data.Relaxation_slope(ind), -Data.Decay_slope(ind), 'k', 'filled','DisplayName','$\rho$ = '+string(rho)+', pval = '+string(pval));
%%
close all

% Period = figure;
% a_p = axes(Period);
% xticks('auto')
% xlabel('Contrac Period')
% ylabel('Ca Period')%('$\sigma_{Period}$')
% xlim([0 2])
% ylim([0 2])
% legend('Location', 'northeast')

Amplitude = figure;
a_a = axes(Amplitude);

Contraction_time = figure;
a_ct = axes(Contraction_time);

Relaxation_time = figure;
a_rt = axes(Relaxation_time);

Contraction_slope = figure;
a_cs = axes(Contraction_slope);

Relaxation_slope = figure;
a_rs = axes(Relaxation_slope);
xticks('auto')

ind = (Data_CTL.Contractile_Particles >0.4);

[rho,pval] = corr(Data_CTL.Contraction_Amplitude(ind),Data_CTL.Calcium_Amplitude(ind));
scatter(a_a, Data_CTL.Contraction_Amplitude(ind), Data_CTL.Calcium_Amplitude(ind), 'k','filled','DisplayName','CTL, N=' +string(sum(ind))+newline+'$\rho$ = '+string(rho)+', pval = '+string(pval));
xlim(a_a,[0 max(Data.Contraction_Amplitude)])
ylim(a_a,[0 max(Data.Calcium_Amplitude)])
xlabel(a_a,'Contraction Amplitude')%('$\sigma_{Amp}$')
ylabel(a_a,'Ca Amplitude')
legend(a_a,'Location', 'northeast')

[rho,pval] = corr(Data_CTL.Contraction_time(ind),Data_CTL.Ascend_time(ind));
scatter(a_ct, Data_CTL.Contraction_time(ind), Data_CTL.Ascend_time(ind), 'k', 'filled','DisplayName','CTL, N=' +string(sum(ind))+newline+'$\rho$ = '+string(rho)+', pval = '+string(pval));
xlim(a_ct,[0 max(Data.Contraction_time)])
ylim(a_ct,[0 max(Data.Ascend_time)])
xlabel(a_ct,'Contraction time (s)')
ylabel(a_ct,'Ascend time (s)')
legend(a_ct,'Location', 'northeast')

[rho,pval] = corr(Data_CTL.Relaxation_time(ind),Data_CTL.Decay_time(ind));
scatter(a_rt, Data_CTL.Relaxation_time(ind), Data_CTL.Decay_time(ind), 'k', 'filled','DisplayName','CTL, N=' +string(sum(ind))+newline+'$\rho$ = '+string(rho)+', pval = '+string(pval));
xlim(a_rt,[0 max(Data.Relaxation_time)])
ylim(a_rt,[0 max(Data.Decay_time)])
xlabel(a_rt,'Relaxation time (s)')
ylabel(a_rt,'Decay time (s)')
legend(a_rt,'Location', 'northeast')

[rho,pval] = corr(Data_CTL.Contraction_slope(ind),Data_CTL.Ascend_slope(ind));
scatter(a_cs, Data_CTL.Contraction_slope(ind), Data_CTL.Ascend_slope(ind), 'k','filled','DisplayName','CTL, N=' +string(sum(ind))+newline+'$\rho$ = '+string(rho)+', pval = '+string(pval));
xlim(a_cs,[0 max(Data.Contraction_slope)])
ylim(a_cs,[0 max(Data.Ascend_slope)])
xlabel(a_cs,'Contraction speed (um/s)')
ylabel(a_cs,'Ascend speed (s$^{-1}$)')
legend(a_cs,'Location', 'northeast')

[rho,pval] = corr(-Data_CTL.Contraction_slope(ind),-Data_CTL.Decay_slope(ind));
scatter(a_rs, -Data_CTL.Relaxation_slope(ind), -Data_CTL.Decay_slope(ind), 'k', 'filled','DisplayName','CTL, N=' +string(sum(ind))+newline+'$\rho$ = '+string(rho)+', pval = '+string(pval));
xlim(a_rs,[0 max(-Data.Relaxation_slope)])
ylim(a_rs,[0 max(-Data.Decay_slope)])
xlabel(a_rs,'Relaxation speed (um/s)')
ylabel(a_rs,'Decay speed (s$^{-1}$)')
legend(a_rs,'Location', 'northeast')
pause

list_strains = unique(Data_MFM.Strain);
colors = {'r', 'g', 'b', 'm', 'c'};
for s=1:length(list_strains)
    ind = (Data_MFM.Contractile_Particles >0.4) & (ismember(Data_MFM.Strain,list_strains{s}));
    
    [rho,pval] = corr(Data_MFM.Contraction_Amplitude(ind),Data_MFM.Calcium_Amplitude(ind));
    scatter(a_a, Data_MFM.Contraction_Amplitude(ind), Data_MFM.Calcium_Amplitude(ind), colors{s},'filled','DisplayName',[list_strains{s} ', N=' num2str(sum(ind)) newline '$\rho$ = ' num2str(rho) ', pval = ' num2str(pval)]);
    xlim(a_a,[0 max(Data.Contraction_Amplitude)])
    ylim(a_a,[0 max(Data.Calcium_Amplitude)])
    xlabel(a_a,'Contraction Amplitude')%('$\sigma_{Amp}$')
    ylabel(a_a,'Ca Amplitude')
    legend(a_a,'Location', 'northeast')

    [rho,pval] = corr(Data_MFM.Contraction_time(ind),Data_MFM.Ascend_time(ind));
    scatter(a_ct, Data_MFM.Contraction_time(ind), Data_MFM.Ascend_time(ind), colors{s}, 'filled','DisplayName',[list_strains{s} ', N=' num2str(sum(ind)) newline '$\rho$ = ' num2str(rho) ', pval = ' num2str(pval)]);
    xlim(a_ct,[0 max(Data.Contraction_time)])
    ylim(a_ct,[0 max(Data.Ascend_time)])
    xlabel(a_ct,'Contraction time (s)')
    ylabel(a_ct,'Ascend time (s)')
    legend(a_ct,'Location', 'northeast')

    [rho,pval] = corr(Data_MFM.Relaxation_time(ind),Data_MFM.Decay_time(ind));
    scatter(a_rt, Data_MFM.Relaxation_time(ind), Data_MFM.Decay_time(ind), colors{s}, 'filled','DisplayName',[list_strains{s} ', N=' num2str(sum(ind)) newline '$\rho$ = ' num2str(rho) ', pval = ' num2str(pval)]);
    xlim(a_rt,[0 max(Data.Relaxation_time)])
    ylim(a_rt,[0 max(Data.Decay_time)])
    xlabel(a_rt,'Relaxation time (s)')
    ylabel(a_rt,'Decay time (s)')
    legend(a_rt,'Location', 'northeast')

    [rho,pval] = corr(Data_MFM.Contraction_slope(ind),Data_MFM.Ascend_slope(ind));
    scatter(a_cs, Data_MFM.Contraction_slope(ind), Data_MFM.Ascend_slope(ind), colors{s},'filled','DisplayName',[list_strains{s} ', N=' num2str(sum(ind)) newline '$\rho$ = ' num2str(rho) ', pval = ' num2str(pval)]);
    xlim(a_cs,[0 max(Data.Contraction_slope)])
    ylim(a_cs,[0 max(Data.Ascend_slope)])
    xlabel(a_cs,'Contraction speed (um/s)')
    ylabel(a_cs,'Ascend speed (s$^{-1}$)')
    legend(a_cs,'Location', 'northeast')

    [rho,pval] = corr(-Data_MFM.Contraction_slope(ind),-Data_MFM.Decay_slope(ind));
    scatter(a_rs, -Data_MFM.Relaxation_slope(ind), -Data_MFM.Decay_slope(ind), colors{s}, 'filled','DisplayName',[list_strains{s} ', N=' num2str(sum(ind)) newline '$\rho$ = ' num2str(rho) ', pval = ' num2str(pval)]);
    xlim(a_rs,[0 max(-Data.Relaxation_slope)])
    ylim(a_rs,[0 max(-Data.Decay_slope)])
    xlabel(a_rs,'Relaxation speed (um/s)')
    ylabel(a_rs,'Decay speed (s$^{-1}$)')
    legend(a_rs,'Location', 'northeast')

    pause
end
%%
close all

[rho,pval] = corr(Data.Ascend_time,Data.Decay_time);
scatter(Data.Ascend_time, Data.Decay_time, 'k','filled', 'DisplayName', '$\rho$ = '+string(rho)+', pval = '+string(pval));
hold on

Data = readtable([PathName 'Ca\Analysis\MFM_Param_Ca.xlsx'], 'Sheet', 'MFM');
list_strains = unique(Data.Strain);
for s=1:length(list_strains)
ind = ismember(Data.Strain,list_strains{s});
[rho,pval] = corr(Data.Ascend_time(ind),Data.Decay_time(ind));
scatter(Data.Ascend_time(ind), Data.Decay_time(ind), colors{s},'filled', 'DisplayName', '$\rho$ = '+string(rho)+', pval = '+string(pval));
hold on
end

xlabel('Ascend time (s)')
ylabel('Decay Time (s)')

legend()%'CTL', list_strains{1}, list_strains{2}, list_strains{3}, list_strains{4}, list_strains{5})

%%
close all
% 
% Period = figure;
% a_p = axes(Period);
% xticks('auto')
% xlabel('Contrac Period')
% ylabel('Ca Period')%('$\sigma_{Period}$')
% xlim([0 2])
% ylim([0 2])
% hold on

Amplitude = figure;
a_a = axes(Amplitude);
xlabel('Contrac Amplitude')%('$\sigma_{Amp}$')
ylabel('Ca Amplitude')
hold on

% 
% Contraction_time = figure;
% a_ct = axes(Contraction_time);
% xlabel('Contraction time (s)')
% ylabel('Ascend time (s)')
% hold on
% 
% Relaxation_time = figure;
% a_rt = axes(Relaxation_time);
% xlabel('Relaxation time (s)')
% ylabel('Decay time (s)')
% hold on
% 
% Contraction_slope = figure;
% a_cs = axes(Contraction_slope);
% xlabel('Contraction speed (um/s)')
% ylabel('Ascend speed (s$^{-1}$)')
% hold on
% 
% Relaxation_slope = figure;
% a_rs = axes(Relaxation_slope);
% xlabel('Relaxation speed (um/s)')
% ylabel('Decay speed (s$^{-1}$)')
% xticks('auto')
% hold on

Data = readtable([PathName 'Comparaison_BF_Ca\Param_contrac_Ca.xlsx'], 'Sheet', 'CTL');


ind = (Data.Contractile_Particles >0.4);

plot(a_a,Data.Contraction_Amplitude(ind),Data.Calcium_Amplitude(ind),'k*'), hold on


%plot(a_a,c(1) + af(1)*cos(t), c(2) + af(2)*sin(t), 'k')
%scatter(a_p, , , 'k','filled');
% scatter(a_a, Data.Contraction_Amplitude(ind), Data.Calcium_Amplitude(ind), 'k','filled');
% scatter(a_ct, Data.Contraction_time(ind), Data.Ascend_time(ind), 'k', 'filled');
% scatter(a_rt, Data.Contraction_time(ind), Data.Decay_time(ind), 'k', 'filled');
% scatter(a_cs, Data.Contraction_slope(ind), Data.Ascend_slope(ind), 'k','filled');
% scatter(a_rs, Data.Relaxation_slope(ind), Data.Decay_slope(ind), 'k', 'filled');