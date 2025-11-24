close all
clear
clc
%% add basic_functions to path

PathName='D:\data1\martinf\results\Yeranuhi\BF\';
%PathName= 'D:\KLT_illustration\';

filenames = dir([PathName 'MuscleMotion\CTL\']);
F=length(filenames);
%%
list_param = {'File_Name', 'Strain', 'MM_Amplitude', 'MM_Contraction_time', 'KLT_Amplitude','KLT_Contraction_time'};
Param_contrac_MM = array2table(zeros(F-2,length(list_param)));
Param_contrac_MM.Properties.VariableNames = list_param;
MM_Amplitude = NaN*ones(F-2,1);
MM_Contraction_time = NaN*ones(F-2,1);
KLT_Amplitude = NaN*ones(F-2,1);
KLT_Contraction_time = NaN*ones(F-2,1);

for f=3:F
    disp([num2str(f) ' | ' num2str(F)]);
    FileName = filenames(f).name;
    MM_params = readtable([PathName 'MuscleMotion\CTL\', FileName, '\Overview-results.txt']);
    
    pos=strfind(FileName,'-');

    Strain{f-2} = FileName(pos(3)+1:pos(4)-1);

    PS_files = dir([PathName, 'Analysis\CTL\' FileName(1:pos(end-1)-1) '*PS.mat']);
    if length(PS_files)>0
        load([PathName, 'Analysis\CTL\', PS_files(1).name])
        if sum(PS.signalarea)>0
            plot(median(MM_params.Var8), PS.Signal.Amp_asc_median, '.k')
            File_Name{f-2} = FileName;
            MM_Amplitude(f-2) = median(MM_params.Var8);
            MM_Contraction_time(f-2) = median(MM_params.Var2);
            KLT_Amplitude(f-2) = PS.Signal.Amp_asc_median;
            KLT_Contraction_time(f-2) = PS.Signal.Asc_time_median;
            hold on
        end
    end
end
%%
Param_contrac_MM.File_Name = File_Name';
Param_contrac_MM.Strain = Strain';
Param_contrac_MM.MM_Amplitude = MM_Amplitude;
Param_contrac_MM.MM_Contraction_time = MM_Contraction_time;
Param_contrac_MM.KLT_Amplitude = KLT_Amplitude;
Param_contrac_MM.KLT_Contraction_time = KLT_Contraction_time;
writetable(Param_contrac_MM,[PathName 'Comparaison_MM_KLT.xlsx'])  ;
%%
close all

Amplitude = figure;
a_a = axes(Amplitude);
xlabel('MM Amplitude')%('$\sigma_{Amp}$')
ylabel('KLT Amplitude')
legend('Location', 'eastoutside')
hold on

Data = readtable([PathName 'Comparaison_MM_KLT.xlsx']);


list_strains = unique(Data.Strain);
colors = {'r', 'g', 'b', 'm', 'c'};
for s=1:5%length(list_strains)-1
    ind = ismember(Data.Strain,list_strains{s}) & ~isnan(Data.MM_Amplitude);
    [rho,pval] = corr(Data.MM_Amplitude(ind),Data.KLT_Amplitude(ind));
    scatter(a_a, Data.MM_Amplitude(ind), Data.KLT_Amplitude(ind), colors{s},'filled','DisplayName',list_strains{s})%C;%list_strains{s}+

    
end
%%
close all

Time = figure;
a_ct = axes(Time);
xlabel('MM contraction time (ms)')%('$\sigma_{Amp}$')
ylabel('KLT contraction time (ms)')
legend('Location', 'eastoutside')
hold on

Data = readtable([PathName 'Comparaison_MM_KLT.xlsx']);

SamplingFrequency =0.0132;

list_strains = unique(Data.Strain);
colors = {'k','r', 'g', 'b', 'm', 'c'};
for s=1:length(list_strains)
    ind = ismember(Data.Strain,list_strains{s}) & ~isnan(Data.MM_Amplitude) & (Data.MM_Contraction_time>0);
    [rho,pval] = corr(Data.MM_Contraction_time(ind),1000*Data.KLT_Contraction_time(ind));
    scatter(Data.MM_Contraction_time(ind)+(randn(sum(ind),1)*100*SamplingFrequency), 1000*Data.KLT_Contraction_time(ind)+(randn(sum(ind),1)*100*SamplingFrequency), colors{s},'filled','DisplayName',list_strains{s})%'$\rho$ = '+string(rho)+', pval = '+string(pval))
    xlim([min(Param_contrac_MM.MM_Contraction_time(~isnan(MM_Contraction_time))) max(Param_contrac_MM.MM_Contraction_time(~isnan(MM_Contraction_time)))])
    ylim(1000*[min(Param_contrac_MM.KLT_Contraction_time(~isnan(MM_Contraction_time))) max(Param_contrac_MM.KLT_Contraction_time(~isnan(MM_Contraction_time)))])
end
