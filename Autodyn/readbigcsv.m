clc
clear
close all;

%%  add basic_functions to path
filePath = matlab.desktop.editor.getActiveFilename;
pos=strfind(filePath,filesep);
pathloc=filePath(1:pos(end-1));
addpath([pathloc,'basic_functions'])

nom='20240124_10h52m57s_SP+Waveform.csv'
dossier='/data1/thoman/Lamia/'
tic
T=readtable([dossier,nom]);
toc
% 30s for reading, it is ok
%%
figure
plot(T.time_ms,T.ch3_mV)

%% param 
    sm=1;
    prop=0.5; 
    type=2;%2 for MEA;1 for the rest
    param_filter=1000;
    col_line=1; %  set to 1, to get each parameter on  a different column in the excels files
    list_param={'N_pks','Amp_mea','Tau_mea','BI','FPD'};

results_foldername=[dossier,filesep,'Results_',nom(1:end-4), filesep];
if ~exist(results_foldername,'file')
    mkdir(results_foldername);
end

matrix_rough_data=T.Variables;
close all;
clc;
PK=AnalysisPeaks(matrix_rough_data,'prop',prop,'type',type,'param_filter',param_filter,'Smoothness',sm,'list_param_name',list_param,'col_line',col_line);





for i=1:PK.number_cells
    
 PK=PeakAnalysis(PK,i,results_foldername);

end

%%


    results_pathname=[results_foldername filesep 'Results_' file(1:end-5),'.xlsx'];
    save([results_foldername filesep 'Results_' file(1:end-5),'.mat'],'PK');
    PK.Save(results_pathname);

