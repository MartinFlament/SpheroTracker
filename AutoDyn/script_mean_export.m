clc
clear
close all

%dossier='/data1/thoman/ownCloud/Albano/results/';
dossier='/data1/thoman/ownCloud/Lamia/Results_17/'
dossier='/data1/thoman/Heloise/transfer_7295319_files_9f0a84df/'
files=dir([dossier,'*_PS.mat']);
F=length(files);

col_line=1; % to get each parameter on  a different column in the excels files if set to 1

%% Create table
Tf= array2table(zeros(0,18));
Tf.Properties.VariableNames = {'frequence_bpm', 'Mean_Amp_mum', 'Std_Amp_mum',...
           'Homogeneity_AU','Tau_contract_moy_s','Tau_contract_std_s','Tau_relax_moy_s','Tau_relax_std_s','Tau_Diast_moy_s','Tau_diast_std_s','BAZ_Tau_Diast_moy_s','BAZ_Tau_diast_std_s','AUC_MOY_mums','AUC_std_mums','Pente_contraction_moy_mumpers','pente_contraction_std_mumpers','Pente_relax_moy_mumpers','pente_relax_std_mumpers'};
       
%%
for ii=1:F
    load([dossier,files(ii).name])
%pause
%mean Amp	Amplitude std	??
A=nanmean(PS.FinalSignal);
B=nanstd(PS.FinalSignal);
T=table(60/PS.Period,A,B,PS.Homogeneity,PS.TaucMedMoy, PS.TaucStdMoy,PS.TaurMedMoy,PS.TaurStdMoy,PS.TaudMedMoy,PS.TaudStdMoy,PS.TaudBazMedMoy,PS.TaudBazStdMoy,PS.AireMedMoy, PS.AireStdMoy,PS.pcMedMoy,PS.pcStdMoy,PS.prMedMoy, PS.prStdMoy);
nom=files(ii).name(1:end-7);
nom=erase(nom,' ');
nom=strrep(nom,'.','_');
nom=strrep(nom,'-','_');
nom=strrep(nom,'µ','mu');
T.Properties.RowNames={nom(1:min(length(nom),63))};
T.Properties.VariableNames=Tf.Properties.VariableNames ;
Tf=[Tf;T];
clear T

    
end

save([dossier,'myresult.mat'],'Tf')
if ~col_line
             
            Tfa = table2array(Tf);
            Tff = array2table(Tfa.');
            Tff.Properties.RowNames = Tf.Properties.VariableNames;
            Tff.Properties.VariableNames = Tf.Properties.RowNames;
writetable(Tff,[dossier,'results.xls'],'WriteRowNames',true)  ;
else
writetable(Tf,[dossier,'results.xls'],'WriteRowNames',true)
end

 