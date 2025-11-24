close all
clear
clc
%% add basic_functions to path
filePath = matlab.desktop.editor.getActiveFilename;
pos=strfind(filePath,filesep);
pathloc=filePath(1:pos(end-1));
addpath([pathloc,'basic_functions'])
makePretty


PathName='D:\data1\martinf\results\Yeranuhi\BF\Analysis\';
%PathName= 'D:\KLT_illustration\';

samples = dir(PathName);
  %% Measuring the data
for S=3%:length(samples)
    PixelSize = 1/0.7722; %Indicate here the pixel size in micrometer !
    SamplingFrequency = 1/0.0132; % Indicate here the framerate (in fps) !
    ExcitationFrequency = 1;
    col_line=1; % to get each parameter on  a different column in the excels files if set to 1
    goodaream=1;% To select manually the good areas for measurements  1 manual, 0 auto
    %% 
    
    files = dir([PathName, samples(S).name, '\*results.mat']);
    F = length(files);

    for f=1:F

        tic

        disp([num2str(f) ' | ' num2str(F)]);
        FileName =  files(f).name
        
        results_pathname=[PathName samples(S).name '\' FileName(1:end-12) '.xlsx'];
        mask_pathname=[PathName samples(S).name '\' FileName(1:end-12) '_mask.mat'];
        param_pathname=[PathName samples(S).name '\' FileName(1:end-12) '_param.mat'];

        [Signal,Time,Var]=getsignal([PathName,samples(S).name,'\'],FileName);%pause
        PS = PeriodicSignalsModel(Signal,'PixelSize',PixelSize,'SamplingFrequency',SamplingFrequency,'ExcitationFrequency',ExcitationFrequency,'results_pathname',results_pathname,'col_line',col_line);
        PS.results_pathname = results_pathname;
        view = PeriodicSignalsView(PS,'SaveFolder',[PathName samples(S).name '\'],'SaveFileName',FileName(1:end-12));

        %% 
        
        PS.determineCenter(Var);
        %PS.determineAngle;
        PS.findFourierPeriod(PS.Amplitude, PS.SamplingFrequency);   
        PS.interpolateSignal(Time, PS.Amplitude, 1);   
        PS.findPhaseShift(PS.InterpTime, PS.InterpSignal, mode(PS.FourierPeriod));
        PS.averageSignal;
        
        PS.ExtractSignal('Smoothlength',20,'prop',0.3);
        % view.SeeSignal;
        view.PhaseShift(1);
        view.FourierAmpl(1);
        %view.Trajectories(1);
        

        if sum(PS.signalarea)>0
            disp('Enough points')
            % view.PeaksValleysSignal('sauv', 1)
            % view.DisplacementsatPeaks(1);
        end
        % save([PathName, samples(S).name, '\', FileName(1:end-12), '_PS.mat'],'PS','-v7.3','-nocompression')
        clear PS
        close all
    
    end
end

%% Saving the data in a table
list_param = {'File_Name', 'Strain','Num_Particles','Contractile_Particles', 'Period','Period_std', 'Amplitude', 'Amplitude_decay','Amplitude_std','Contraction_time', 'Relaxation_time', 'Contraction_slope', 'Relaxation_slope','Relaxation_slope_00_50', 'Relaxation_slope_50_100', 'Resting_time', 'SigNoise', 'Homogeneity', 'AUC_median', 'Pressure', 'Shear', 'Enstrophy'};%'Manip', 'Strain', 'Batch', 'Spheroid'
for S=3%length(samples)
    clear Manip Strain Batch Spheroid File
    files = dir([PathName, samples(S).name, '\*PS.mat']);
    F = length(files);

    Param_contrac = array2table(zeros(F,length(list_param)));
    Param_contrac.Properties.VariableNames = list_param;
    Num_Particles = NaN*ones(F,1);
    Contractile_Particles = NaN*ones(F,1);
    Period = NaN*ones(F,1);
    Period_std = NaN*ones(F,1);
    Amplitude = NaN*ones(F,1);
    Amplitude_dec = NaN*ones(F,1);
    Amplitude_std = NaN*ones(F,1);
    Contraction_time = NaN*ones(F,1);
    Relaxation_time = NaN*ones(F,1);
    Resting_time = NaN*ones(F,1);
    Contraction_slope = NaN*ones(F,1);
    Relaxation_slope = NaN*ones(F,1); 
    Relaxation_slope_00_50 = NaN*ones(F,1);
    Relaxation_slope_50_100 = NaN*ones(F,1);
    SigNoise = NaN*ones(F,1);
    Homogeneity = NaN*ones(F,1);
    AUC_median = NaN*ones(F,1);
    Shearing = NaN*ones(F,1);
    Pressure = NaN*ones(F,1);
    Enstrophy = NaN*ones(F,1);

    for f=1:F
        disp([num2str(f) ' | ' num2str(F)]);
        FileName =  files(f).name

        load([PathName, samples(S).name, '\', FileName])

        mask_pathname=[PathName samples(S).name '\' FileName(1:end-7) '_mask.mat'];
        param_pathname=[PathName samples(S).name '\' FileName(1:end-7) '_param.mat'];
        
        pos=strfind(FileName,'-');
        Param_contrac.Properties.RowNames{f} = FileName(1:end-6);

        % Manip{f} = FileName(1:pos(1)-1);
        St = FileName(pos(3)+1:pos(4)-1);
        if length(strfind(St, ' '))>0
            pp = strfind(St, ' ');
            Strain{f} = St(1:pp(1)-1);
        %     Batch{f} = FileName(pos(3)+pp(end)+1:pos(5)-6);
        else
            Strain{f} = St;
        %     Batch{f} = FileName(pos(4)+1:pos(5)-6);
        end
        %     Spheroid{f} = FileName(pos(end-1)+2:end-7);
        % 
        % 
        % 
        if ~ismember(samples(S).name, 'CTL')
            Strain{f} = St(1:min(5,length(St)));
        end

        File{f} = FileName(1:end-7);
        Num_Particles(f) = sum(PS.signalarea);
        Contractile_Particles(f) = sum(PS.signalarea)/PS.NumberOfParticles;


        if sum(PS.signalarea)/PS.NumberOfParticles>0.4
       Period(f) = PS.Signal.Period_median;
            Period_std(f) = PS.Signal.Period_std;
            Amplitude(f) = PS.Signal.Amp_asc_median;
            Amplitude_dec(f) = PS.Signal.Amp_decay_median;
            Amplitude_std(f) = PS.Signal.Amp_asc_std;
            Contraction_time(f) = PS.Signal.Asc_time_median;
            Relaxation_time(f) = PS.Signal.Decay_time_median;
            Resting_time(f) = PS.Signal.Taud_median;
            Contraction_slope(f) = PS.Signal.Asc_slope_median;
            Relaxation_slope(f) = - PS.Signal.Decay_slope_median;
            Relaxation_slope_00_50(f) = - PS.Signal.Decay_slope_0_50_median;
            Relaxation_slope_50_100(f) = - PS.Signal.Decay_slope_50_100_median;
            AUC_median(f) = PS.Signal.AUC_median;

            % [pkss, locss] = findpeaks(nanmedian(abs(PS.Shear),1), 'SortStr','descend');
            % [pksp, locsp] = findpeaks(nanmedian(abs(PS.Pressure),1), 'SortStr','descend');
            % [pksv, locsv] = findpeaks(0.5*mean(PS.Vorticity.^2), 'SortStr','descend');
            % 
            %  Shearing(f) = pkss(1);
            %  Pressure(f) = pksp(1);
            % Enstrophy(f)=pksv(1);

            SigNoise(f) = snr(PS.Signal.FinalSignal);
            Homogeneity(f) = PS.Signal.Homogeneity;
        end
      end

    % Param_contrac.Manip = Manip';
    Param_contrac.Strain = Strain';
    % Param_contrac.Batch = Batch';
    % Param_contrac.Spheroid = Spheroid';
    Param_contrac.File_Name = File';
    Param_contrac.Num_Particles = Num_Particles;
    Param_contrac.Contractile_Particles = Contractile_Particles;
    Param_contrac.Period = Period;
    Param_contrac.Period_std = Period_std;
    Param_contrac.Amplitude = Amplitude;
    Param_contrac.Amplitude_decay = Amplitude_dec;
    Param_contrac.Amplitude_std = Amplitude_std;
    Param_contrac.Contraction_time = Contraction_time;
    Param_contrac.Relaxation_time = Relaxation_time;
    Param_contrac.Resting_time = Resting_time;
    Param_contrac.Contraction_slope = Contraction_slope;
    Param_contrac.Relaxation_slope = Relaxation_slope;
    Param_contrac.Relaxation_slope_00_50 = Relaxation_slope_00_50;
    Param_contrac.Relaxation_slope_50_100 = Relaxation_slope_50_100;
    Param_contrac.SigNoise = SigNoise;
    Param_contrac.Homogeneity = Homogeneity;
    Param_contrac.AUC_median = AUC_median;
    Param_contrac.Shear = Shearing;
    Param_contrac.Pressure = Pressure;
    Param_contrac.Enstrophy = Enstrophy;
    writetable(Param_contrac,[PathName 'MFM_Param_contrac_resphased.xlsx'], 'Sheet',samples(S).name)  ;
end
%%
close all


T_CTL = readtable([PathName 'MFM_Param_contrac_resphased.xlsx'], 'Sheet', 1);
[rho,pval] = corr(T_CTL.Contraction_Amplitude(ind),Data_CTL.Calcium_Amplitude(ind));
scatter(a_a, Data_CTL.Contraction_Amplitude(ind), Data_CTL.Calcium_Amplitude(ind), 'k','filled','DisplayName','CTL, N=' +string(sum(ind))+newline+'$\rho$ = '+string(rho)+', pval = '+string(pval));
xlim(a_a,[0 max(Data.Contraction_Amplitude)])
ylim(a_a,[0 max(Data.Calcium_Amplitude)])
xlabel(a_a,'Contraction Amplitude')%('$\sigma_{Amp}$')
ylabel(a_a,'Ca Amplitude')
legend(a_a,'Location', 'northeast')

