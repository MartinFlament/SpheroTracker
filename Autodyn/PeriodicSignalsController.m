close all
clear
clc
%% add basic_functions to path
filePath = matlab.desktop.editor.getActiveFilename;
pos=strfind(filePath,filesep);
pathloc=filePath(1:pos(end-1));
addpath([pathloc,'basic_functions'])
makePretty


PathName='D:\data1\martinf\results\Gabriel\test\';
PathName='D:\data1\martinf\results\Yeranuhi\BF\Analysis\CTL\';

%PathName='D:\data1\martinf\test Fourier\';

files = dir([PathName, 'M039-CT-YH-Ukki026b_2.lif - 10-2*results.mat']);
F = length(files);
    PixelSize = 1/0.7722; %Indicate here the pixel size in micrometer !
    SamplingFrequency = 1/0.0132; % Indicate here the framerate (in fps) !
col_line=1; % to get each parameter on  a different column in the excels files if set to 1
goodaream=1;% To select manually the good areas for measurements  1 manual, 0 auto
  %% 
for f=1:F    
    tic

    disp([num2str(f) ' | ' num2str(F)]);
    FileName = files(f).name;
    results_pathname=[PathName,FileName(1:end-12),'.xlsx'];
    [Signal, Time, Var]=getsignal(PathName,FileName);
    load([PathName,FileName(1:end-12), '_mask.mat'])
    load([PathName,FileName(1:end-12), '_param.mat'])
    %%
    PS = PeriodicSignalsModel(Signal,'PixelSize',PixelSize,'SamplingFrequency',SamplingFrequency,'results_pathname',results_pathname,'col_line',col_line);
        % load([PathName FileName])
    view = PeriodicSignalsView(PS,'SaveFolder',PathName,'SaveFileName',FileName(1:end-12));
    
    
    PS.determineCenter();
    PS.findFourierPeriod(PS.Amplitude, PS.SamplingFrequency);   
    PS.interpolateSignal(Time, PS.Amplitude, 1);   
    PS.findPhaseShift(PS.InterpTime, PS.InterpSignal, mode(PS.FourierPeriod));
    PS.averageSignal;
    PS.ExtractSignal('Smoothlength',20,'prop',0.2);

    if sum(PS.signalarea)>0
                disp('Enough points')

    %        figure;hold on
    % plot(PS.FourierPeriod(PS.signalarea), PS.FourierAmplitude(PS.signalarea), 'o')
    % 
    % pause 
    % close
    % end

    %     plot(f,PS.Signal.Homogeneity, 'o');hold on
    % end
    view.SeeSignal;
    view.PhaseShift(1);

    view.DisplacementsatPeaks(1);
    view.PeaksValleysSignal('sauv', 1);
    
    % parts = find(view.PS.signalarea);
    % N_part = length(parts);
    % for p=1:N_part
    %     ind = (PS.FinalTimeScaled >= PS.Signal.xmc(1)) & (PS.FinalTimeScaled <= PS.Signal.xmr(1));
    %     x_traj = view.PS.xPosition(parts(p),ind);
    %     x_traj = x_traj - x_traj(1);
    %     y_traj = view.PS.yPosition(parts(p),ind);
    %     y_traj = y_traj - y_traj(1);
    %     N=length(x_traj);
    %     ell = fit_ellipse(x_traj, y_traj);
    %     if ~isempty(ell.a)
    %         a_r(p) = min(ell.a,ell.b)/max(ell.a,ell.b);
    %     end
        % figure;plot(x_traj, y_traj);hold on
        % ellipse(ell.a,ell.b, -ell.phi, ell.X0_in, ell.Y0_in, 'r')
        % pause
    % end
    % swarmchart(f*ones(length(a_r)),a_r);hold on
    % aspect_ratio(f) = median(a_r);
    %%
     view.AverageSignal(1);
    
    
    
    %% param cutoff for peaks detection
  
    %PS.Area=[];
    %nareas=2;
    %PS.makeAreas(nareas);  
    %PS.calculateParameters('Smoothlength',10,'prop',0.1);
   
    %% keeping only good areas
    %view.PeaksValleys('sauv',1,'keep',1);
    
    %PS=view.PS;

    %%
    % PS.GetFinalParameters;
    end
    
    %view.SeeArea(1);
    %view.Trajectories
 %pause

      save([PathName,FileName(1:end-12),'_PS.mat'],'PS','-v7.3','-nocompression')%(1:end-12),'_PS.mat'
% 
%      
%     
     clear Time Signal PS a b a_r
     close all
end

