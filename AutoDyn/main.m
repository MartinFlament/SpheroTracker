close all force
clear
clc
%% add basic_functions to path
filePath = matlab.desktop.editor.getActiveFilename;
pos=strfind(filePath,filesep);
pathloc=filePath(1:pos(end-1));
addpath([pathloc,'basic_functions'])

% warning('off','MATLAB:str2func:invalidFunctionName');

%% all files
disp('Please choose the folder with your folder movies')
dossier = uigetdir;

dirs = dir(dossier);

disp('Please choose the folder where to save the traces')
savedirr = uigetdir;
savedirr = [savedirr '\'];

%% parameters
prompt = {'\bf \fontsize{12} Do you use the GPU ? (0/1):',...
        '\bf \fontsize{12} Are your movies a .tif stack ? (0/1):'...
    };

    dlgtitle = 'General parameters';
    dims = [1 88];
    definput = {'',''};
    opts.Interpreter = 'tex';
    answers = inputdlg(prompt,dlgtitle,dims,definput, opts);

param.gpu = str2num(answers{1,1});
param.stack = str2num(answers{2,1});
param.minDist = 6;           % minimum distance between the tracked points
param.mask = 1;               % use a mask (1 = yes, 0 = no)
param.pyramidLevel = 3;       % pyramid level used in the KLT program
param.winSize = 1;            % window size used in the KLT program
param.dt = 1;                 % calculate displacement between two frames dt apart
param.method = 0;             % compare image i with image 1 (0) or image i with i+dt (1)
param.removedrift = 1;        % remove drift (1 = yes, 0 = no)
param.maxIteration = 20;      % maximum number of iterations in KLT program
param.threshold = 0.1;        % threshold for convergence in KLT program
param.auto =  0;              % run GUI automatically (1) or manual (0)
param.Ca = 0;

prompt = {'\bf \fontsize{12} Pixel size (µm/px):',...
        '\bf \fontsize{12} Sampling frequency:'...
    };

    dlgtitle = 'General parameters';
    dims = [1 88];
    definput = {'',''};
    opts.Interpreter = 'tex';
    answers = inputdlg(prompt,dlgtitle,dims,definput, opts);

param.PixelSize = str2num(answers{1,1}); %Indicate here the pixel size in micrometer !
param.SamplingFrequency = str2num(answers{2,1}); % Indicate here the framerate (in fps) !
param.col_line=1; % to get each parameter on  a different column in the excels files if set to 1
%% make masks
disp('make masks:');
for dd=1:length(dirs)
g=GUImask(dirs{dd},savedir,param, 0);
waitfor(g); clear g
end

%% choose points parameters for KLT

disp('determine Parameters:');
if param.stack
    for dd=1:length(dirs)
        list = listTiffs(dirs{dd});
        g = GUIparam(dirs{dd}, 'SaveDir', savedir, 'Parameters', param,'num',1);
         waitfor(g)
        
    end
else
    for d = 1:length(dirs)
        disp([num2str(d) ' | ' num2str(length(dirs))]);
        paramf(d) = GUIparam(dirs{d}, 'SaveDir', savedir, 'Parameters', param);
        param=paramf(d);
        
    end
    
end
%% then run klt for all times and all stacks
disp('run KLT:');
if param.stack
    for dd=1:length(dirs)
        list = listTiffs(dirs{dd});
        for d = 1:length(list)
             disp([num2str(d) ' | ' num2str(length(list))]);
             FileName = list{d};
             load([savedir FileName(1:end-4) '_param.mat']);
             paramf.auto=1;
             g = GUIklt_App(dirs{dd}, 'SaveDir', savedir, 'Parameters', paramf, 'num', d);
                waitfor(g); clear g
        end
    end
    
else
    for d = 1:length(dirs)
        disp([num2str(d) ' | ' num2str(length(dirs))]);
        [g, paramf(d)] = GUIklt(dirs{d}, 'SaveDir', savedir, 'Parameters', paramf(d));
        waitfor(g); clear g
        
        
    end
    
end

%% Finally run PeriodicSignalsController for all stacks for which a result was saved
for dd=1:length(dirs)
     g = GUIcontroller_App(dirs{dd}, 'SaveDir', savedir, 'PixelSize', param.PixelSize, 'SamplingFrequency', param.SamplingFrequency);
     waitfor(g);clear g
end

%% Saves Xcel Trace of measured contraction parameters
list_param = {'File_Name', 'Num_Particles','Contractile_Particles', 'Period','Period_std', 'Amplitude', 'Amplitude_decay','Amplitude_std','Contraction_time', 'Relaxation_time', 'Contraction_slope', 'Relaxation_slope','Relaxation_slope_00_50', 'Relaxation_slope_50_100','Homogeneity'};
files = dir([savedir, '*PS.mat']);
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
Contraction_slope = NaN*ones(F,1);
Relaxation_slope = NaN*ones(F,1); 
Relaxation_slope_00_50 = NaN*ones(F,1);
Relaxation_slope_50_100 = NaN*ones(F,1);
Homogeneity = NaN*ones(F,1);

for f=1:F
    disp([num2str(f) ' | ' num2str(F)]);
    FileName =  files(f).name;

    load([files(f).folder, '\', FileName])
    
    pos=strfind(FileName,'-');
    Param_contrac.Properties.RowNames{f} = FileName(1:end-6);

    File{f} = FileName(1:end-7);
    Num_Particles(f) = sum(PS.signalarea);
    Contractile_Particles(f) = sum(PS.signalarea)/PS.NumberOfParticles;


    if sum(PS.signalarea)/PS.NumberOfParticles>0

        Period(f) = PS.Signal.Period_median;
        Period_std(f) = PS.Signal.Period_std;
        Amplitude(f) = PS.Signal.Amp_asc_median;
        Amplitude_dec(f) = PS.Signal.Amp_decay_median;
        Amplitude_std(f) = PS.Signal.Amp_asc_std;
        Contraction_time(f) = PS.Signal.Asc_time_median;
        Relaxation_time(f) = PS.Signal.Decay_time_median;
        Contraction_slope(f) = PS.Signal.Asc_slope_median;
        Relaxation_slope(f) = PS.Signal.Decay_slope_median;
        Relaxation_slope_00_50(f) = -PS.Signal.Decay_slope_0_50_median;
        Relaxation_slope_50_100(f) = -PS.Signal.Decay_slope_50_100_median;
        Homogeneity(f) = 1/PS.Signal.Homogeneity;
    end
  end

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
Param_contrac.Contraction_slope = Contraction_slope;
Param_contrac.Relaxation_slope = Relaxation_slope;
Param_contrac.Relaxation_slope_00_50 = Relaxation_slope_00_50;
Param_contrac.Relaxation_slope_50_100 = Relaxation_slope_50_100;
Param_contrac.Homogeneity = Homogeneity;
writetable(Param_contrac,[savedir 'Contraction_parameters.xlsx'])  ;