close all force
clear
clc

warning('off','MATLAB:str2func:invalidFunctionName');
%% add basic_functions to path
filePath = matlab.desktop.editor.getActiveFilename;
pos=strfind(filePath,filesep);
pathloc=filePath(1:pos(end-1));
addpath([pathloc,'basic_functions'])
%% all files
disp('Please choose the folder with your folder movies')
dossier = uigetdir;

dirs = dir(dossier);

disp('Please choose the folder where to save the traces')
savedirr = uigetdir;
savedirr = [savedirr '\'];
%% parameters
prompt = {'\bf \fontsize{12} Are your movies a .tif stack ? (0/1):'...
    };

    dlgtitle = 'General parameters';
    dims = [1 88];
    definput = {''};
    opts.Interpreter = 'tex';
    answers = inputdlg(prompt,dlgtitle,dims,definput, opts);

param.stack = str2num(answers{1,1});
param.minDist =20;           % minimum distance between the tracked points
param.maxFeatures = 10000;      % maximum number of tracked points
param.mask = 1;               % use a mask (1 = yes, 0 = no)
param.pyramidLevel =3;       % pyramid level used in the KLT program
param.winSize = 8;            % window size used in the KLT program
param.dt = 1;                 % calculate displacement between two frames dt apart
param.method = 0;             % compare image i with image 1 (0) or image i with i+dt (1)
param.removedrift = 0;        % remove drift (1 = yes, 0 = no)
param.maxIteration = 20;      % maximum number of iterations in KLT program
param.threshold = 0.1;        % threshold for convergence in KLT program
param.auto = 1;              % run GUI automatically (1) or manual (0)
param.Ca=1;                   % put to 1 if you are doing mask on Calcium data
param.fsize=21;             %size for stdfilter

prompt = {'\bf \fontsize{12} Sampling frequency:'...
    };

    dlgtitle = 'General parameters';
    dims = [1 88];
    definput = {''};
    opts.Interpreter = 'tex';
    answers = inputdlg(prompt,dlgtitle,dims,definput, opts);
param.SamplingFrequency = str2num(answers{1,1}); % Indicate here the framerate (in fps) !
%% Filter images
for d=1:length(dirs)
    disp([num2str(d-3) ' | ' num2str(length(dirs)-3)])
    dossierf=[dossier,'\',dirs(d).name,'\Data_Ca\'];
    if exist(dossierf, 'dir')
        files = dir([dossierf, '*.tif']);
        for k = 1 : length(files)
          baseFileName = files(k).name;
          delete([dossierf files(k).name]);
        end
        rmdir(dossierf)
    end
    mkdir(dossierf)

    files=listTiffs([dossier,'\',dirs(d).name,filesep]);
    for ff=1:length(files)
        disp([num2str(ff) ' | ' num2str(length(files))]);
        clear st
        info=imfinfo([dossier,'\',dirs(d).name,filesep,files{ff}]);
        %read first 100 images to filter
        for ii=1:length(info)
            st(:,:,ii)=double(imread([dossier,'\',dirs(d).name,filesep,files{ff}],ii));
        end
        
        stf=Kalman_Stack_Filter(st,0.8);
        stf=stf-min(stf(:));
        stf=stf/max(stf(:))*256;
        
        for ii=1:length(info)
        imwrite(uint8(stf(:,:,ii)),[dossierf,files{ff}(1:end-4),'_f.tif'],'compression','none','writemode','append')
        end
    
    end
end

%% make masks
for dd=1:length(dirs)
disp([num2str(dd) ' | ' num2str(length(dirs))]);
g=GUImask([dossier,'\',dirs(dd).name,'\Data_Ca\'],savedirr,param);
waitfor(g); clear g
end

%%
for dd=1:length(dirs)
    g=GUIcontrollerCa_App([dossier,'\',dirs(dd).name,'\Data_Ca\'], 'SaveDir', savedirr, 'PixelSize', 1, 'SamplingFrequency', param.SamplingFrequency);
    waitfor(g)
end
%% Saves Xcel Trace of measured calcium parameters
list_param = {'File_Name', 'Period','Period_std', 'Amplitude','Amplitude_std','Ascend_time','Decay_time', 'Homogeneity'};
files = dir([savedir, '*PS.mat']);
F = length(files);

Param_contrac = array2table(zeros(F,length(list_param)));
Param_contrac.Properties.VariableNames = list_param;
Num_Particles = NaN*ones(F,1);
Contractile_Particles = NaN*ones(F,1);
Period = NaN*ones(F,1);
Period_std = NaN*ones(F,1);
Amplitude = NaN*ones(F,1);
Amplitude_std = NaN*ones(F,1);
Ascend_time = NaN*ones(F,1);
Decay_time = NaN*ones(F,1);
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
        Amplitude_std(f) = PS.Signal.Amp_asc_std;
        Ascend_time(f) = PS.Signal.Asc_time_median;
        Decay_time(f) = PS.Signal.Decay_time_median;
        Homogeneity(f) = PS.Signal.Homogeneity;
    end
  end

Param_contrac.File_Name = File';
Param_contrac.Num_Particles = Num_Particles;
Param_contrac.Contractile_Particles = Contractile_Particles;
Param_contrac.Period = Period;
Param_contrac.Period_std = Period_std;
Param_contrac.Amplitude = Amplitude;
Param_contrac.Amplitude_std = Amplitude_std;
Param_contrac.Ascend_time = Ascend_time;
Param_contrac.Decay_time = Decay_time;
Param_contrac.Homogeneity = Homogeneity;
writetable(Param_contrac,[savedir 'Ca_parameters.xlsx'])  ;