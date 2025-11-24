close all force
clear
clc

warning('off','MATLAB:str2func:invalidFunctionName');
%% add basic_functions to path
filePath = matlab.desktop.editor.getActiveFilename;
pos=strfind(filePath,filesep);
pathloc=filePath(1:pos(end-1));
addpath([pathloc,'basic_functions'])
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

prompt = {
        '\bf \fontsize{12} Sampling frequency:'...
    };

    dlgtitle = 'General parameters';
    dims = [1 88];
    definput = {''};
    opts.Interpreter = 'tex';
    answers = inputdlg(prompt,dlgtitle,dims,definput, opts);
param.SamplingFrequency = str2num(answers{1,1}); % Indicate here the framerate (in fps) !
%% all files
dossier = 'D:\data1\martinf\Films\Yeranuhi\Ca\CTL\';
dossier = uigetdir;


dirs = dir(dossier);

savedirr='D:\data1\martinf\results\Yeranuhi\Ca\Analysis\CTL\';
savedirr = uigetdir;
savedirr = [savedirr '\'];
%% Filter images
for d=41%:length(dirs)
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
for dd=3:length(dirs)
disp([num2str(dd) ' | ' num2str(length(dirs))]);
g=GUImask([dossier,'\',dirs(dd).name,'\Data_Ca\'],savedirr,param);
waitfor(g); clear g
end

%%
for dd=4:length(dirs)
    g=GUIcontrollerCa_App([dossier,'\',dirs(dd).name,'\Data_Ca\'], 'SaveDir', savedirr, 'PixelSize', 1, 'SamplingFrequency', param.SamplingFrequency);
    waitfor(g)
end
