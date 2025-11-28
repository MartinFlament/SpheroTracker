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
PathName_film='D:\data1\martinf\Films\Yeranuhi\BF\';
%PathName= 'D:\KLT_illustration\';

samples = dir(PathName);
for S=3:length(samples)
%% viewing the amplitude measured with respect to the contrast of the region
    files = dir([PathName, samples(S).name, '\*PS.mat']);
    F = length(files);
    for f=1:F

        im1 = mat2gray(imread([PathName_film, samples(S).name, '\', FileName(1:pos(5)-6), '\', FileName(1:end-7), '.tif'],20));
        im1=im1(2:end,1:end);
        load([PathName, samples(S).name, '\', FileName(1:end-6), 'results.mat'])
        load([PathName, samples(S).name, '\', FileName(1:end-6), 'param.mat'])
        clear Contrast
        for p=1:PS.NumberOfParticles
            im_p = im1(max(1,results(1).pts1(p,1)-param.winSize/2):min(size(im1,1),results(1).pts1(p,1)+param.winSize/2), ...
                max(1,results(1).pts1(p,2)-param.winSize/2):min(size(im1,2),results(1).pts1(p,2)+param.winSize/2));
            contrast = std(im_p);%(max(im_p, [], 'all') - min(im_p, [], 'all'))/(max(im_p, [], 'all')+min(im_p, [], 'all'));
            if PS.signalarea(p)==0
                plot(contrast, max(PS.Amplitude(p,:))-min(PS.Amplitude(p,:)), 'xc')
            else
                plot(contrast, max(PS.Amplitude(p,:))-min(PS.Amplitude(p,:)), 'om')
            end
            hold on
        end
        saveas(gcf,[PathName,FileName(1:end-6),'Resolution.fig'])
        saveas(gcf,[PathName,FileName(1:end-6),'Resolution.png'])
        close
    end
end