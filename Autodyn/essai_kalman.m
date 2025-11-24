clear;
close all;
clc;
%% add basic_functions to path
filePath = matlab.desktop.editor.getActiveFilename;
pos=strfind(filePath,filesep);
pathloc=filePath(1:pos(end-1));
addpath([pathloc,'basic_functions'])
%%
dossier = uigetdir('/data1/thoman/Martin/Data_Ca');
file=listTiffs(dossier);
info=imfinfo([dossier,filesep,file{1}]);
%read first 100 images to filter
for ii=1:length(info)
    st(:,:,ii)=double(imread([dossier,filesep,file{1}],ii));
end

stf=Kalman_Stack_Filter(st,0.8);
%%
figure
imagesc(stf(:,:,end))
colormap(gray)

figure
imagesc(st(:,:,end))
colormap(gray)

counts = imhist(stf(:,:,end)/256, 256);
th=otsuthresh(counts);
BW=imbinarize(stf(:,:,end)/256,th);
figure
imagesc(BW)
%%
imfilt =im2uint8(mat2gray( stdfilt(stf(:,:,end),ones(21,21))));
%MA.Mask = im2bw(imfilt,MA.Th/256);
figure
imagesc(imfilt)
counts = imhist(imfilt/256, 256);
th=otsuthresh(counts);
BW=imbinarize(imfilt/256,th);
figure
imagesc(BW)
%%
        [B,L] = bwboundaries(BW,'noholes');

        boundary=[];
        for k = 1:length(B)
            boundary = cat(1,boundary,B{k});
        
        end
[z,R]= fitcircle(boundary);
%%
figure
hold on
imagesc(stf(:,:,end))
colormap(gray)
            viscircles([z(2) z(1)],R)
            set(gca,'Ydir','reverse')
        axis equal
        axis off