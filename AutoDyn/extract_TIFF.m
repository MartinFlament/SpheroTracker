%%
clear all
close all force
clc

%%
currentfolder = cd;
foldername = ('/Users/am3183/Documents/');
reader = bfGetReader('/Users/am3183/MyCore CNRS_Nextcloud/Autobeats/test_ tiff_files/AG08 p56 j30 TA CTL 051022.czi');
omeMeta = reader.getMetadataStore();
Series = reader.getSeriesCount;

filename = '1';
mainfolder = [foldername filename];
%mkdir(mainfolder);

foldername = mainfolder;

if ~strcmp(foldername(end),'/')
    foldername=[foldername '/'];
end

%%

h = waitbar(0, 'processing images');
z = Series;


for j = 1:z
        
    %create appropriate folder
    b = num2str(j);
    zone_number = sprintf('%02s', b);
    workingfoldername = [filename '_z' zone_number];
    workingfolder = [foldername workingfoldername];
    mkdir([workingfolder])
    
    %extract and transform TIF from videos
    tic
    transformvideos(j, reader, omeMeta, currentfolder, workingfolder, filename)
    waitbar(j/z);
    toc
end

close all force
disp('Files transformed !');