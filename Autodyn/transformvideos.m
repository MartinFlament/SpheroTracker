function transformvideos(j, reader, omeMeta, currentfolder, workingfolder, filename)

    iSeries = j;
    a = (iSeries - 1);
    reader.setSeries(a); % set the serie to work with
    stackSizeT = omeMeta.getPixelsSizeT(a).getValue(); % number of frames
    stackSizeX = omeMeta.getPixelsSizeX(0).getValue(); % image width, pixels

cd (workingfolder);

b = num2str(iSeries);
zone_number = sprintf('%02s', b);

finalfilename = [filename '_z' zone_number '.TIF'];
if isfile(finalfilename)
    delete(finalfilename)
end
for i = 1:stackSizeT
    if stackSizeX==2048
        temp_frame_numb = num2str(i);
        frame_number = sprintf('%04s', temp_frame_numb);
        I = bfGetPlane(reader, i);
        I2 = uint8(I/256); % convert to 8 bits
        I3 = imfilter(I2, [0.25 0.25; 0.25 0.25], 'replicate');
        I4 = I3(1:2:end,1:2:end); % I3 + I4 = do a 2x2 binning.
        I5 = imfilter(I4, [0.25 0.25; 0.25 0.25], 'replicate');
        I6 = I5(1:2:end,1:2:end);
        blackimage = sum(I6(:));
            if blackimage > 0
            imwrite(I6, finalfilename,'compression','none','WriteMode','append'); % save the image if the image is not black
            end
        clear I I2 I3 I4 I5 I6 temp_frame_numb frame_number finalfilename blackimage
   else
        temp_frame_numb = num2str(i);
        frame_number = sprintf('%04s', temp_frame_numb);
        finalfilename = [filename '_z' zone_number '_t' frame_number '.TIF'];
        I = bfGetPlane(reader, i);
        I2 = uint8(I/256); % convert to 8 bits
        I3 = imfilter(I2, [0.25 0.25; 0.25 0.25], 'replicate');
        I4 = I3(1:2:end,1:2:end); % I3 + I4 = do a 2x2 binning.
        blackimage = sum(I4(:));
            if blackimage > 0
            imwrite(I4, finalfilename,'compression','none','WriteMode','append'); % save the image if the image is not black
            end
        clear I I2 I3 I4 I5 I6 temp_frame_numb frame_number finalfilename blackimage
    end
end

cd (currentfolder);