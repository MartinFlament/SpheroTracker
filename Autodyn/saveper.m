function saveper(PathName,FileName,PS)
   save([PathName FileName(1:end-12) '_period.mat'], 'PS');