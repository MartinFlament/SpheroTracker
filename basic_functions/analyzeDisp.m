function [Emean, E] = analyzeDisp(cropeddir)

d=dir(strcat(cropeddir,'\*.jpg'));
totcordt = transpose([1;0;0]);
for ix=1:length(d)
  
  fn=d(ix).name
  [pathstr,name,ext] = fileparts(fn)
  fnI = strcat(cropeddir,'\' , name , '.jpg');
  

  rgb = imread(strcat(fnI));
  % figure
  % imshow(rgb)
  gray1 = rgb2gray(rgb);
  % imshow(gray1)
  I = gray1;
  if (ix == 1)
      ref = I;
  end
  
  E(:,:,ix) = disparity(ref, I, 'UniquenessThreshold', 30);
  E2 = filter2(fspecial('average',5),E(:,:,ix))/255;
  E3 = E2 <= 0 | E2 >= 0.55;
  E(:,:,ix) = imcomplement(E3);
  Emean(1,2*ix - 1) = double(mean(mean(E(:,:,ix))))
  Emean(1,2*ix) = double(mean(mean(E(:,:,ix))));


end
[minval, index] = findpeaks(Emean, 'minpeakdistance', 8, 'minpeakheight', .7*(max(Emean)));

plot(Emean); 
% hold on; plot(index, minval, 'k^', 'markerfacecolor', [1 0 0])