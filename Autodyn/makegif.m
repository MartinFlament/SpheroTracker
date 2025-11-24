 function makegif(f,nomsauv,delaytime,map)
 if nargin==2
     delaytime=0.1;
 end
if nargin<4
    [~,map] =rgb2ind(f(end).cdata,256,'nodither');

end
 
    sz=size(f(1).cdata);
for i=1:length(f)
i
%im(:,:,1,i) = uint8(rgb2ind(f(i).cdata(1:sz(1),1:sz(2),:),map));
if sum(size(f(i).cdata)-sz)~=0
    tata=imresize(f(i).cdata,'OutputSize',sz(1:2));
    titi=uint8(rgb2ind(tata,map));
else
    titi=uint8(rgb2ind(f(i).cdata,map));
end
im(:,:,1,i) = titi;
end
imwrite(im,map,nomsauv,'DelayTime',delaytime,'LoopCount',inf)   
