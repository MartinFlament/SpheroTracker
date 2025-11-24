function [dirs,subdirs]=dirbranchendsh(folder)
% find the complete path names of all sub-folders inside the main folder


if ismac||isunix
if ~strcmp(folder(end),'/')
    folder=[folder '/'];
end
else
    if ~strcmp(folder(end),'\')
    folder=[folder '\'];

    end
end
dirs = {folder};
if ismac||isunix
subdirs = {'/'};
else
    subdirs = {''};
end
cont=1;
while cont==1
    D=length(dirs);
    cont=0;
    for d=1:D
  
        dirinfo = dir(dirs{d});
        dirinfo(~[dirinfo.isdir]) = [];  %remove non-directories
       
        dirinfo = dirinfo(~ismember({dirinfo.name},{'.','..'})); %remove current and parent directory
 
  
         if  ~isempty({dirinfo.name})
    
   
             if ismac||isunix
                    
            dirs{d}=strcat(dirs{d},{dirinfo.name},'/');%replace dir with subdirs
            subdirs{d}=strcat(subdirs{d},{dirinfo.name},'/');
            cont=1;
             else
                      
            dirs{d}=strcat(dirs{d},{dirinfo.name},'\');%replace dir with subdirs
            subdirs{d}=strcat(subdirs{d},{dirinfo.name},'\');
    
            
            
            cont=1;
       end
         
        else
     
%             dirs{d}=strcat(dirs{d},{dirinfo.name});
%             subdirs{d}=strcat(subdirs{d},{dirinfo.name});
            dirs{d}=dirs(d);
            subdirs{d}=subdirs(d);       
            %if dir doesnt have subdirs, keep dir
        end
    end
     dirs=[dirs{:}]';
     subdirs=[subdirs{:}]';
    clear dirinfo
end 


