function [Signal,Time, Var]=getsignal(PathName,FileName)
 
    load([PathName FileName]);

     a = logical([results.conver]) ; 
     
        b = any(a == 0, 2);
        
      
       NumberOfParticles = length(b) - sum(b);
    for i = 1 : length(results)
       
        Signal(1 : NumberOfParticles, i, 1) = results(i).pts2(~b, 1);
        Signal(1 : NumberOfParticles, i, 2) = results(i).pts2(~b, 2);

        Time(1 : NumberOfParticles, i) = i;
        %Var(i,:) = max(results(i).var(:,:), [], 'all');
        %h = stretchHistogram(results(i).var(:),0.5);
        Var(i) = mean(abs(results(i).var(:)));
        %figure;imshow(results(i).var);
    end
    end
