for ii=1:length(results)
    mm(ii)=max(results(ii).sp(:,2));
    titi(ii)=results(ii).pts2(100,1);
end

figure
plot(mm)
figure
plot(titi)