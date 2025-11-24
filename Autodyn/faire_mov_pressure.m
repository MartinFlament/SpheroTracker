close all
clear 
clc

% scrsz=get(groot,'ScreenSize');

%% image information

dossier='D:\data1\martinf\results\MFM\';
dossier='D:\data1\martinf\results\CTL\';
dossier='D:\data1\martinf\results\Gabriel\Coeur\';
dossier='D:\data1\martinf\results\Yeranuhi\BF\Analysis\CTL\';
dosr=dir([dossier,'M039-CT-YH-Ukki026b_2.lif - 8-2*_results.mat']);
D = length(dosr);
dosps=dir([dossier,'M039-CT-YH-Ukki026b_2.lif - 8-2*_PS.mat']);
dosm=dir([dossier,'M039-CT-YH-Ukki026b_2.lif - 8-2*_mask.mat']);
dosp=dir([dossier,'M039-CT-YH-Ukki026b_2.lif - 8-2*_param.mat']);

%%
for d=1:D
    clear PS results 
nom=dosr(d).name;
nomPS=dosps(d).name;
nommask=dosm(d).name;
nompam=dosp(d).name;
load([dossier,nom])
load([dossier,nomPS])
load([dossier,nommask])
load([dossier,nompam])
sz=size(mask);
pas=param.minDist;
[X,Y]=meshgrid(1:pas:sz(2),1:pas:sz(1));
 %% calculate pressure
 if PS.Signal.PointNumb/PS.NumberOfParticles>0%.4

for ii=1:length(results)

%for ii=1:10
   


ux=griddata(PS.Center(:,1),PS.Center(:,2),PS.xPosition(:,ii),X,Y,'cubic');
uy=griddata(PS.Center(:,1),PS.Center(:,2),PS.yPosition(:,ii),X,Y,'cubic');

        h=figure;
        fact=10;
        quiver(X,Y,ux*fact,uy*fact,0, 'Color','b');%hold on; imagesc()
        %scatter(PS.Center(:,1)+PS.xPosition(:,ii),PS.Center(:,2)+PS.yPosition(:,ii), 'b.');%hold on; imagesc()
            set(gca,'XTick',[],'YTick',[]);
        axis equal
        f(ii)=getframe(h);
       close

     %% calcul pressure
    [epsxx,epsxy]=gradientN_h2(ux,'optimal5',0.1);   
    [epsyx,epsyy]=gradientN_h2(uy,'optimal5',0.1);

    epsxy=1/2*(epsxy+epsyx);
    epsxx(isnan(epsxx)) = 0;
    epsyy(isnan(epsyy)) = 0;
    epsxy(isnan(epsxy)) = 0;
    [vec, Dt, Pc] = arrayfun(@eigenstrain2Df,epsxx,epsyy,epsxy,'UniformOutput',0); %Pc composante symétrique Dt valeur propre (shear)
    P=cell2mat(Pc); 
%    P(reshape(ind2,size(ux)))=NaN;
   %  h2=figure;
   %  imagesc(P)
   % 
   %   set(gca,'Clim',[-5 5])
   %   set(gca,'Ydir','reverse')
   %  colorbar
   %  axis off
   %  g(ii)=getframe(h2);
   % 
   % 
   % close

   % D=cell2mat(Dt);
   % h2=figure;
   % imagesc(-D)
   % set(gca,'Clim',[-5 5])
   % set(gca,'Ydir','reverse')
   % colorbar
   % axis off
   % m(ii)=getframe(h2);
   % close
%pause
end
 makegif(f,[dossier,[nom(1:end-12),'_beating_movie2.gif']],0.01)

 %makegif(g,[dossier,[nom(1:end-12),'_pressure_movie2.gif']],0.01)

 %makegif(m,[dossier,[nom(1:end-12),'_shear_movie2.gif']],0.01)
% makegif(hh,[dossier,'meanpressure_movie1.gif'],0.01)

% figure
% hold on
% plot(P,'b')
% 
% 
% xlabel('Frame #')
%         saveas(gcf,[dossier,'meanpressure_1_2.png'])
%         saveas(gcf,[dossier,'meanpressure_1_2.fig'])

 end
end