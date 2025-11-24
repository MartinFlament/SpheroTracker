close all
clear 
clc
filePath = matlab.desktop.editor.getActiveFilename;
pos=strfind(filePath,filesep);
pathloc=filePath(1:pos(end-1));
addpath([pathloc,'basic_functions'])
% scrsz=get(groot,'ScreenSize');

%% image information

dossier='/data1/thoman/Heloise/ringresults/'
dossierim='/data1/thoman/Heloise/transfer_8013188_files_ce2f8450/'
dosr=dir([dossier,'*_results.mat']);
nom=dosr.name;
dosps=dir([dossier,'*_PS.mat']);
nomPS=dosps.name;
dosm=dir([dossier,'*_mask.mat']);
nommask=dosm.name;

load([dossier,nom])
load([dossier,nomPS])
load([dossier,nommask])
sz=size(mask);
pas=10;
[X,Y]=meshgrid(1:pas:sz(2),1:pas:sz(1));

%% Make a circle to calculate the deformation only in the pillar
files=listTiffs(dossierim);
num=1;
% im=imread([dossierim,files{num}],1);
% figure
% imagesc(im)
% colormap(gray)
% axis equal 
% axis off
% title('Draw a Circular ROI');
% 
% % Use imellipse to draw an ellipse (we will constrain it to be a circle)
% h = imellipse(gca, [50 50 100 100]); % Initial position and size
% api = iptgetapi(h);
% 
% % Constrain the aspect ratio to make it a circle
% api.setFixedAspectRatioMode(true);
% 
% % Wait for the user to double-click to confirm the position
% position = wait(h);
% 
% % Display the position of the circular ROI
% disp('Circular ROI Position:');
% disp(position);
% 
% % Create a mask from the circular ROI
% mask = createMask(h);
% BW=createMask(h);
% 
% save([dossier,files{1}(1:end-4),'_mask_pillar.mat'],'BW')
%%
load([dossier,files{1}(1:end-4),'_mask_pillar.mat']);




 %% calculate pressure
for ii=1:length(results)
    
% for ii=1:10
    ii
   


ux=griddata(PS.Center(:,1),PS.Center(:,2),PS.xPosition(:,ii),X,Y,'cubic');
uy=griddata(PS.Center(:,1),PS.Center(:,2),PS.yPosition(:,ii),X,Y,'cubic');


        h=figure;
        fact=10;
        quiver(X,Y,ux*fact,uy*fact,0)
            set(gca,'ydir','reverse','XTick',[],'YTick',[]);
        axis equal
        f(ii)=getframe(h);
%        close

     %% calcul pressure
    indxy=sub2ind(size(BW),Y(:),X(:));
    ind=BW(indxy)==0;
   ux(ind)=NaN;
   uy(ind)=NaN;
           h=figure;
        fact=10;
        quiver(X,Y,ux*fact,uy*fact,0)
            set(gca,'ydir','reverse','XTick',[],'YTick',[]);
        axis equal

   
  %% Calculate solid rotation u
  indnn=~isnan(ux);
fixedPoints=[X(indnn)+ux(indnn),Y(indnn)+uy(indnn)];
movingPoints=[X(indnn), Y(indnn)];% moving points== points to be transformed
t = cp2tform(movingPoints,fixedPoints,'nonreflective similarity');

u = [0 1]; 
v = [0 0]; 
[x, y] = tformfwd(t,u,v); 
dx = x(2) - x(1)
dy = y(2) - y(1) 
angle = (180/pi) * atan2(dy, dx) 
theta=atan2(dy, dx) 
scale = 1 / sqrt(dx^2 + dy^2)




[Xt,Yt]=tformfwd(t,X(indnn),Y(indnn));
dX=Xt-X(indnn);
dY=Yt-Y(indnn);

% Find the center of rotation
% Assuming we have enough data points
%https://stackoverflow.com/questions/2627073/formula-for-best-approximation-for-center-of-2d-rotation-with-small-angles
%donne le même résultat que cp2tform permet de comprendre
%   vec1=zeros(2*length(X(indnn)),1);
%   vec2=vec1;
%   vec3=vec1;
%   vec4=vec1;
%   B=vec1;
%   vec1(1:2:end)=X(indnn);
%   vec1(2:2:end)=Y(indnn);
%   vec2(1:2:end)=Y(indnn);
%   vec2(2:2:end)=-X(indnn);
%   vec3(1:2:end)=1;
%   vec4(2:2:end)=1;
%   B(1:2:end)=ux(indnn)+X(indnn);
%   B(2:2:end)=uy(indnn)+Y(indnn);
%   A=[vec1 vec2 vec3 vec4];
%   
%   
%   
%   titi=A\B;
% sc=sqrt(titi(1)^2+titi(2)^2)
% dx=titi(3)/sc
% dy=titi(4)/sc
% theta=atan2(titi(2),titi(1))% probleme ici j'obtiens l'inverse de l'angle je ne sais pas pourquoi


% center(1) = 0.5*(dx - sin(theta)*dy/(1-cos(theta)))
% center(2)= 0.5*(dy + sin(theta)*dx/(1-cos(theta)))
% Plot the vector field and the center of rotation
h=figure;
hold on;
quiver(X, Y, ux, uy,0,'b');
quiver(X(indnn),Y(indnn),dX,dY,0,'r')

%  plot(center(1), center(2), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
title('Vector Field with Solid Rotation');
%  set(gca,'ydir','reverse','XTick',[],'YTick',[]);
xlabel('X');
ylabel('Y');
axis equal;
grid on;


   
    [epsxx,epsxy]=gradientN_h2(ux,'optimal5',0.1);   
    [epsyx,epsyy]=gradientN_h2(uy,'optimal5',0.1);

    epsxy=1/2*(epsxy+epsyx);
    
    omega=1/2*(epsxx-epsyy);
    bin=linspace(-1,1,10);
    figure
    histogram(omega,bin)
    Omega= nanmean(omega(:))
    Omega=nanmedian(omega(:))
    V = Omega * (Y(indnn) - center(2));
    U = -Omega * (X(indnn) - center(1));
    figure(h)
    quiver(X(indnn),Y(indnn),U,V,'r')
    
    
    pause
% 
%    
%     epsxx(isnan(epsxx)) = 0;
%     epsyy(isnan(epsyy)) = 0;
%     epsxy(isnan(epsxy)) = 0;
%     [vec, Dt, Pc] = arrayfun(@eigenstrain2Df,epsxx,epsyy,epsxy,'UniformOutput',0);
%     P=cell2mat(Pc);
% %    P(reshape(ind2,size(ux)))=NaN;
%     h2=figure;
%     imagesc(P)
%     
%      set(gca,'Clim',[-5 5])
%     colorbar
%     axis off
%     g(ii)=getframe(h2);
%     
%     h3=figure;
%     plot_strain(X(indo),Y(indo),def_full(tt).th(indo),def_full(tt).def(indo),fact,rmicr,0)
%     
%    close
% %pause
end
 makegif(f,[dossier,'beating_movie2.gif'],0.01)

 makegif(g,[dossier,'pressure_movie2.gif'],0.01)
% makegif(hh,[dossier,'meanpressure_movie1.gif'],0.01)

% figure
% hold on
% plot(P,'b')
% 
% 
% xlabel('Frame #')
%         saveas(gcf,[dossier,'meanpressure_1_2.png'])
%         saveas(gcf,[dossier,'meanpressure_1_2.fig'])

