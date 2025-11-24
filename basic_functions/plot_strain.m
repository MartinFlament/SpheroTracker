function plot_strain(Xg,Yg,angle,eps,fact,rmicr,map)

U=fact*eps/2.*cos(angle);
V=fact*eps/2.*sin(angle);

if map
pas=1;
[X,Y]=meshgrid(-max(Xg(:)):pas:max(Xg(:)),-max(Yg(:)):pas:max(Yg(:)));

Neps=griddata(Xg,Yg,eps,X,Y);
Neps((X(:).^2+Y(:).^2)<(nanmax(rmicr)).^2)=NaN;
Zg=100*ones(size(Xg));

surf(X,Y,Neps,'edgecolor','none')
shading interp
a=colorbar;
a.Label.Interpreter='Latex';
a.Label.String = '$||Def||$';
quiver3(Xg,Yg,Zg,U,V,zeros(size(Xg)),0,'ShowArrowHead','off','color','k')
quiver3(Xg,Yg,Zg,-U,-V,zeros(size(Xg)),0,'ShowArrowHead','off','color','k')
else
    Zg=100*ones(size(Xg));
    quiver3(Xg,Yg,Zg,U,V,zeros(size(Xg)),0,'ShowArrowHead','off','color','g','LineWidth',2)
quiver3(Xg,Yg,Zg,-U,-V,zeros(size(Xg)),0,'ShowArrowHead','off','color','g','LineWidth',2)
end


