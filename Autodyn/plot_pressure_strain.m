close all
clear
clc

% scrsz=get(groot,'ScreenSize');

%% image information

dossier='D:\data1\martinf\results\MFM\';
dossier='D:\data1\martinf\results\CTL\';
PathName='D:\data1\martinf\results\Yeranuhi\BF\Analysis\';
dossier= [PathName 'CTL\'];

dosr=dir([dossier,'*_results.mat']);
dosps=dir([dossier,'*_PS.mat']);
dosm=dir([dossier,'*_mask.mat']);
dosp=dir([dossier,'*_param.mat']);

F = length(dosr);

list_param = {'File_Name','Strain', 'Shear' 'Pressure', 'Enstrophy'};
Param_strain = array2table(zeros(F,length(list_param)));
Param_strain.Properties.VariableNames = list_param;
Shearing = NaN*ones(F,1);
Pressure = NaN*ones(F,1);
Enstrophy = NaN*ones(F,1);
%%
for f=1:F  
    disp([num2str(f) ' | ' num2str(F)]);
    clear Press Shear V
    nom=dosr(f).name;
    nomPS=[dosr(f).name(1:end-12) '_PS.mat'];
    nommask=[dosr(f).name(1:end-12) '_mask.mat'];
    nomparam=[dosr(f).name(1:end-12) '_param.mat'];
    load([dossier,nom])
    load([dossier,nomPS])
    load([dossier,nommask])
    sz=size(mask);
    pas=10;
    [X,Y]=meshgrid(1:pas:sz(2),1:pas:sz(1));
        u_mask = griddata(PS.Center(:,1),PS.Center(:,2),double(PS.signalarea'),X,Y,'linear');
    ind_mask = ~isnan(u_mask);
    ind_signal = logical(u_mask(ind_mask));

    
File_Name{f} = dosr(f).name(1:end-12);


pos=strfind(nom,'-');

    St = nom(pos(3)+1:pos(4)-1);
    if length(strfind(St, ' '))>0
        pp = strfind(St, ' ');
        Strain{f} = St(1:pp(1)-1);
    else
        Strain{f} = St;
    end

% PS.rPosition = sqrt((PS.Center(:,1)-PS.SpheroCenter(1)).^2 + (PS.Center(:,2)-PS.SpheroCenter(2)).^2);
% PS.thPosition=atan2((PS.Center(:,2)-PS.SpheroCenter(2)),(PS.Center(:,1)-PS.SpheroCenter(1)));
% ur=griddata(PS.Center(PS.Signalarea,1),PS.Center(PS.Signalarea,2),PS.rPosition(PS.Signalarea),X,Y,'cubic');
% ut=griddata(PS.Center(PS.Signalarea,1),PS.Center(PS.Signalarea,2),PS.thPosition(PS.Signalarea),X,Y,'cubic');
%     Ur=reshape(ur(~isnan(ur)), sum(~isnan(ur), "all"), 1);
%      Ut=reshape(ut(~isnan(ut)), sum(~isnan(ut), "all"), 1);
 %% calculate pressure
 if sum(PS.signalarea)/PS.NumberOfParticles>0.4
    %for ii=1:length(results)
    
    PS.measureStrain([dossier,nommask], [dossier,nomparam])
    save([dossier,nomPS],'PS','-v7.3','-nocompression')%(1:end-12),'_PS.mat'
    % ux=griddata(PS.Center(:,1),PS.Center(:,2),PS.xPosition(:,ii),X,Y,'cubic');
    % uy=griddata(PS.Center(:,1),PS.Center(:,2),PS.yPosition(:,ii),X,Y,'cubic');
    % 
    %      %% calcul pressure
    %     [epsxx,epsxy]=gradientN_h2(ux,'optimal5',0.1);   
    %     [epsyx,epsyy]=gradientN_h2(uy,'optimal5',0.1);
    % 
    %     epsxy=1/2*(epsxy+epsyx);
    %     epsxx(isnan(epsxx)) = 0;
    %     epsyy(isnan(epsyy)) = 0;
    %     epsxy(isnan(epsxy)) = 0;
    %     [vec, Dt, Pc] = arrayfun(@eigenstrain2Df,epsxx,epsyy,epsxy,'UniformOutput',0); %Pc composante symétrique Dt valeur propre (shear)
    %     P=cell2mat(Pc); 
    %     % Press(:, ii) = reshape(P(~isnan(ut)), sum(~isnan(ut), "all"), 1);
    %     P = P(ind_mask);
    % Press(:, ii) = P(ind_signal);
    % %    P(reshape(ind2,size(ux)))=NaN;
    % 
    %    D=cell2mat(Dt);
    %    D = D(ind_mask);
    %    % Shear(:, ii) = reshape(D(~isnan(ut)), sum(~isnan(ut),"all"), 1);
    %    Shear(:, ii) = D(ind_signal);
    % 
    %    [curlz, cav] = curl(X,Y,ux,uy);
    %    cur = curlz(ind_mask);
    %    V(:,ii)=cur(ind_signal);
    %end
 
    
[pkss, locss] = findpeaks(nanmedian(abs(PS.Strain.Shear),1), 'SortStr','descend');
[pksp, locsp] = findpeaks(nanmedian(abs(PS.Strain.Pressure),1), 'SortStr','descend');
[pksv, locsv] = findpeaks(0.5*nanmean(PS.Strain.Vorticity.^2), 'SortStr','descend');

ux=griddata(PS.Center(:,1),PS.Center(:,2),PS.xPosition(:,locsp(1)),X,Y,'cubic');
    uy=griddata(PS.Center(:,1),PS.Center(:,2),PS.yPosition(:,locsp(1)),X,Y,'cubic');
[epsxx,epsxy]=gradientN_h2(ux,'optimal5',0.1);   
        [epsyx,epsyy]=gradientN_h2(uy,'optimal5',0.1);
    
        epsxy=1/2*(epsxy+epsyx);
        epsxx(isnan(epsxx)) = 0;
        epsyy(isnan(epsyy)) = 0;
        epsxy(isnan(epsxy)) = 0;
        [vec, Dt, Pc] = arrayfun(@eigenstrain2Df,epsxx,epsyy,epsxy,'UniformOutput',0); %Pc composante symétrique Dt valeur propre (shear)
        P=cell2mat(Pc); 
    figure;
    clim = max(abs(P), [], 'all');
scale = (-clim):2*clim/sum(P~=0, 'all'):clim;
map = turbo(length(scale));
imagesc(P.*ind_mask)
colormap(map)
set(gca,'Ydir','reverse')
saveas(gca,[dossier,nom(1:end-12),'_pressure_sign.png'])
close

 Shearing(f) = pkss(1);
 Pressure(f) = pksp(1);
Enstrophy(f)=pksv(1);
%%
close all
Pressfig = figure;
a_p = axes(Pressfig);
xlabel('Time (s)')%('$\sigma_{Amp}$')
ylabel('Pressure')
hold on

Shearfig = figure;
a_s = axes(Shearfig);
xlabel('Time (s)')%('$\sigma_{Amp}$')
ylabel('Shear')
hold on

Vorticity = figure;
a_v = axes(Vorticity);
xlabel('Time (s)')%('$\sigma_{Amp}$')
yyaxis(a_v, 'left')
ylabel('Vorticity')
yyaxis(a_v, 'right')
ylabel('Enstrophy')
hold on

% map=parula(32);
% map2=[parula(16);flip(parula(16))];
% Urn=(Ur-min(Ur(:)))*32/max((Ur(:)-min(Ur(:))));
% Uthn=(Ut-min(Ut(:)))*32/max((Ut(:)-min(Ut(:))));

plot(a_s, nanmedian(abs(PS.Strain.Shear),1));

plot(a_p, nanmedian(abs(PS.Strain.Pressure),1));

yyaxis(a_v, 'left')
plot(a_v, nanmedian(PS.Strain.Vorticity,1));
yyaxis(a_v, 'right')
plot(a_v, 0.5*nanmean(PS.Strain.Vorticity.^2,1));


 saveas(a_s,[dossier,nom(1:end-12),'_shear.png'])
 saveas(a_p,[dossier,nom(1:end-12),'_pressure.png'])
  saveas(a_v,[dossier,nom(1:end-12),'_vorticity.png'])

  close all
 end

end
% for ii=1:100:size(Shear,1)
% plot(a_s, Shear(ii,:),'color',map(max(1,floor(nanmean(Uthn(ii,:)))),:));
% plot(a_p, Press(ii,:),'color',map(max(1,floor(nanmean(Uthn(ii,:)))),:));
% end

%%
Param_strain.File_Name = File_Name';
Param_strain.Strain = Strain';
Param_strain.Shear = Shearing;
Param_strain.Pressure = Pressure;
Param_strain.Enstrophy = Enstrophy;
writetable(Param_strain,[PathName 'Strain_data.xlsx'], 'Sheet','CTL');

%%
Data_CTL = readtable([PathName 'Strain_data.xlsx'], 'Sheet',1);
Data = readtable([PathName 'Strain_data.xlsx'], 'Sheet','MFM');

list_strains = unique(Data.Strain);
list_CTL = unique(Data_CTL.Strain);
colors = {'r', 'g', 'b', 'm', 'c'};

%% Assessing the homogeneity of control
close all

Press = figure;
a_p = axes(Press);
xticks('auto')
ylabel('Pressure')
hold on

Shear = figure;
a_s = axes(Shear);
xticks('auto')
ylabel('Shear')
hold on

Enstrophy = figure;
a_e = axes(Enstrophy);
ylabel('Enstrophy')%('$\sigma_{Amp}$')
xticks('auto')
hold on

for s=1:length(list_CTL)
 ind = ismember(Data_CTL.Strain,list_CTL{s}) & ~isnan(Data_CTL.Shear);
    swarmchart(a_p, (s)*ones(sum(ind)), Data_CTL.Pressure(ind), 'filled');
    swarmchart(a_s, (s)*ones(sum(ind)), Data_CTL.Shear(ind), 'filled');
    swarmchart(a_e, (s)*ones(sum(ind)), Data_CTL.Enstrophy(ind), 'filled');
    a_p.XTickLabel{s} = list_CTL{s};a_s.XTickLabel{s} = list_CTL{s};a_e.XTickLabel{s} = list_CTL{s};
end
a_p.XTick = 1:s;a_s.XTick = 1:s;a_e.XTick = 1:s;
%% Measuring differences between mutations and control
close all

Press = figure;
a_p = axes(Press);
xticks('auto')
ylabel('Pressure')
hold on

Shear = figure;
a_s = axes(Shear);
xticks('auto')
ylabel('Shear')
hold on

Enstrophy = figure;
a_e = axes(Enstrophy);
ylabel('Enstrophy')%('$\sigma_{Amp}$')
xticks('auto')
hold on

ind = ~isnan(Data_CTL.Shear);
    a_p.XTickLabel{1} = 'WT';a_s.XTickLabel{1} = 'WT';a_e.XTickLabel{1} = 'WT';
    swarmchart(a_p, ones(sum(ind)), Data_CTL.Pressure(ind), 'b', 'filled');
    swarmchart(a_s, ones(sum(ind)), Data_CTL.Shear(ind), 'b', 'filled');
    swarmchart(a_e, ones(sum(ind)), Data_CTL.Enstrophy(ind), 'b', 'filled');

str = ["LAYV" "MUJO" "PC173" "PC177" "PC179"];
mut = ["E439K" "P419H" "D214-E245del" "E245D" "S46Y"];
str2mut = dictionary(str, mut);

for s=1:length(list_strains)
 ind = ismember(Data.Strain,list_strains{s}) & ~isnan(Data.Shear);
    swarmchart(a_p, (s+1)*ones(sum(ind)), Data.Pressure(ind), 'filled');
    swarmchart(a_s, (s+1)*ones(sum(ind)), Data.Shear(ind), 'filled');
    swarmchart(a_e, (s+1)*ones(sum(ind)), Data.Enstrophy(ind), 'filled');
    a_p.XTickLabel{s+1} = str2mut(list_strains{s});a_s.XTickLabel{s+1} = str2mut(list_strains{s});a_e.XTickLabel{s+1} = str2mut(list_strains{s});
end
a_p.XTick = 1:s+1;a_s.XTick = 1:s+1;a_e.XTick = 1:s+1;

%% Looking for correlations
figure
ind = ~isnan(Data_CTL.Shear);
scatter(Data_CTL.Shear(ind), Data_CTL.Pressure(ind), 'k','filled','DisplayName','CTL')
hold on
for s=1:length(list_strains)
    ind = ismember(Data.Strain,list_strains{s}) & ~isnan(Data.Shear);
    %[rho,pval] = corr(Data.MM_Amplitude(ind),Data.KLT_Amplitude(ind));
    scatter(Data.Shear(ind), Data.Pressure(ind), colors{s},'filled','DisplayName',list_strains{s})%'$\rho$ = '+string(rho)+', pval = '+string(pval));%list_strains{s}+
end
legend('Location','eastoutside')
xlabel('Shear')
ylabel('Pressure')

figure
ind = ~isnan(Data_CTL.Shear);
scatter(Data_CTL.Shear(ind), Data_CTL.Enstrophy(ind), 'k','filled','DisplayName','CTL')
hold on
for s=1:length(list_strains)
    ind = ismember(Data.Strain,list_strains{s}) & ~isnan(Data.Shear);
    %[rho,pval] = corr(Data.MM_Amplitude(ind),Data.KLT_Amplitude(ind));
    scatter(Data.Shear(ind), Data.Enstrophy(ind), colors{s},'filled','DisplayName',list_strains{s})%'$\rho$ = '+string(rho)+', pval = '+string(pval));%list_strains{s}+
    hold on
end
legend('Location','eastoutside')
xlabel('Shear')
ylabel('Enstrophy')
