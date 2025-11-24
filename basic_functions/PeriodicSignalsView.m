classdef PeriodicSignalsView < handle
    properties (SetAccess = private)
        PS
        SaveFolder
        SaveFileName
        couleur={'y','m','c','r','g','b','w','k'};
    end
    methods
        function view = PeriodicSignalsView(varargin)
            
            p = inputParser;
            addRequired(p, 'PS', @(x) isa(x,'PeriodicSignalsModel'));
            addOptional(p, 'SaveFolder', @(x) isa(x,'String'));
            addOptional(p, 'SaveFileName', @(x) isa(x,'String'));
            %addParameter(p, 'Save', 'No', @(x) any(validatestring(x, {'Yes','No'})));
            parse(p, varargin{:});
            
            view.PS = p.Results.PS;
            view.SaveFolder = p.Results.SaveFolder;
            view.SaveFileName = p.Results.SaveFileName;
            %view.Save = p.Results.Save;
        end
        function Position(view)
            figure;
            for p = 2 : round(view.PS.NumberOfParticles/10) : view.PS.NumberOfParticles
                subplot(1, 2, 1); plot(view.PS.xPosition(p, :)); hold on; xlabel('FrameNumber'); ylabel('$x - x_0$');
                subplot(1, 2, 2); plot(view.PS.yPosition(p, :)); hold on; xlabel('FrameNumber'); ylabel('$y - y_0$');
            end
            if ~isempty(view.PS.Period)
                subplot(1, 2, 1); axis([0 5 * view.PS.Period min(view.PS.xPosition(:)) max(view.PS.xPosition(:))]);
                subplot(1, 2, 2); axis([0 5 * view.PS.Period min(view.PS.yPosition(:)) max(view.PS.yPosition(:))]);
            end
            drawnow; set(get(handle(gcf),'JavaFrame'),'Maximized',1);
            
            figure; hold on
            for p = 2 : round(view.PS.NumberOfParticles/10) : view.PS.NumberOfParticles
                plot(view.PS.Amplitude(p, :));
            end
            if ~isempty(view.PS.Period)
                axis([0 5 * view.PS.Period min(view.PS.Amplitude(:)) max(view.PS.Amplitude(:))]);
            end
            xlabel('FrameNumber'); ylabel('$|P|$');
            drawnow; set(get(handle(gcf),'JavaFrame'),'Maximized',1);
        end
        function FourierPeriod(view)
            figure; hold on
            plot(view.PS.FourierPeriod, '.');
            xlabel('Particle number'); ylabel('Fourier Period');
            drawnow; set(get(handle(gcf),'JavaFrame'),'Maximized',1);
        end
        
        
        
        function PhaseShiftold(view)
            %             figure;
            %             for p = 1 : round(view.PS.NumberOfParticles/10) : view.PS.NumberOfParticles
            %                 subplot(1, 2, 1); hold on; plot(view.PS.InterpTime, view.PS.InterpSignal(p, :)); xlabel('$t$'); ylabel('$|P|$');
            %                 subplot(1, 2, 2); hold on; plot(view.PS.ShiftTime, view.PS.ShiftSignal(p,:));  xlabel('$t - t_{start}$'); ylabel('$|P|$');
            %             end
            %             if ~isempty(view.PS.Period)
            %                 subplot(1, 2, 1); axis([0 1 * view.PS.Period min(view.PS.InterpSignal(:)) max(view.PS.InterpSignal(:))]);
            %                 subplot(1, 2, 2); axis([0 1 * view.PS.Period min(view.PS.ShiftSignal(:)) max(view.PS.ShiftSignal(:))]);
            %             end
            %             drawnow;
            %             set(get(handle(gcf),'JavaFrame'),'Maximized',1);
            
            figure;
            N = 15;
            
            P = zeros(N);
            compt=zeros(N);
            for p = 1 : view.PS.NumberOfParticles
                x = round(view.PS.Center(p, 1)/max(view.PS.Center(:, 1)) * (N - 1)) + 1;
                y = round(view.PS.Center(p, 2)/max(view.PS.Center(:, 2)) * (N - 1)) + 1;
                P(x, y) = P(x,y)+view.PS.PhaseShift(p);
                compt(x,y)=compt(x,y)+1;
            end
            P=P./compt;
            % gevaarlijk
            m = nanmean(P(:)); s = nanstd(P(:));
            P(P < (m - 4*s) | P > (m + 4*s)) = NaN;
            %gevaarlijk
            h=imagesc(P);
            %             h = imshow(P, [], 'InitialMagnification', 32000/N);
            %             alpha = ones(size(P));
            %             alpha(isnan(P)) = 0;
            %             set(h, 'AlphaData', alpha);
            colorbar
            set(gca,'clim',[-pi pi])
            colormap(gca, 'parula');
            title('PhaseShift');

            
        end
        function Angle(view)
            figure;
            N = 10;
            P = zeros(N);
            compt=zeros(N);
            for p = 1 : view.PS.NumberOfParticles
                x = round(view.PS.Center(p, 1)/max(view.PS.Center(:, 1)) * (N - 1)) + 1;
                y = round(view.PS.Center(p, 2)/max(view.PS.Center(:, 2)) * (N - 1)) + 1;
                P(x, y) = P(x,y)+view.PS.Angle(p);
                compt(x,y)=compt(x,y)+1;
            end
            P=P./compt;
            %           h = imshow(P, [], 'InitialMagnification', 32000/N);
            h = imagesc(P);
            %             alpha = ones(size(P));
            %             alpha(isnan(P)) = 0;
            %             set(h, 'AlphaData', alpha);
            caxis([-pi pi])
            %             %cmap = colorcet('C2');
            %             cmap = colorcet('C4');
            %             colormap(gca, cmap);
            colorbar
            
            title('Angle');
            
        end
        
        function AverageSignal(view,sauv)
            figure; hold on
            plot(view.PS.FinalTimeScaled,view.PS.FinalSignalScaled)
            xlabel('time')
            ylabel('displacement')
            drawnow;
            set(gcf,'units','normalized','outerposition',[0 0 1 1]);
            if sauv
                saveas(gcf,[view.SaveFolder,view.SaveFileName,'_Av_signal.png'])
            end
            close
        end
        
        function SeeSignal(view)
            MX=max(view.PS.Center(:,1));
            MY=max(view.PS.Center(:,2));
            figure;hold on
            scatter(view.PS.Center(view.PS.signalarea,2),view.PS.Center(view.PS.signalarea,1),100,ColorSpecToRGB('m'),'filled')
            scatter(view.PS.Center(~view.PS.signalarea,2),view.PS.Center(~view.PS.signalarea,1),100,ColorSpecToRGB('c'),'filled')
            %   scatter(view.PS.Center(:,2),view.PS.Center(:,1),100,tabcol,'filled')
            
            set(gca,'Xdir','reverse','Xlim',[1, MX+10],'Ylim',[1,MY+10])
            
            axis square
            axis off
            
            saveas(gcf,[view.SaveFolder,view.SaveFileName,'_signal_areas.png'])
            close
        end

        function SeeArea(view,sauv)
            
            ind=~isnan(view.PS.T(:));

            for ii=1:length(view.PS.T(:))
               
                if ~isnan(view.PS.T(ii))
                    tabcol(ii,:)=ColorSpecToRGB(view.couleur{view.PS.T(ii)});
                else
                    tabcol(ii,:)=NaN;
                end
            end
         
            MX=max(view.PS.Center(1:view.PS.TotalPoints-1,1));
            MY=max(view.PS.Center(1:view.PS.TotalPoints-1,2));
            figure
            scatter(view.PS.Center(ind,1),view.PS.Center(ind,2),100,tabcol(ind,:),'filled')
            %   scatter(view.PS.Center(:,2),view.PS.Center(:,1),100,tabcol,'filled')
            
            set(gca,'Xdir','reverse','Xlim',[1 MX+10],'Ylim',[1,MY+10])
            
            axis square
            axis off
            
            if sauv
                saveas(gcf,[view.SaveFolder,view.SaveFileName,'_diff_areas.png'])
            end
            
        end
        
        function SlidingFFT(view,save)
            if(isfield(view.PS.SFFT,'box'))
                
                box=view.PS.SFFT.box;
                Ft=view.PS.FinalTime(box/2+1:(length(view.PS.FinalSignal)-box/2));
                length(view.PS.FinalTime)
                Fs=view.PS.FinalSignal(box/2+1:(length(view.PS.FinalSignal)-box/2));
                as=view.PS.SFFT.as;
             
                rate=view.PS.SFFT.rate;
                
                figure
                hold on
                plot(view.PS.FinalTime,view.PS.FinalSignal,'g')
                plot(Ft,Fs,'b')
                plot(Ft(as),Fs(as),'+r')
                xlim([0 max(view.PS.FinalTime)])
                ylim([0 max(view.PS.FinalSignal)])
                title(['asynchrony rate : ',num2str(rate)])
            else
                for kk=1:length(view.PS.Area)
                    box=view.PS.Area(kk).SFFT.box;
                    ind=isnan(view.PS.Area(kk).FinalSignaln);
                    Time=view.PS.FinalTime(~ind);
                    Ft=Time(box/2+1:(length(view.PS.Area(kk).FinalSignal)-box/2));
                    Fs=view.PS.Area(kk).FinalSignal(box/2+1:(length(view.PS.Area(kk).FinalSignal)-box/2));
                    as=view.PS.Area(kk).SFFT.as;
                    rate=view.PS.Area(kk).SFFT.rate;
                    
                    figure
                    hold on
                    plot(view.PS.FinalTime,view.PS.Area(kk).FinalSignaln,'g')
                    plot(Ft,Fs,'b')
                    plot(Ft(as),Fs(as),'+r')
                    title(['asynchrony rate for area ',num2str(kk),' : ',num2str(rate)])
                    xlim([0 max(view.PS.FinalTime)])
                    ylim([0 max(view.PS.Area(kk).FinalSignaln)])
                    if nargin==2
                        if save
                            hgexport(gcf,[view.SaveFolder view.SaveFileName '_Slidingfft_',num2str(kk),'.png'],hgexport('factorystyle'),'Format','png');
                            
                        end
                    end
                    
                end
            end
        end
        
        function PhaseShift(view,sauv)
                    pas=10;
        [X,Y]=meshgrid(1:pas:max(view.PS.Center(:,1)),1:pas:max(view.PS.Center(:,2)));
        
        FP=scatteredInterpolant(view.PS.Center(view.PS.signalarea,2),view.PS.Center(view.PS.signalarea,1),view.PS.PhaseShift(view.PS.signalarea)');
        FP.Method='nearest';
        FP.ExtrapolationMethod='none';
    
        PHASE=FP(X,Y);
        figure;
        imagesc(PHASE)
        colorbar
        axis equal
        axis off
                  if sauv
                    saveas(gcf,[view.SaveFolder,view.SaveFileName,'_Phase_Shift.png'])
                  end
            close
        end

        function FourierAmpl(view,sauv)
                    pas=10;
        [X,Y]=meshgrid(1:pas:max(view.PS.Center(:,1)),1:pas:max(view.PS.Center(:,2)));
        
        FP=scatteredInterpolant(view.PS.Center(view.PS.signalarea,2),view.PS.Center(view.PS.signalarea,1),view.PS.FourierAmplitude(view.PS.signalarea)');
        FP.Method='nearest';
        FP.ExtrapolationMethod='none';
    
        AMP=FP(X,Y);
        figure;
        imagesc(AMP)
        colorbar
        axis equal
        axis off
                  if sauv
                    saveas(gcf,[view.SaveFolder,view.SaveFileName,'_Fourier_Ampl.png'])
                  end
            close
        end

        function SlidingFFTth(view)
            
            if isfield(view.PS.SFFT,'val')
                th=0.75*(nanmedian(view.PS.SFFT.val(:)));
                val=view.PS.SFFT.val;
                box=view.PS.SFFT.box;
                Ft=view.PS.FinalTime(box/2+1:(length(view.PS.FinalSignal)-box/2));
                figure
                hold on
                plot(Ft,val)
                plot([Ft(1) Ft(end)],[th th],'r')
            
                
            else
                for kk=1:length(view.PS.Area)
                    
                    th=0.75*(nanmedian(view.PS.Area(kk).SFFT.val(:)));
                    val=view.PS.Area(kk).SFFT.val;
                               box=view.PS.Area(kk).SFFT.box;
              Ft=Time(box/2+1:(length(view.PS.Area(kk).FinalSignal)-box/2));
                    figure
                    hold on
                    plot(Ft,val)
                    plot([Ft(1) Ft(end)],[th th],'r')
                    
                    title(['th for area ',num2str(kk)])
                end
            end
        end
        
        function PeaksValleysSignal(view,varargin)
                p = inputParser;

            addOptional(p, 'sauv',1);
            parse(p, varargin{:});
            sauv=p.Results.sauv;
                
                pos=view.PS.Signal.posper(:,:,1);
                M=view.PS.Signal.Mper(:,:,1);
                
               figure
                hold on
               plot(view.PS.Signal.FinalTimeScaled,view.PS.Signal.FinalSignalScaled)
               % plot(view.PS.FinalTime,view.PS.Area(kk).FinalSignalScaled)
                plot(view.PS.Signal.xmc,view.PS.Signal.mmvg,'+c')
                plot(view.PS.Signal.posM,view.PS.Signal.M,'+k')
                plot(view.PS.Signal.xmr,view.PS.Signal.mmvd,'+g')
      
                
               
                plot(pos(:),M(:),'+m','linewidth',2)
                
                
                
                xlabel('time')
                ylabel('displacement')
%                 set(gca,'Xlim',[0,max(view.PS.FinalTimeScaled)],'Ylim',[min(view.PS.Area(kk).FinalSignalScaled),max(view.PS.Area(kk).FinalSignalScaled)])

                if sauv
                    saveas(gcf,[view.SaveFolder,view.SaveFileName,'_Peaks_Signal.png'])
                end
                close
                
            end
 
        function PeaksValleys(view,varargin)
                p = inputParser;

            addOptional(p, 'sauv',1);
            addOptional(p, 'keep',0);
            parse(p, varargin{:});
            sauv=p.Results.sauv;
            keep=p.Results.keep;
            nareas=length(view.PS.Area);
            for kk=1:nareas
                
                pos=view.PS.Area(kk).posper(:,:,1);
                M=view.PS.Area(kk).Mper(:,:,1);
                
               figure
                hold on
               plot(view.PS.Area(kk).FinalTimeScaled,view.PS.Area(kk).FinalSignalScaled)
               % plot(view.PS.FinalTime,view.PS.Area(kk).FinalSignalScaled)
                plot(view.PS.Area(kk).xmc,view.PS.Area(kk).mmvg,'+c')
                plot(view.PS.Area(kk).posM,view.PS.Area(kk).M,'+k')
                plot(view.PS.Area(kk).xmr,view.PS.Area(kk).mmvd,'+g')
      
                
               
                plot(pos(:),M(:),'+m','linewidth',2)
                
                
                
                xlabel('time')
                ylabel('displacement')
%                 set(gca,'Xlim',[0,max(view.PS.FinalTimeScaled)],'Ylim',[min(view.PS.Area(kk).FinalSignalScaled),max(view.PS.Area(kk).FinalSignalScaled)])
               if keep&&nareas>1
                 aa=mvdlg('Should we keep this area (1 for yes, 0 for no)?','ok',[.4 .2 .2 .2]);
%                  aa=inputdlg('Should we keep this area (1 for yes, 0 for no)?')
                 view.PS.goodar(kk)=str2num(aa{1});
               end
               if nareas==1
                   view.PS.goodar(1)=1;
               end
    
                if view.PS.goodar(kk)
                    titre='signal';
                else
                    titre='noise';
                end
                title(['Area   ',num2str(kk), ' : ',titre],'color',ColorSpecToRGB(view.couleur{kk}))
                
                
          
                if sauv
                    saveas(gcf,[view.SaveFolder,view.SaveFileName,'_Peaks_Area',num2str(kk),'.png'])
                end
                close
            end
        end
        
        function Trajectories(view,varargin)
            if ~isempty(varargin)
                sauv=varargin{1};
            else
                sauv=0;
            end
            ind=logical(view.PS.signalarea);
            MX=max(view.PS.Amplitude(ind));
            
            figure;hold on
            %xlim([-max(MX) MX]); ylim([-max(MX) MX]);
            parts = find(view.PS.signalarea);
            N_part = length(parts);
            for p=1:10:N_part
                clear Trajectories

                %[~,I] = sort(view.PS.Signal.xmr(~isnan(view.PS.Signal.xmr)) - view.PS.Signal.xmc(~isnan(view.PS.Signal.xmc)), 'descend');
                %N = sum((view.PS.FinalTimeScaled > view.PS.Signal.xmc(I(1))) & (view.PS.FinalTimeScaled < view.PS.Signal.xmr(I(1))));
                % for pp=1:length(I)
                %    mask = (view.PS.FinalTimeScaled >= view.PS.Signal.xmc(I(pp))) & (view.PS.FinalTimeScaled <= view.PS.Signal.xmr(I(pp)));
                %    time_traj = view.PS.FinalTimeScaled(mask);
                %    x_traj = interp1(time_traj, view.PS.xPosition(parts(p),mask), linspace(time_traj(1), time_traj(end), N));
                %    y_traj = interp1(time_traj, view.PS.yPosition(parts(p),mask), linspace(time_traj(1), time_traj(end), N));
                %    Trajectories(pp,:,:) = [x_traj' ,  y_traj'];
                %end
                %Trajectories = median(Trajectories, 1);

                ind = (view.PS.FinalTimeScaled >= view.PS.Signal.xmc(1)) & (view.PS.FinalTimeScaled <= view.PS.Signal.xmr(1));
                ind=[1:100];
                x_traj = 50*view.PS.xPosition(parts(p),ind)+view.PS.Center(p,1);
                y_traj = 50*view.PS.yPosition(parts(p),ind)+view.PS.Center(p,2);
                Trajectories = [x_traj' ,  y_traj'];
                N=length(x_traj);
                 map=parula(N-1);
                 for t=1:N-1
                     plot([Trajectories(t, 1),Trajectories(t+1, 1)],[Trajectories(t, 2),Trajectories(t+1, 2)], 'color',map(t,:))
                 end
                
                % t=(1:N)';
                % i = size(t);   % t is a column vector
                % 
                % 
                % a0 = [10 10];
                % c = [mean(x_traj) mean(y_traj)];
                % f = @(a) ((x_traj-c(1)).^2)/a(1).^2 + ((y_traj-c(2)).^2)/a(2).^2 -1;
                % options = optimset('Display','iter');
                % af = lsqnonlin(f, a0, [], [], options);
                % af = abs(af);
                % 
                % t_rad = (1:N)*2*pi/N;
                % xnew = c(1) + af(1)*cos(t_rad);
                % ynew = c(2) + af(2)*sin(t_rad);
                % %plot(xnew, ynew ,'r')
                % scatter(min(af(1),af(2)), min(af(1),af(2))./max(af(1),af(2)), 'o'); hold on
            end

            axis equal
            colorbar
            set(gca,'clim',[0 100])
            colormap(gca, 'parula');
     
             if sauv
                saveas(gcf,[view.SaveFolder,view.SaveFileName,'_traj.png'])
             end
             close
        end

        function DisplacementsatPeaks(view,varargin)
            if ~isempty(varargin)
                sauv=varargin{1};
            else
                sauv=0;
            end
            ind=logical(view.PS.signalarea);
            amp=median(view.PS.Amplitude(ind,:),1);
            
            AA=quantile(amp,0.99);
            [pks,Locs]=findpeaks(amp','MinPeakHeight',AA/2);
                    %figure
                    %hold on
                    %plot(amp)
                    %plot(Locs,pks,'+r')
            for uu=1: length(Locs)
                
                ii=round(Locs(uu));
                
                xPos(:,uu)=view.PS.xPosition(:,ii);
                yPos(:,uu)=view.PS.yPosition(:,ii);
                
            end
            figure
            if(exist('xPos','var'))
                hold on
            quiver(view.PS.Center(ind,2),view.PS.Center(ind,1),mean(xPos(ind,:),2),mean(yPos(ind,:),2), 'm')
            quiver(view.PS.Center(~ind,2),view.PS.Center(~ind,1),mean(xPos(~ind,:),2),mean(yPos(~ind,:),2), 'c')
            set(gca, 'Ydir', 'reverse')
            axis equal
            axis off
            
             if sauv
                saveas(gcf,[view.SaveFolder,view.SaveFileName,'_dispatpeak.png'])
             end
             close
            end
        end
    end
    end
    function color = ColorSpecToRGB(s)
    color=[];
    switch s
        case 'y'
            color = [1 1 0];
        case 'm'
            color = [1 0 1];
        case 'c'
            color = [0 1 1];
        case 'r'
            color = [1 0 0];
        case 'g'
            color = [0 1 0];
        case 'b'
            color = [0 0 1];
        case 'w'
            color = [1 1 1];
        case 'k'
            color = [0 0 0];
    end
    
    end