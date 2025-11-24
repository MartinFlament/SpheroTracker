classdef PeriodicSignalsCaView < handle
    properties (SetAccess = private)
        PS
        SaveFolder
        SaveFileName
        couleur={'y','m','c','r','g','b','w','k'};
    end
    methods
        function view = PeriodicSignalsCaView(varargin)
            
            p = inputParser;
            % addRequired(p, 'PS', @(x) isa(x,'PeriodicSignalsModel'));
            addRequired(p, 'PS');
            addOptional(p, 'SaveFolder', @(x) isa(x,'String'));
            addOptional(p, 'SaveFileName', @(x) isa(x,'String'));
            %addParameter(p, 'Save', 'No', @(x) any(validatestring(x, {'Yes','No'})));
            parse(p, varargin{:});
            
            view.PS = p.Results.PS;
            view.SaveFolder = p.Results.SaveFolder;
            view.SaveFileName = p.Results.SaveFileName;
            %view.Save = p.Results.Save;
        end
        
        function circle(view,varargin)
            if ~isempty(varargin)
                sauv=varargin{1};
            else
                sauv=0;
            end
            
            im1=double(imread([view.PS.PathName,view.PS.FileName]));
            
            figure
            hold on
            imagesc(double(im1))
            plot(view.PS.boundary(:,2), view.PS.boundary(:,1), '+g', 'LineWidth', 2)
           
            colormap(gray)
            viscircles([view.PS.z(2) view.PS.z(1)],view.PS.R);
            set(gca,'Ydir','reverse')
            
            axis equal
            axis off
            if sauv
                saveas(gcf,[view.SaveFolder,view.SaveFileName,'_circle.png'])
                close
            end
            
        end
        
        
        
        function Boxes(view,varargin)
            if ~isempty(varargin)
                sauv=varargin{1};
            else
                sauv=0;
            end
            im=imread([view.PS.PathName,view.PS.FileName]);
            map=jet(length(view.PS.vecr));
            figure
            
            hold on
            imagesc(im)
            colormap('gray')
            axis equal
            axis tight
            axis off
            set(gca,'Ydir','reverse')
            cc=1;
            for r=view.PS.vecr
                dtheta=view.PS.Box/r;
                for theta=0:dtheta:2*pi
                    cbox=[view.PS.z(2)+r*cos(theta) view.PS.z(1)+r*sin(theta)];
                    posrect=[cbox(1)-view.PS.Box/2 cbox(2)-view.PS.Box/2 view.PS.Box view.PS.Box];
                    if cc==1
                        h=rectangle('Position',posrect);
                        h.EdgeColor=map(cc,:);
                    elseif cc==length(view.PS.vecr)
                        h=rectangle('Position',posrect);
                        h.EdgeColor=map(cc,:);
                    end
                end
                cc=cc+1;
            end
            
            if sauv
                saveas(gcf,[view.SaveFolder,view.SaveFileName,'_Boxes.png'])
                close
            end
            
            
        end
        
        
        function FourierGraph(view,varargin)
            if ~isempty(varargin)
                sauv=varargin{1};
            else
                sauv=0;
            end
            figure('Name','Fourier')
            %figToolbarFix
            hold on
            plot(view.PS.fvec(2:(length(view.PS.F)-1)/2),abs(view.PS.F(2:(length(view.PS.F)-1)/2)))
            plot(view.PS.fvec(view.PS.pk1),abs(view.PS.F(view.PS.pk1)),'r+')
            plot(view.PS.fvec(view.PS.pk2),abs(view.PS.F(view.PS.pk2)),'r+')
            ylabel('Intensity')
            xlabel('frequency')
            set(gca,'Xlim',[0 40])
            
            
            if sauv
                saveas(gcf,[view.SaveFolder,view.SaveFileName,'_Fourier.png'])
                close
            end
        end
        
        
        function IntensAngular(view,varargin)
            if ~isempty(varargin)
                sauv=varargin{1};
            else
                sauv=0;
            end
            
            map=gray(length(unique(view.PS.intensth)));
           h=figure;
           
            hold on
            ylabel('time (s)')
            xlabel('Angular coordinate (rad)')
            %set(gca,'Xlim',[0 6.1])
            cp=1;
            M=size(view.PS.intensth,1);
            imagesc(view.PS.intensth')
            
            % legend(leg,'Location','eastoutside')
            set(gca, 'Ydir', 'reverse')
            disp(size(view.PS.intensth))
            set(gca, 'Xlim', [0.5 size(view.PS.intensth, 1)+0.5])
            set(gca, 'Ylim', [0 size(view.PS.intensth, 2)])
            
            set(get(handle(gcf),'JavaFrame'),'Maximized',1);
            %pause
            if sauv
                saveas(gcf,[view.SaveFolder,view.SaveFileName(1:end-5),'_Angular_intensity.png'])
                close
            end
            
        end

        function IntensRadial(view,varargin)
            if ~isempty(varargin)
                sauv=varargin{1};
            else
                sauv=0;
            end
            
            map=jet(length(view.PS.vecr));
           h=figure;
           
            hold on
            xlabel('time (s)')
            ylabel('Normalized Intensity')
            %set(gca,'Xlim',[0 6.1])
            cp=1;
            M=size(view.PS.intensr,1);
            for ii=1:2:M
                figure(h)
                plot(view.PS.Time,view.PS.intensr(ii,:),'color',map(ii,:))
                leg{cp}=['R=',num2str(round(view.PS.vecr(ii)*view.PS.PixelSize))];
                
                cp=cp+1;
            end
            
            legend(leg,'Location','eastoutside')
            set(get(handle(gcf),'JavaFrame'),'Maximized',1);
            %pause
            if sauv
                saveas(gcf,[view.SaveFolder,view.SaveFileName,'_Radial_intensity.png'])
                close
            end
            
        end
        function PlotPhaseRadial(view,varargin)
            if ~isempty(varargin)
                sauv=varargin{1};
            else
                sauv=0;
            end

            figure('Units','Normalized','Position',[0.2 0.2 0.5 0.8 ])
  
            hold on
            %yyaxis left
            plot(view.PS.vecr*view.PS.PixelSize,view.PS.phase1)
            ylabel('$\Phi_1$')
            %yyaxis right
            %plot(view.PS.vecr*view.PS.PixelSize,view.PS.phase2)
            xlabel('r ($\mu$m)')
            %ylabel('$\Phi_2$')
             set(gca,'Xlim',[0,max(view.PS.vecr*view.PS.PixelSize)])
            %pause
            if sauv
                saveas(gcf,[view.SaveFolder,view.SaveFileName,'_Fourier_phase.png'])
                close
            end
            
        end
        
        function PlotIntPkFourier(view,varargin)
            if ~isempty(varargin)
                sauv=varargin{1};
            else
                sauv=0;
            end
            figure('Units','Normalized','Position',[0.2 0.2 0.5 0.8 ])
            
           % figure('Position',[0 0 0.5 0.5])
            hold on
            yyaxis left
            plot(view.PS.vecr*view.PS.PixelSize,view.PS.Intpk2)
            ylabel('Peak Amplitude')
            yyaxis right
            plot(view.PS.vecr*view.PS.PixelSize,view.PS.Intpk1)
             xlabel('r ($\mu$m)')
             ylabel('Peak Amplitude')
             set(gca,'Xlim',[0,max(view.PS.vecr*view.PS.PixelSize)])
            if sauv
                saveas(gcf,[view.SaveFolder,view.SaveFileName,'_Fourier_int.png'])
                close
            end
        end
     
                function PeaksValleys(view,varargin)
                p = inputParser;

            addOptional(p, 'sauv',1);
            addOptional(p, 'keep',0);
            parse(p, varargin{:});
            sauv=p.Results.sauv;
            keep=p.Results.keep;
            nareas=length(view.PS.Area);
           
                 kk=1;
                pos=view.PS.Signal.posper(:,:,1);
                M=view.PS.Signal.Mper(:,:,1);
               
               figure
                hold on
               plot(view.PS.Time,view.PS.Global1)
               % plot(view.PS.FinalTime,view.PS.Area(kk).FinalSignalScaled)
                plot(view.PS.Signal.xmc,view.PS.Signal.mmvg,'+c')
                plot(view.PS.Signal.posM,view.PS.Signal.M,'+k')
                plot(view.PS.Signal.xmr,view.PS.Signal.mmvd,'+g')
      
                
               
                plot(pos(:),M(:),'+m','linewidth',2)
                
                
                
                xlabel('time')
                ylabel('displacement')
%                 set(gca,'Xlim',[0,max(view.PS.FinalTimeScaled)],'Ylim',[min(view.PS.Area(kk).FinalSignalScaled),max(view.PS.Area(kk).FinalSignalScaled)])
            
              
                
                
                if sauv
                    saveas(gcf,[view.SaveFolder,view.SaveFileName,'_Peaks.png'])
                    close
                end
                
            
        end
        
        
    end
end
