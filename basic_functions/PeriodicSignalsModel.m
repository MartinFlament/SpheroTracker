classdef PeriodicSignalsModel < handle
    % properties (SetAccess = private)
    properties
        Position
        NumberOfParticles
        wel
        Center
        xPosition
        yPosition
        Amplitude
        Amplitude1
        Amplitude2
        FourierPeriod
        FourierAmplitude
        Period
        SFFT
        Step=1
        InterpTime
        InterpSignal
        ShiftTime
        ShiftSignal
        StartPoint
        FinalSignal
        FinalTimeScaled
        FinalSignalScaled
        FinalTime
        PhaseShift
        SamplingFrequency
        ExcitationFrequency
        PixelSize
        signalarea
        Signal
        Noise
        Strain
        T
        good
        goodar
        Area
        framerate
        Homogeneity
        TotalPoints
        TaucMedMoy
        TaucStdMoy
        TaurMedMoy
        TaurStdMoy
        TausMedMoy
        TausStdMoy
        TausBazMedMoy
        TausBazStdMoy
        TaudMedMoy
        TaudStdMoy
        TaudBazMedMoy
        TaudBazStdMoy
        AireMedMoy
        AireStdMoy
        pcMedMoy
        pcStdMoy
        prMedMoy
        prStdMoy
        hrMedMoy
        hrStdMoy
        hcMedMoy
        hcStdMoy
        MMedMoy
        MStdMoy
        minMedMoy
        minStdMoy
        RestMedMoy
        RestStdMoy
        Pressure
        Shear
        Vorticity
        Thampl
        amp
        results_pathname
        col_line
  
    end
    methods
        function PS = PeriodicSignalsModel(varargin)
            p = inputParser;
            addRequired(p, 'Position', @(x) ndims(x) == 3);
            addOptional(p, 'PixelSize', 1);
            addOptional(p, 'SamplingFrequency', 1);
            addOptional(p, 'ExcitationFrequency', 1);
            addOptional(p, 'results_pathname','');
            addOptional(p, 'col_line',0);
            parse(p, varargin{1:end});
            PS.Position = p.Results.Position;
            PS.PixelSize = p.Results.PixelSize;
            PS.results_pathname=p.Results.results_pathname;
            PS.NumberOfParticles = size(PS.Position, 1);
            PS.SamplingFrequency=p.Results.SamplingFrequency;
            PS.ExcitationFrequency=p.Results.ExcitationFrequency;
            PS.NumberOfParticles = size(PS.Position, 1);
            PS.framerate = (PS.SamplingFrequency)/1000; %fps in ms
            PS.col_line=p.Results.col_line;
        end
        
        function PS = determineCenter(PS, varargin)            

            pos_x = PS.Position(:,:,1);
            pos_y = PS.Position(:,:,2);

           %  gradient_x  = diff(pos_x,1,2);
           %  gradient_y  = diff(pos_y,1,2);
           %  gradient_r  = diff(PS.rPosition(PS.rPosition(:,1)>quantile(PS.rPosition(:,1),0.7),:),1,2);
           %  norm = gradient_x.^2 +gradient_y.^2;
           %  %norm = sum(gradient_r.^2,1);
           %  ind=PS.rPosition(:,1)>quantile(PS.rPosition(:,1),0.8);
           % 
           %  ind1 = PS.rPosition(:,1)>quantile(PS.rPosition(:,1),0.8)&thPosition(:,1)>0;
           % 
           %  ind2 = PS.rPosition(:,1)>quantile(PS.rPosition(:,1),0.8)&thPosition(:,1)<0;
           %  % rot1 = (pos_x(ind,1:end-1) - median(median(pos_x,1))).*gradient_x(ind,:) + (pos_y(ind,1:end-1) - median(median(pos_y,1))).*gradient_y(ind,:);
           %  % rot2=(pos_x(ind1,1:end-1) - median(median(pos_x,1))).*gradient_y(ind1,:) - (pos_y(ind1,1:end-1) - median(median(pos_y,1))).*gradient_x(ind1,:);
           %  % rot3=(pos_x(ind2,1:end-1) - median(median(pos_x,1))).*gradient_y(ind2,:) - (pos_y(ind2,1:end-1) - median(median(pos_y,1))).*gradient_x(ind2,:);
           %  % figure;plot(abs(sum(rot1,1)),'r');hold on; plot(abs(sum(rot2,1)),'b'),plot(abs(sum(rot3,1)),'g');pause;close;
           %  % figure; hold on; plot(abs(sum(gradient_x(ind,:),1)),'g');plot(abs(sum(gradient_y(ind,:),1)),'k');
           % 
           % % [pks, ind_max] = findpeaks(abs(smooth(sum(rot,1))),'Npeaks', round(size(PS.Position,2)*PS.ExcitationFrequency/PS.SamplingFrequency),'MinPeakDistance',0.9*PS.SamplingFrequency/PS.ExcitationFrequency);
           %  [mm1,ind1_max]=max(abs(sum(rot1,1)));
           %  [mm2,ind2_max]=max(abs(sum(rot2,1)));
           %  if mm2>mm1
           %      ind_max=ind2_max;
           %  else
           %      ind_max=ind1_max;
           %  end
            % xCenterMedian = median(PS.Position(:, round((ind_max(2:end)+ind_max(1:end-1))/2), 1),2);
            % yCenterMedian = median(PS.Position(:, round((ind_max(2:end)+ind_max(1:end-1))/2), 2),2);
            
            if ~isempty(varargin)
                Var=varargin{1};
                [pk,ind_max] = max(Var(10:end));

                %figure;plot(Var);hold on;scatter(ind_max+9,pk, 'r+');scatter(ind_max,Var(ind_max), 'g+');pause;close
                xCenterMedian = median(PS.Position(:, ind_max:ind_max+9, 1),2);
                yCenterMedian = median(PS.Position(:, ind_max:ind_max+9, 2),2);
                PS.Center = [xCenterMedian, yCenterMedian];
                PS.xPosition = PS.Position(:, :, 1) - PS.Center(:, 1);
                PS.yPosition = PS.Position(:, :, 2) - PS.Center(:, 2);
                PS.Amplitude = sqrt((PS.Position(:, :, 1) - PS.Center(:, 1)).^2 + (PS.Position(:,:, 2) - PS.Center(:, 2)).^2);
            else
                for p = 1: PS.NumberOfParticles
                    %figure
                    [N, C] = hist3(squeeze(PS.Position(p, :, :)));
                    % imagesc(N)
                    [~, IND] = max(N(:));
                    [I, J] = ind2sub(size(N), IND);
                    xmin = C{1}(I) - mode(gradient(C{1}))/2;
                    xmax = C{1}(I) + mode(gradient(C{1}))/2;
                    ymin = C{2}(J) - mode(gradient(C{2}))/2;
                    ymax = C{2}(J) + mode(gradient(C{2}))/2;
                    
                    
                    PS.wel(p, :) = PS.Position(p, :, 1) > xmin & PS.Position(p, :, 1) < xmax & PS.Position(p, :, 2) > ymin & PS.Position(p, :, 2) < ymax;
                    
                    %                 xCenterMean = mean(PS.Position(p, logical(PS.wel(p, :)), 1));
                    %                 yCenterMean = mean(PS.Position(p, logical(PS.wel(p, :)), 2));
                    
                     xCenterMedian = median(PS.Position(p,logical(PS.wel(p, :)), 1));
                     yCenterMedian = median(PS.Position(p, logical(PS.wel(p, :)), 2));
                     PS.Center(p, :) = [xCenterMedian, yCenterMedian];
                    
                    
                    PS.xPosition(p, :) = PS.Position(p, :, 1) - PS.Center(p, 1);
                    PS.yPosition(p, :) = PS.Position(p, :, 2) - PS.Center(p, 2);
                    PS.Angle(p) = median(atan2(PS.yPosition(p, ~PS.wel(p, :)), PS.xPosition(p, ~PS.wel(p, :))));
                    u(1)=cos(PS.Angle(p));
                    u(2)=sin(PS.Angle(p));
                    PS.Amplitude(p, :) = sqrt((PS.Position(p, :, 1) - PS.Center(p, 1)).^2 + (PS.Position(p,:, 2) - PS.Center(p, 2)).^2);
                    PS.Amplitude1(p, :) = PS.xPosition(p,:)*u(1)+PS.yPosition(p,:)*u(2);
                    PS.Amplitude2(p, :) = PS.xPosition(p,:)*u(2)-PS.yPosition(p,:)*u(1);
                    
                    
                end
            end
        end
        
        function PS = determineCenter_zebra(PS)
            figure
            axis equal
            hold on
            set(gca,'Ydir','reverse')
            for p = 1 :10: PS.NumberOfParticles
                
                plot(PS.Position(p, :, 2),PS.Position(p,:,1),'.')
            end       
        end

        function PS = findFourierPeriod(PS, varargin)
            Signal = varargin{1};
            SamplingFrequency = varargin{2};
            for p = 1 : size(Signal,1)
                
%                 figure
%                 hold on
%                 plot(Signal(p,:))
                FA = FourierAnalysisModel(Signal(p, :), 'SamplingFrequency', SamplingFrequency);
                FA.performFft;
                PS.PhaseShift(p)=FA.PhaseShift;
                PS.FourierPeriod(p) = FA.Period;
                PS.FourierAmplitude(p) = FA.Amplitude;
            end
            %plot(PS.FourierPeriod, PS.FourierAmplitude, 'o');pause;close
            PS.Period = mode(PS.FourierPeriod);
            % PS.Homogeneity=nanstd(PS.PhaseShift);
        end
        
        
        function PS = interpolateSignal(PS, varargin)
            Time = varargin{1};
            Signal = varargin{2};
            PS.Step = varargin{3};
            
            %TimeGrad = gradient(Time);
            %PS.Step = max([prctile(TimeGrad(:), 10), 0.001]); %Is this a good way?
            PS.InterpTime = 0 : PS.Step : max(Time(:));
            PS.InterpSignal = NaN(PS.NumberOfParticles, length(PS.InterpTime));
            
            for p = 1 : 1 : size(Signal,1)
                clear F NT NS ia
                
                [NT, ia, ~] = unique(Time(p, :));
                NS = Signal(p, ia);
                
                F = griddedInterpolant(NT, NS);
                PS.InterpSignal(p, :) = F(PS.InterpTime);
            end
        end
        function PS = interpolateSignalF(PS, varargin)
            Time = varargin{1};
            Signal = varargin{2};
            PS.Step = varargin{3};
            
            %TimeGrad = gradient(Time);
            %PS.Step = max([prctile(TimeGrad(:), 10), 0.001]); %Is this a good way?
            PS.InterpTime = 0 : PS.Step : max(Time(:));
            PS.InterpSignal = NaN(PS.NumberOfParticles, length(PS.InterpTime));
            L=length(PS.InterpTime);
            
            
            for p = 1 : size(Signal,1)
                ll=L-length(Signal);
%                 figure
%                 plot(Signal(p,:))
                                SignalF=fftshift(fft(Signal(p,:)));
%                 SignalF=(fft(Signal(p,:)));
                if mod(ll,2)==0
                    Signal_pad=[zeros(1,ll/2),SignalF(:),zeros(1,ll/2)];
                    %                 Signal_pad2=[zeros(ll/2,1);Signal(p,:);zeros(ll/2)];
                else
                    Signal_pad=[zeros(1,(ll-1)/2),SignalF,zeros(1,(ll+1)/2)];
                    
                    %                 Signal_pad2=[zeros((ll-1)/2,1);Signal(p,:,2);zeros((ll+1)/2)];
                    
                end
                           PS.InterpSignal(p,:)= real(ifft((ifftshift(Signal_pad))))*10;

%                 figure
%                 plot(PS.InterpSignal(p,:))
%                 pause
              
            end
        end
        
        
        function PS = findPhaseShift(PS, varargin)
            Time = varargin{1};
            Signal = varargin{2};
            if length(varargin{3}) == 1
                Period = repmat(varargin{3}, PS.NumberOfParticles);
            else
                Period = varargin{3};
            end
            
            
            for p = 1 : 1 : size(Signal,1)
                % clear tau  x v locs
                if abs(PS.Period-PS.FourierPeriod(p))/PS.Period>0.05
                    
                    FA = FourierAnalysisModel(PS.Amplitude(p,:), 'SamplingFrequency', PS.SamplingFrequency);
                    FA.performFft(PS.Period);
                    %FA.findBasePeriod;
                    PS.PhaseShift(p)=FA.PhaseShift;
                    %               PS.FourierPeriod(p) = FA.Period;
                    %   PS.FourierAmplitude(p) = FA.Amplitude;
                end
                PS.StartPoint(p)=round(PS.PhaseShift(p).*PS.Period/PS.Step*PS.SamplingFrequency/2/pi);
                %   PS.StartPoint(p)=round(PS.PhaseShift(p).*PS.Period/2/pi);
                ll=length(Signal(p,:));
                if PS.StartPoint(p)>0
                    PS.ShiftSignal(p, :) = [NaN(1,PS.StartPoint(p)) Signal(p, 1 : ll-PS.StartPoint(p))];
                else
                    if ~isnan(PS.StartPoint(p))
                        PS.ShiftSignal(p, :) = [Signal(p,  -PS.StartPoint(p)+1 :end),NaN(1,-PS.StartPoint(p))];
                    else
                        PS.ShiftSignal(p, :) =NaN;
                    end
                end
                
                
                PS.ShiftTime = Time;
            end
            PS.Homogeneity=nanstd(PS.PhaseShift);
        end
        
        % function PS = measure_strain_pressure_vorticity(PS,varargin)
        %     mask = varargin{1};
        %     params = varargin{2};
        %     sz=size(mask);
        %     pas=param.minDist;
        % 
        %     [X,Y]=meshgrid(1:pas:sz(2),1:pas:sz(1));

        function PS = ExtractSignal(PS, varargin)
            Per_thresh = otsuthresh(PS.FourierPeriod);
            PS.signalarea = PS.FourierPeriod > Per_thresh;
            
            if sum(PS.signalarea)>0
                PS_shift = PS.PhaseShift - PS.PhaseShift';
                PS_shift(PS_shift<-pi) = mod(PS_shift(PS_shift<-pi), pi);
                PS_shift(PS_shift>pi) = mod(PS_shift(PS_shift>pi), -pi);
                [PS.Homogeneity,I] = min(nanstd(PS_shift),[],2);
                PS_tot = PS_shift(I,:);
                % PS_plus = mod(PS.PhaseShift-m, pi);
                % PS_moins = mod(PS.PhaseShift-m, -pi);
                % ind=abs(PS_plus)>abs(PS_moins);
                % PS_tot=PS_plus;
                % PS_tot(ind)= PS_moins(ind);

                PS.Signal.FinalSignaln=median(PS.ShiftSignal(PS.signalarea,:),1);
                PS.Signal.FinalSignal=PS.Signal.FinalSignaln(~isnan(PS.Signal.FinalSignaln));
                % PS.Signal.FinalSignal=mean(PS.ShiftSignal(PS.signalarea,:),1);
                PS.Signal.RawFinalSignal=median(PS.Amplitude(PS.signalarea,:),1);
                PS.Signal.Period=mode(mode(PS.FourierPeriod(PS.signalarea)));
                PS.Signal.PhaseShift=median(PS.PhaseShift(PS.signalarea));
                PS.Signal.StdPhaseShift=std(PS.PhaseShift(PS.signalarea))/sqrt(sum(PS.signalarea));
                PS.Signal.PointNumb = sum(PS.signalarea);
                PS.Signal.FinalSignalScaled=PS.Signal.FinalSignal*PS.PixelSize;
                PS.Signal.PointNumbPercent = round(((sum(PS.signalarea))*100)./(PS.NumberOfParticles));
                tata=PS.FinalTimeScaled(~isnan(PS.Signal.FinalSignaln));
                titi=PS.FinalTime(~isnan(PS.Signal.FinalSignaln));
                PS.Signal.FinalTimeScaled=tata-min(tata);
                PS.Signal.FinalTime=titi-min(titi);
                
                % figure;
                % hist(PS_tot, 50)
                % figure;
                % hist(PS_tot(PS.signalarea), 50)
                % pause

                PS.Signal.Homogeneity=nanstd(PS_tot(PS.signalarea));
                
                
                p=inputParser;
                addOptional(p, 'Smoothlength', 80);
                addOptional(p, 'prop', 0.05);
                parse(p, varargin{1:end});
                prop=p.Results.prop;
                sl=p.Results.Smoothlength;
                list_param_name={'N_pks','Period','Asc_time','Decay_time',...
                            'Decay_time_95','Decay_time_90','Decay_time_70','Decay_time_50','Decay_time_30','Decay_time_20',...
                            'Taud','Baz_taud','AUC','Asc_slope','Decay_slope',...
                            'Amp_asc','Amp_decay','Maxima','Minima'};
                list_calc =logical([0,1,1]);% for median and standard calculation
                
                    MM(:,1)=PS.Signal.FinalTimeScaled;
                    MM(:,2)=PS.Signal.FinalSignal;
                    PK=AnalysisPeaks(MM,'list_param_name',list_param_name,'list_calc', list_calc,'PixelSize',PS.PixelSize,'smoothness',sl,'prop',prop);
                    PK.matrix_filtered_fluorescences(:,1)=MM(:,2);
                    PK.PeriodF=PS.Signal.Period;
                    PK.Homogeneity=PS.Signal.Homogeneity;
                    PK.col_line=PS.col_line;
                    PK.CalculateParameters(1);
                    PS.Signal = PK.Save_in_struct(PS.Signal);
                    PK.sheet=1;
                    PK.Save_sheet(PS.results_pathname);
                
     
                clear MM
                
                if sum(~PS.signalarea)>0
                PS.Noise.FinalSignaln=median(PS.ShiftSignal(~PS.signalarea,:),1);
                PS.Noise.FinalSignal=PS.Noise.FinalSignaln(~isnan(PS.Noise.FinalSignaln));
                % PS.Noise.FinalSignal=mean(PS.ShiftSignal(PS.Noisearea,:),1);
                PS.Noise.RawFinalSignal=median(PS.Amplitude(~PS.signalarea,:),1);
                PS.Noise.Period=mode(mode(PS.FourierPeriod(~PS.signalarea)));
                PS.Noise.PhaseShift=median(PS.PhaseShift(~PS.signalarea));
                PS.Noise.StdPhaseShift=std(PS.PhaseShift(~PS.signalarea))/sqrt(sum(~PS.signalarea));
                PS.Noise.PointNumb = sum(~PS.signalarea);
                PS.Noise.FinalSignalScaled=PS.Noise.FinalSignal*PS.PixelSize;
                PS.Noise.PointNumbPercent = round(((sum(~PS.signalarea))*100)./(PS.NumberOfParticles));
                tata=PS.FinalTimeScaled(~isnan(PS.Noise.FinalSignaln));
                titi=PS.FinalTime(~isnan(PS.Noise.FinalSignaln));
                PS.Noise.FinalTimeScaled=tata-min(tata);
                PS.Noise.FinalTime=titi-min(titi);
                PS.Noise.Homogeneity=nanstd(PS_tot(~PS.signalarea));
                
                
                p=inputParser;
                addOptional(p, 'Smoothlength', 80);
                addOptional(p, 'prop', 0.05);
                parse(p, varargin{1:end});
                prop=p.Results.prop;
                sl=p.Results.Smoothlength;
                list_param_name={'N_pks','Period','Asc_time','Decay_time',...
                            'Decay_time_95','Decay_time_90','Decay_time_70','Decay_time_50','Decay_time_30','Decay_time_20',...
                            'Taud','Baz_taud','AUC','Asc_slope','Decay_slope',...
                            'Amp_asc','Amp_decay','Maxima','Minima'};
                list_calc =logical([0,1,1]);% for median and standard calculation
                
          
                    MM(:,1)=PS.Noise.FinalTimeScaled;
                    MM(:,2)=PS.Noise.FinalSignal;
                    PK=AnalysisPeaks(MM,'list_param_name',list_param_name,'list_calc', list_calc,'PixelSize',PS.PixelSize,'smoothness',sl,'prop',prop);
                    PK.matrix_filtered_fluorescences(:,1)=MM(:,2);
                    PK.PeriodF=PS.Noise.Period;
                    PK.Homogeneity=PS.Noise.Homogeneity;
                    PK.col_line=PS.col_line;
                    PK.CalculateParameters(1);
                    PS.Noise = PK.Save_in_struct(PS.Noise);
                    PK.sheet=1;
                    PK.Save_sheet(PS.results_pathname);
                
     
                clear MM
                end
            end
        end

        function PS=measureStrain(PS, varargin)
            mask_pathname=varargin{1};
            param_pathname=varargin{2};

            load(mask_pathname);
            load(param_pathname);
            sz=size(mask);
            pas=param.minDist;
            [X,Y]=meshgrid(1:pas:sz(2),1:pas:sz(1));
            r = sqrt((PS.Center(:,1)-median(PS.Center(:,1))).^2 + (PS.Center(:,2)-median(PS.Center(:,2))).^2);
            th = atan2(PS.Center(:,2), PS.Center(:,1));
            ur=griddata(PS.Center(:,1),PS.Center(:,2),r,X,Y,'cubic');
            uth=griddata(PS.Center(:,1),PS.Center(:,2),th,X,Y,'cubic');
            
            u_mask = griddata(PS.Center(:,1),PS.Center(:,2),double(PS.signalarea'),X,Y,'linear');
            ind_mask = ~isnan(u_mask);
            ind_signal = logical(u_mask(ind_mask));
            PS.Strain.mask = ind_signal;
            PS.Strain.rPosition = ur(~isnan(ur));
            PS.Strain.thPosition = uth(~isnan(uth));
            
            for ii=1:size(PS.Amplitude,2)
                ux=griddata(PS.Center(:,1),PS.Center(:,2),PS.xPosition(:,ii),X,Y,'cubic');
                uy=griddata(PS.Center(:,1),PS.Center(:,2),PS.yPosition(:,ii),X,Y,'cubic');
        
                [epsxx,epsxy]=gradientN_h2(ux,'optimal5',0.1);   
                [epsyx,epsyy]=gradientN_h2(uy,'optimal5',0.1);
            
                epsxy=1/2*(epsxy+epsyx);
                epsxx(isnan(epsxx)) = 0;
                epsyy(isnan(epsyy)) = 0;
                epsxy(isnan(epsxy)) = 0;
                [vec, Dt, Pc] = arrayfun(@eigenstrain2Df,epsxx,epsyy,epsxy,'UniformOutput',0); %Pc composante symétrique Dt valeur propre (shear)
                P = cell2mat(Pc);
                P = P(ind_mask);
                Press(:, ii) = P(ind_signal);
            
               D=cell2mat(Dt);
               D = D(ind_mask);
               Shear(:, ii) = D(ind_signal);
        
               [curlz, cav] = curl(X,Y,ux,uy);
               cur = curlz(ind_mask);
               V(:,ii)=cur(ind_signal);
            end
            PS.Strain.Pressure=Press;
            PS.Strain.Shear=Shear;
            PS.Strain.Vorticity = V;
           
        end

        function PS=makeAreas(PS,varargin)
            
            
            if isempty(varargin)
                nareas=5;
            else
                nareas=varargin{1};
            end
            
            
            
            V(:,1)=PS.Center(:,2);
            V(:,1)=(V(:,1)-nanmean(V(:,1)))/nanstd(V(:,1));
            V(:,2)=PS.Center(:,1);
            V(:,2)=(V(:,2)-nanmean(V(:,2)))/nanstd(V(:,2));
            V(:,3)=PS.PhaseShift';
            factor=5;% factor which gives the relative importance of phaseshift for clustering
            
            V(:,3)=(V(:,3)-nanmean(V(:,3)))/nanstd(V(:,3))*factor;
            
            %     V(:,4)=PS.FourierPeriod';
            %     V(:,4)=(V(:,4)-mean(V(:,4)))/std(V(:,4));
            V(:,5)=max(PS.Amplitude,[],2);
            %V(:,5)=median(PS.Amplitude,2);
            %             figure
            %             histogram(V(:,5))
            %             th=0.25*max(V(:,5)); % set the minimum displacement for taking into account a point.
            th=0.1;
            PS.Thampl=th;
            ind=V(:,5)<th;
            % ind=True(1:length(V));
            V(:,5)=(V(:,5)-nanmean(V(:,5)))/nanstd(V(:,5));
            V(ind,3)=NaN;
            
            PS.T=kmeans(V(:,[1 2 3 ]),nareas);
            PS.TotalPoints = sum(~isnan(PS.T));
            PS.goodar=false(1,nareas);
            cc=1;
            for kk=1:nareas
                ind=PS.T==kk;
                if sum(ind)>3
                    PS.Area(cc).FinalSignaln=median(PS.ShiftSignal(ind,:),1);
                    PS.Area(cc).FinalSignal=PS.Area(cc).FinalSignaln(~isnan(PS.Area(cc).FinalSignaln));
                    % PS.Area(cc).FinalSignal=mean(PS.ShiftSignal(ind,:),1);
                    PS.Area(cc).RawFinalSignal=median(PS.Amplitude(ind,:),1);
                    PS.Area(cc).Period=mode(mode(PS.FourierPeriod(ind)));
                    PS.Area(cc).PhaseShift=median(PS.PhaseShift(ind));
                    PS.Area(cc).StdPhaseShift=std(PS.PhaseShift(ind))/sqrt(sum(ind));
                    PS.Area(cc).PointNumb = sum(ind);
                    PS.Area(cc).FinalSignalScaled=PS.Area(cc).FinalSignal*PS.PixelSize;
                    PS.Area(cc).PointNumbPercent = round(((sum(ind))*100)/(PS.TotalPoints));
                    tata=PS.FinalTimeScaled(~isnan(PS.Area(cc).FinalSignaln));
                    titi=PS.FinalTime(~isnan(PS.Area(cc).FinalSignaln));
                    PS.Area(cc).FinalTimeScaled=tata-min(tata);
                    PS.Area(cc).FinalTime=titi-min(titi);
                    PS.Area(cc).Homogeneity=nanstd(PS.PhaseShift(ind));    
                    
%                     diff_freq=std(PS.FourierPeriod(ind))/ PS.Area(cc).Period
                    diff_freq=nanmedian(abs(PS.FourierPeriod(ind)-PS.Period)/ PS.Period);
                    
                    crit=diff_freq<0.2&&PS.Area(cc).Period>30;%% Be careful this number will highly depend on the experiment!!!!Ideally to be changed
                    
                    
   
                
                    PS.T(ind)=cc;
                    PS.good(ind)=crit;
                    PS.goodar(kk)=crit;
                    cc=cc+1;
                end
            end
            
            
        end
        
        function PS = averageSignal(PS)
            
            PS.FinalSignal = nanmedian(PS.ShiftSignal, 1);
            PS.FinalTime = PS.ShiftTime;
            PS.FinalTimeScaled = PS.FinalTime/PS.SamplingFrequency;
            PS.FinalSignalScaled = PS.FinalSignal*PS.PixelSize;
            
        end
        
        function PS=SlidingFFT(PS,varargin)
            % argument varargin peut contenir le num?ro de la zone sur
            % laquelle faire la sliding fft
            
            if isempty(varargin)
                Signal=PS.FinalSignal;
                ffm=1./(PS.Period)*PS.Step/PS.SamplingFrequency;
            else
                narea=varargin{1};
                Signal=PS.Area(narea).FinalSignal;
                ffm=1./(PS.Area(narea).Period)*PS.Step/PS.SamplingFrequency;
            end
            pp=nextpow2(length(Signal))+3;
            kkp=round((ffm*2^pp));
            if kkp>1000
                np=6;
            else
                np=3;
            end
            box=round(1./ffm*np);
            %% make box even
            box=box+mod(box,2);
            
            larg=round(kkp/4);
            
            range=kkp-larg:kkp+larg;
            val=NaN(1,length(Signal)-box);
            
            for ii=1:1:length(Signal)-box
                
                
                
                %% First step finding precisely the position of the peak
                tata=Signal(ii:ii+box);
                %                 figure
                %                 plot(tata)
                pf1=abs(fft(tata,2^pp)).^2;
                
                lls=length(pf1);
                pf=pf1(1:lls/2);
                pf(1:10)=0;
                [pk,loc]=findpeaks(pf(range));
                %                 figure
                %                 plot(pf(range))
                %                 pause
                [pks,I]=sort(pk);
                locs=loc(I);
                if isempty(loc)
                    mm(ii)=0;
                    pos(ii)=0;
                    val(ii)=0;
                else
                    if locs(1)==1
                        mm(ii)=pks(2);
                        pos(ii)=locs(2);
                    else
                        mm(ii)=pks(1);
                        pos(ii)=locs(1);
                    end
                end
                
                val(ii)=sum(pf(range))/sum(pf(:));
                clear pk loc I pks locs
                
                
            end
            
            ff=kkp+pos-larg-1;
            th=0.75*(nanmedian(val(:)));
            as=val<0.75*(nanmedian(val(:)));
            rate=sum(as)/length(as)*100;
            
            if isempty(varargin)
                PS.SFFT.rate=rate;
                PS.SFFT.mm=mm;
                PS.SFFT.val=val;
                PS.SFFT.as=as;
                PS.SFFT.box=box;
            else
                PS.Area(narea).SFFT.rate=rate;
                PS.Area(narea).SFFT.mm=mm;
                PS.Area(narea).SFFT.val=val;
                PS.Area(narea).SFFT.as=as;
                PS.Area(narea).SFFT.box=box;
            end
            
            
            
        end
        
        
        function PS = calculateParameters(PS,varargin)
            
            % We will create a caldyn structure PK that will contain the
            % signals of the different areas
            p=inputParser;
            addOptional(p, 'Smoothlength', 80);
            addOptional(p, 'prop', 0.05);
            parse(p, varargin{1:end});
            prop=p.Results.prop;
            sl=p.Results.Smoothlength;
            list_param_name={'N_pks','Period','Asc_time','Decay_time',...
                        'Decay_time_95','Decay_time_90','Decay_time_70','Decay_time_50','Decay_time_30','Decay_time_20',...
                        'Taud','Baz_taud','AUC','Asc_slope','Decay_slope',...
                        'Amp_asc','Amp_decay','Maxima','Minima'};
            list_calc =logical([0,1,1]);% for median and standard calculation
            
      
            for kk=1:length(PS.Area)
                MM(:,1)=PS.Area(kk).FinalTimeScaled;
                MM(:,2)=PS.Area(kk).FinalSignal;
                PK=AnalysisPeaks(MM,'list_param_name',list_param_name,'list_calc', list_calc,'PixelSize',PS.PixelSize,'smoothness',sl,'prop',prop);
                PK.matrix_filtered_fluorescences(:,1)=MM(:,2);
                PK.PeriodF=PS.Area(kk).Period;
                PK.Homogeneity=PS.Area(kk).Homogeneity;
                PK.col_line=PS.col_line;
                PK.CalculateParameters(1);
                PS.Area(kk) =PK.Save_in_struct(PS.Area(kk));
                PK.sheet=1;
                PK.area=kk;
                PK.Save_sheet_area(PS.results_pathname);
                
     
                clear MM
            end
         
            
        end
        
        function PS = GetFinalParameters(PS)
            ind=PS.goodar;
            PS.TaucMedMoy = dot([PS.Area(ind).Asc_time_median],[PS.Area(ind).PointNumb])/PS.TotalPoints;
            PS.TaucStdMoy = dot([PS.Area(ind).Asc_time_std],[PS.Area(ind).PointNumb])/PS.TotalPoints;
            PS.TaurMedMoy = dot([PS.Area(ind).Decay_time_median],[PS.Area(ind).PointNumb])/PS.TotalPoints;
            PS.TaurStdMoy = dot([PS.Area(ind).Decay_time_std],[PS.Area(ind).PointNumb])/PS.TotalPoints;
            PS.TaudMedMoy = dot([PS.Area(ind).Taud_median],[PS.Area(ind).PointNumb])/PS.TotalPoints;
            PS.TaudStdMoy = dot([PS.Area(ind).Taud_std],[PS.Area(ind).PointNumb])/PS.TotalPoints;
            PS.TaudBazMedMoy = dot([PS.Area(ind).Baz_taud_median],[PS.Area(ind).PointNumb])/PS.TotalPoints;
            PS.TaudBazStdMoy = dot([PS.Area(ind).Baz_taud_std],[PS.Area(ind).PointNumb])/PS.TotalPoints;
            PS.AireMedMoy = dot([PS.Area(ind).AUC_median],[PS.Area(ind).PointNumb])/PS.TotalPoints;
            PS.AireStdMoy = dot([PS.Area(ind).AUC_std],[PS.Area(ind).PointNumb])/PS.TotalPoints;
            PS.pcMedMoy = dot([PS.Area(ind).Asc_slope_median],[PS.Area(ind).PointNumb])/PS.TotalPoints;
            PS.pcStdMoy = dot([PS.Area(ind).Asc_slope_std],[PS.Area(ind).PointNumb])/PS.TotalPoints;
            PS.prMedMoy = dot([PS.Area(ind).Decay_slope_median],[PS.Area(ind).PointNumb])/PS.TotalPoints;
            PS.prStdMoy = dot([PS.Area(ind).Decay_slope_std],[PS.Area(ind).PointNumb])/PS.TotalPoints;
            PS.hrMedMoy = dot([PS.Area(ind).Amp_decay_median],[PS.Area(ind).PointNumb])/PS.TotalPoints;
            PS.hrStdMoy = dot([PS.Area(ind).Amp_decay_std],[PS.Area(ind).PointNumb])/PS.TotalPoints;
            PS.hcMedMoy = dot([PS.Area(ind).Amp_asc_median],[PS.Area(ind).PointNumb])/PS.TotalPoints;
            PS.hcStdMoy = dot([PS.Area(ind).Amp_asc_std],[PS.Area(ind).PointNumb])/PS.TotalPoints;
            PS.MMedMoy = dot([PS.Area(ind).Maxima_median],[PS.Area(ind).PointNumb])/PS.TotalPoints;
            PS.MStdMoy = dot([PS.Area(ind).Maxima_std],[PS.Area(ind).PointNumb])/PS.TotalPoints;
            PS.minMedMoy = dot([PS.Area(ind).Minima_median],[PS.Area(ind).PointNumb])/PS.TotalPoints;
            PS.minStdMoy = dot([PS.Area(ind).Minima_std],[PS.Area(ind).PointNumb])/PS.TotalPoints;

%             if (isfield(PS.Area(1),'MedRest'))
%                 PS.RestMedMoy = dot([PS.Area(:).MedRest],[PS.Area(:).PointNumb])/PS.TotalPoints;
%                 PS.RestStdMoy = dot([PS.Area(:).StdRest],[PS.Area(:).PointNumb])/PS.TotalPoints;
%             else
%                 PS.RestMedMoy =NaN;
%                 PS.RestStdMoy =NaN;
%             end
        end
        
        
    end
    
end
