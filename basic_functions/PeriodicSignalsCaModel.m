classdef PeriodicSignalsCaModel < handle
    % properties (SetAccess = private)
    properties
        FileName
        PathName
        BW
        boundary
        PixelSize
        SamplingFrequency
        framerate
        z
        R
        numim
        vecr
        vecmth
        intensr
        intensth
        Box=50
        vect
        pk1
        pk2
        Intpk1
        Intpk2
        phase1
        phase2
        F
        fvec
        Global1
        Global2
        T
        Area
        Signal
        Homogeneity
        TotalPoints
        MMedMoy
        MStdMoy
        mmvgMedMoy
        mmvgStdMoy
        mmvdMedMoy
        mmvdStdMoy
        RestMedMoy
        RestStdMoy
        Thampl
        Tauc
        Taud
        Taur
        Taus
        TausBaz
        TaudBaz
        TauBaz
        Aire
        pc
        pr
        hr
        hc
        mmvg
        mmgd
        mmvd
        mmgg
        Rest
        Period
        Time
        M
        xmc
        xMc
        xmr
        xMr
        locc
        locr
        MedTauc
        StdTauc
        MedTaur
        StdTaur
        MedTaud
        StdTaud
        MedTaus
        StdTaus
        MedTausBaz
        StdTausBaz
        MedTaudBaz
        StdTaudBaz
        MedAire
        StdAire
        Medpc
        Stdpc
        Medpr
        Stdpr
        Medhc
        Stdhc
        Medhr
        Stdhr
        MedM
        StdM
        Medmmvg
        Stdmmvg
        Medmmvd
        Stdmmvd
        MedRest
        StdRest
    end
    methods
        function PS = PeriodicSignalsCaModel(varargin)
            p = inputParser;
            addRequired(p, 'FileName');
            addRequired(p, 'Mask');
            addRequired(p,'PathName')
            addOptional(p, 'PixelSize', 1);
            addOptional(p, 'SamplingFrequency', 1);
            parse(p, varargin{1:end});
            PS.FileName = p.Results.FileName;
            PS.BW=p.Results.Mask;
            PS.PathName=p.Results.PathName;
            PS.PixelSize = p.Results.PixelSize;
            PS.SamplingFrequency=p.Results.SamplingFrequency;
            PS.framerate = (PS.SamplingFrequency)/1000; %fps in ms
        end
        
        function PS = determineCenter(PS)
            [B,~] = bwboundaries(PS.BW,'noholes');
            
            boundary=[];
            for k = 1:length(B)
                boundary = cat(1,boundary,B{k});
                
            end
            PS.boundary=boundary;            
            [PS.z,PS.R]= fitcircle(boundary);

            
        end
        
        function PS=getinfo(PS)
            info=imfinfo([PS.PathName,PS.FileName]);
            PS.numim=length(info);
            PS.PixelSize=1/info(1).XResolution;%pixel per µm
            infomulti=info(1).ImageDescription;
            posch1=strfind(infomulti,'finterval=');
            posch2=strfind(infomulti,'loop');
            posch1=posch1+10;
            posch2=posch2-1;
    
            dt=str2num(infomulti(posch1:posch2));
            if (~isempty(dt))
            PS.SamplingFrequency=1/dt;
            else
            posch2=strfind(infomulti,'fps');            
            posch2=posch2-1;
            dt=str2num(infomulti(posch1:posch2));
             PS.SamplingFrequency=1/dt;
            end
            PS.framerate = (PS.SamplingFrequency)/1000; %fps in ms
            
        end
        
        
        function PS = MeasureIntensity(PS,varargin)
            if ~isempty(varargin)
                PS.Box=varargin{1};
            end
            vecr=(0:5:PS.R);
            PS.vect=(0:PS.numim-1)/PS.SamplingFrequency;
            dtheta=PS.Box/2/PS.R;
            vecmth=0:dtheta:2*pi;
            
            intens=nan(length(vecr),length(vecmth));
            
            PS.vecr=vecr;
            PS.vecmth=vecmth;
            
            for tt=1:PS.numim
                
                im=double(imread([PS.PathName,PS.FileName],tt));
                Global2(tt)=nanmean(im(PS.BW(:)));
                cc=1;
                for r=vecr
                    
                    
                    if r==0
                        posrect=[PS.z(2)-PS.Box/2 PS.z(1)-PS.Box/2 PS.Box PS.Box];
                        imc=imcrop(im,posrect);
                        
                        
                        intens(cc,1)=sum(imc(:));
                        cc=cc+1;
                        
                    else
                        dtheta=PS.Box/r;
                        cc2=1;
                        for theta=0:dtheta:2*pi
                            cbox=[PS.z(2)+r*cos(theta) PS.z(1)+r*sin(theta)];
                            posrect=[cbox(1)-PS.Box/2 cbox(2)-PS.Box/2 PS.Box PS.Box];
                            
                            imc=imcrop(im,posrect);
                            
                            if size(imc)==[PS.Box+1 PS.Box+1]
                                
                                intens(cc,cc2)=sum(imc(:));
                                
                            end
                            cc2=cc2+1;
                        end
                        
                        cc=cc+1;
                    end
                    intensr(:,tt)=squeeze(nanmean(intens,2));

                    intensth(:,tt)=squeeze(nanmean(intens,1));
                end
            end
            
            for uu=1:size(intensr,1)
                intensr(uu,:)=detrend(intensr(uu,:),2);
            end
            M=max(intensr,[],2);
            for ii=1:length(M)
                intensr(ii,:)=intensr(ii,:)/M(ii);
            end

            for uu=1:size(intensth,1)
                intensth(uu,:)=detrend(intensth(uu,:),2);
            end
            M=max(intensth,[],2);
            for ii=1:length(M)
                intensth(ii,:)=intensth(ii,:)/M(ii);
            end
            
            PS.intensth=resize(intensth, [sum(~isnan(intensth), 'all')/size(intensth,2) size(intensth,2)]);
            PS.intensr=intensr;
            PS.Global2=Global2;
            PS.Global1=nanmean(PS.intensr,1);
%             figure
%             hold on
%             plot(PS.Global1)
%             plot(detrend(PS.Global2))
            
            
        end
       
        
        
        function PS = FourierAnalysis(PS, varargin)
            %find both big peaks
            
            
            %% find position of pk2
            F=abs(fft(PS.intensr(1,:))).^2;
            F(1)=0;F(end)=0;
            F=F/sum(F(:));
            PS.fvec=[0:1/length(F):(length(F)-1)/length(F)]*PS.SamplingFrequency;
            [pks,locs]=findpeaks(F(1:round(end/10)),'sort','descend');
            PS.pk1=locs(1);
            PS.Period=1/PS.fvec(locs(1));
            PS.F=F;
            
            
            [~,pos]=max(abs(F(50:round(length(F)/2))));
            pk2=pos+49;
            PS.pk2=pk2;
            
            PS.Time=[1:size(PS.intensr,2)]/PS.SamplingFrequency;
            
            %%
            M=size(PS.intensr,1);
            %% fourier calculation
            for ii=1:M
                F=fft(PS.intensr(ii,:));
                Fa=abs(F).^2;
                Fa(1)=0;
                Fa(end)=0;
                Fa=Fa/sum(Fa);
                phase1(ii)=atan2(imag(F(PS.pk1)),real(F(PS.pk1)));
                phase2(ii)=atan2(imag(F(PS.pk2)),real(F(PS.pk2)));
                Intpk2(ii)=abs(Fa(PS.pk2));
                Intpk1(ii)=abs(Fa(PS.pk1));
                
            end
            PS.phase1=phase1;
            PS.phase2=phase2;
            PS.Intpk1=Intpk1;
            PS.Intpk2=Intpk2;
            
            
            
            
            
            
        end
        
        
        
        
        
      function PS=GetCaractPeak(PS,varargin) 
          
          
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
             MM(:,1)=PS.Time;
             MM(:,2)=PS.Global1;
            PK=AnalysisPeaks(MM,'list_param_name',list_param_name,'list_calc', list_calc,'PixelSize',PS.PixelSize,'smoothness',sl,'prop',prop);
            PK.matrix_filtered_fluorescences(:,1)=MM(:,2);
            PK.CalculateParameters(1);           
            PS.Signal=PK.Save_in_struct(PS.Signal);
          
      end
        
    end
    
end
