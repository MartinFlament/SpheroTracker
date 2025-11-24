classdef FourierAnalysisModel < handle
    properties (SetAccess = private)
        Signal
        SamplingFrequency
        X
        Fourier
        Frequency
        Period
        Amplitude
        Phase
        PhaseShift
    end
    methods
        function FA = FourierAnalysisModel(varargin)
            p = inputParser;
            addRequired(p, 'Signal', @isvector);
            addParameter(p, 'SamplingFrequency', 1, @isnumeric);
            parse(p, varargin{:});
            
            if iscolumn(p.Results.Signal)
                FA.Signal = p.Results.Signal';
            else
                FA.Signal = p.Results.Signal;
            end
            FA.SamplingFrequency = p.Results.SamplingFrequency;
        end
        function FA = performFft(FA,varargin)
            
            SignalLength = length(FA.Signal);
            Window = (cos(linspace(-pi/2, pi/2, SignalLength))).^2;
            N = 2^(nextpow2(SignalLength)+3);
            
            df = FA.SamplingFrequency / N; %frequency resolution
            sampleIndex = 0 : N - 1; %ordered index for FFT plot
            
  
            
            FA.X = 1/N * fft(Window.* FA.Signal, N);
            FA.Fourier = abs(FA.X);
            FA.Frequency = sampleIndex * df;
            if ~isempty(varargin)
                PP=varargin{1}*FA.SamplingFrequency;
                kkp=round((N/PP));
                larg=40;
                range=max(1,kkp-larg):min(kkp+larg,length(FA.Fourier));
            else
                larg=0;
                range=1:round(length(FA.Fourier)/2);
            end
        
            
            [pks, locs] = findpeaks(FA.Fourier(range), 'SortStr', 'descend');
             %  figure
             % hold on
             % plot(FA.Fourier(range))
             % plot(locs(1),pks(1),'r+')
             % 
             %  figure
             % hold on
             % plot(FA.X(range))
             % %plot(locs(1),pks(1),'r+')
             % disp(size(FA.X))
             % pause

          
         
            % if ~isempty(locs)
            % % if locs(1)<100 
            % %     pks(1)=[];
            % %     locs(1)=[];
            % % end
            % end
            if~isempty(varargin)
                if ~isempty(pks)
                    %kkp=kkp-larg+locs(1)-1;
                    FA.Period = N/kkp/FA.SamplingFrequency;
                    FA.PhaseShift=atan2(imag(FA.X(kkp)),real(FA.X(kkp)));
                    FA.Amplitude=pks(1);
%                                 figure
%                                 hold on
%                                 plot(FA.Fourier(range))
%                                  plot(locs(1),pks(1),'+')
%                                  plot(locs(2),pks(2),'+')
%                                 pause
                else
                    FA.Period = NaN;
                    FA.PhaseShift=atan2(imag(FA.X(kkp)),real(FA.X(kkp)));
                    FA.Amplitude=FA.Fourier(kkp);
                end
            else
                if ~isempty(pks)

                FA.Period=1/FA.Frequency(locs(1));
                FA.PhaseShift=atan2(imag(FA.X(locs(1))),real(FA.X(locs(1))));
                FA.Amplitude = pks(1);
                else
                FA.Period=NaN;
                FA.PhaseShift=NaN;
                FA.Amplitude = NaN;
                end
            end
            
            %    FA.PhaseShift=angle(FA.X(locs(1)));
            
        end
        
        
        
    end
end