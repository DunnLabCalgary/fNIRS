function [Coh,fyy]=fftCoh2(y,SamplingRate,lowfz,highfz);
%Calculated the Coherence of the Ref signal with everything else in the
%frequency domain lowfz to highfz

new_band='y';
Itt=0;
iii=0;
fm=lowfz;
fM=highfz;

for idx = 1:size(y,2)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Calculation of coherence %%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    K=8; %number of sections %8 segments
    D=0.5; % overlap %50% overlap
    x=y(:,idx); % REF signal
    L=round(length(x)/(1+(K-1)*(1-D)));  %Use REF for lengths  / 4.5
    BF=Bfilter(ceil(L/2));
    if size(BF,2)~=L
        BF = BF(1:L);
    end
    U=hamming(L)'*hamming(L); %yey hamming window
    while 2^iii<L
        iii=iii+1;
    end
    nFFT=2^(iii);
    %fft(x) computes discrete fourier transform having as many samples
    %as the original signal x whereas fft(x,n) will have n points in dft. Usually n is taken as integral power of 2 to fasten the computation. In order to recover original signal from fft, number of samples in fft must be greater than or equal to number of samples in original signal . As n increases resolution in frequency domain increases i.e the distance between two consecutive frequency bins decreases .
    
    
    
    for i=1:size(y,2)   % number of channels... exclude?
        Sxy=0;
        Sxx=0;
        Syy=0;
        Itr=0;
        for ii=1:round(D*L):length(x)-round(L)-1 % counter for start of each Ref segment
            Xx=x(ii:ii+L-1)'; % pull x vector of segment
            Xx=(Xx-min(Xx))/(max(Xx)-min(Xx)); %normalize the time segment
            slope=(Xx(length(Xx))-Xx(1))/length(Xx); %one number
            bias=Xx(1);
            t=1:length(Xx);
            Xx=Xx-slope*t-bias;%take trend away
            Xx=Xx-mean(Xx); %take mean away
            Yy=y(ii:ii+L-1,i)';  %USING segment of Sigs -> temporal segments!
            %(but these are only 7... 1:664, 332:996,
            %664,996,1328,1660,1992,2324 (7*332 which is over ii
            %threshold 2321)
            %frequency cohrence in the certain timeinterval
            Yy=(Yy-min(Yy))/(max(Yy)-min(Yy));
            slope=(Yy(length(Yy))-Yy(1))/length(Yy);
            bias=Yy(1);
            Yy=Yy-slope*t-bias;
            Yy=Yy-mean(Yy);
            Xh=BF.*Xx; %Bfilter =smoothes the edges of the temporal segment
            Yh=BF.*Yy;
            Fx=fft(Xh,nFFT); %frequency domain of Ref segment 1024... why is that?
            Fy=fft(Yh,nFFT);%frequency domain of normal SD pair segment
            
            %                 %To see the frequency domain signal
%                                             f = SamplingRate*(0:(nFFT/2))/nFFT;
%                                             P = abs(fyy(1,:)/nFFT);
%                                             plot(f,P(1:nFFT/2+1))
%                                             title('Gaussian Pulse in Frequency Domain')
%                                             xlabel('Frequency (f)')
%                                             ylabel('|P(f)|')
            %                 %
            
            Sxy=Sxy+Fx.*conj(Fy)/(U); %outputs the multiplication between the frequencies ?
            %accentuates the higher vs lower ones (lowers more lower, higher more higher)
            %and adds each timesegments Sxy on top of it!
            Sxx=Sxx+Fx.*conj(Fx)/(U);
            Syy=Syy+Fy.*conj(Fy)/(U);
            Itr=Itr+1;
            
        end
        fxy(i,:)=Sxy/Itr;  %%%%% Calculating mean %%%%%%%%% This is the mean coherence over the segments over all frequencies
        fxx(i,:)=Sxx/Itr;  %%%%% Calculating mean %%%%%%%%%
        fyy(i,:)=Syy/Itr;  %%%%% Calculating mean %%%%%%%%%
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Changing frequency to sample %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    FF=size(fxy,2)/2;
    ffm=round(fm*FF/(SamplingRate/2));
    ffM=round(fM*FF/(SamplingRate/2));
    if ffm==0
        ffm=1;
    end
    if ffM==0
        ffM=1;
    end
    for i=1:size(fxx,1)
        fxx1=mean(fxx(i,ffm:ffM)); %output is 1 number. coherence of specific frequency over the segment average coherence
        fyy1=mean(fyy(i,ffm:ffM));
        fxy1=mean(fxy(i,ffm:ffM));
        Coh(idx,i)=abs(fxy1.^2)./(fxx1.*fyy1);
    end
    
    %Importnat output: Coh
    
end
end