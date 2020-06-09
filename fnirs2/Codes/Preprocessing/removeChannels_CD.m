% function [mlact r] = removeChannels(d,t,SD,fCard,threshCard,fNoise,threshNoise)
%
% This function will remove channels based on cardiac peak and noise floor.
% Either criteria must be met for all wavelengths for a given source detector pair.
%
% INPUTS
% d - intensity data (#time points x #channels)
% t - time vector (#time points x 1)
% SD - the SD structure
% fCard - frequency window containing cardiac peak [fmin fmax].
%         It is important that this window be large but the peak must be
%         the cardiac peak.
% threshCard - A number??? below which he channel is rejected.
% fNoise - frqeuency window for estimating the noise floor.
% threshNoise - ???
%
% OUTPUTS
% mlact - the active measurement list contains a 1 for an active channel
%         and a 0 for a channel that did not pass the criteria. (#Measurements x 1)
% r - ???
%
% DEPENDENCIES
% None
%
% TO DO
% This presently assumes that MeasList is ordered by first wavelength and
% then second. Handle SD.
% generalize to N wavelengths (don't assume N=2)


function [mlact r mnsn] = removeChannels_CD(d,t,SD,fCard,threshCard,fNoise,threshNoise,allchans,figs)

% frequency window where to look for heartbeat
%fwin = [100/60 250/60];
%fwnoise = [200/60 240/60];
% treshold for r, remove channels with lower r
%rtr = 4.0; % Treshold for heartbeat detection
%rtn = 0.15; % Treshold for noise detection

nCh = size(d,2)/2; %number of channels per wavelength
[d]=hmrIntensity2OD(d);% not sure why but i believe code wants dod not d, look into converting and resaving .nirs with dod not d

fs = 1/(t(2)-t(1));

np = 2^(nextpow2(size(d(1:end,:),1))-1); % define number of points in fft
fftdat = abs(fft(d(1:end,:),np)); % absolute power spectra

f = fs/2 * linspace(0,1,np/2); % frequency vector

fistart = find(f > fCard(1), 1); % start of cardiac frequency 
fistop = find(f > fCard(2), 1); % end of cardiac frequency

fnstart = find(f > fNoise(1), 1); % start of frequency noise
fnstop = find(f > fNoise(2), 1); % end of frequency noise

if isempty(fnstop)
    fnstop = length(f);
end

if isempty(fnstart)
fnstart = fistop;
end

% calculate first score ra
mns = mean(fftdat(fistart:fistop,:)); % take mean of each channel power spectra across cardiac range 
stds = std(fftdat(fistart:fistop,:)); % std value of for each channel power spectra
%[maxs maxidx] = max(fftdat(fistart:fistop,:));
[fftsort fftsidx] = sort(fftdat(fistart:fistop,:),'descend'); % sorted fft output with vector describing rearrangement
%maxs = (mean((fftsort(1:3,:)-ones(3,size(fftsort,2))*diag(mns)))).^0.5;
maxs = mean(fftsort(1:3,:)); %mean of top 3 fft values
ra = (maxs-mns)./stds; % subtract cardiac floor from peak, then divide by cardiac window std. High number for more defined peak, accounting for cardiac window mean and std. Cardiac Peak power value

% fix to remove channels with clippling
ra = ra .* (fftsort(1,:)<8); 

maxidx = fftsidx(1:3,:) + fistart -1; % find index of 3 peaks for each channel
stdidx = std(reshape(fftsidx,[],1)); % creates vector of all channel index arrays and then takes std. Gives one value that represents the standard deviation of peak idx distribution
%mnidx = mean(mean(maxidx));
tmp = mean(maxidx); % take mean of 3 peak idx
mnidx = ra(1:nCh)*tmp(1:nCh)'/sum(ra(1:nCh));   % SPECIFIES ONE WAVELENGTH. Matrix multiply Cardiac Peak power score & peak idx. 
mnidxdev = mean(abs(fftsidx(1:3,:) + fistart -1 - mnidx)); % basically figures out how far each channels peaks deviate from the norm.

%maxidx = maxidx + fistart -1;
rb = 0.5 * mnidxdev./stdidx; % normalized deviation of where each channels frequency peaks are

mnsn = mean(fftdat(fnstart:fnstop,:)); %mean of noise spectra
stdn = std(fftdat(fnstart:fnstop,:)); % std of noise spectra


r = ra-rb; % peak power score - deviation of peak frequencies
%ch = find(r < threshCard);
mlact = ~(mnsn>threshNoise | r <= threshCard); % exclude if mean is above threshNoise and r below cardiac peak thresh
% mlact(1:nCh) =(r(1:nCh) >= threshCard);
% mlact(nCh+1:2*nCh) = (stdn(nCh+1:2*nCh)<threshNoise | r(nCh+1:2*nCh) >= threshCard2);

% make sure both wavelength of one SD pair are deactivated
% THIS ASSUMES THAT MeasList IS ORDERED BY FIRST THEN SECOND WAVELENGTH
nSD = length(mlact)/2;
mlact(1:nSD) = mlact(1:nSD) & mlact(nSD+1:end);
mlact(nSD+1:end) = mlact(1:nSD) & mlact(nSD+1:end);

if figs == 1
%%% Chris Edit
Ch = size(fftdat,2);
% % Define axes stretched past the default positions
m = 10; % 9 rows
n = ceil(Ch/m); % # of columns determined by # of channels to plot
l = m*n; % total # of plots (SNR counts as entire row of plots)
Rect = [0.05, 0.05, 0.93, 0.9];
AxisPos = moPlotPos(n,m, Rect);


% Plot cardiac frequency spectra
% Figure handle
h2=figure('name',['Cardiac Frequency'],...
    'NumberTitle','off',...
    'Color','White',...
    'units','normalized',...
    'outerposition',[0 0 1 1],...
    'DefaultaxesPosition',Rect);
for idx=1:size(fftdat,2)
    subplot(m,n,idx);
    if mlact(idx) == 1
        plot(f(fistart:fistop),fftdat(fistart:fistop,idx),f(maxidx(:,idx)),fftdat(maxidx(:,idx),idx),'o'); % plot spectra in cardiac range and circle top 3 peaks
    else
        plot(f(fistart:fistop),fftdat(fistart:fistop,idx),'r',f(maxidx(:,idx)),fftdat(maxidx(:,idx),idx),'o'); % if channel is cancelled make it red
    end
    set(gca,'XLim',fCard);
    if ~isnan(mnidx)
        line([f(round(mnidx)) f(round(mnidx))],get(gca,'YLim'),'Color','red','LineStyle','--');
        if idx <= nCh
        title(['CH' num2str(idx) ' r=' num2str(r(idx),2) ' mn=' num2str(mnsn(idx),2)]);
        else 
        title(['CH' num2str(idx-nCh) ' r=' num2str(r(idx),2) ' mn=' num2str(mnsn(idx),2)]);
        end
    end
    % cross out previously pruned channels 
    if ~ismember(idx,allchans)
        xl=get(gca,'XLim');yl=get(gca,'YLim');
        line(xl,yl,'Color','red','LineStyle','--');
    end
%         if idx <= nCh
%             title(['r=' num2str(r(idx),2) ' ra=' num2str(ra(idx),2) ' rb=' num2str(rb(idx),2)]);
%         else
%             title(['r=' num2str(r(idx),2) ' st=' num2str(stdn(idx),2) ' mn=' num2str(mnsn(idx),2)]);
%         end
%     if idx<=idx/2
%         title(num2str(idx));
%     else
%         title(num2str(idx-(idx/2)));
%     end
end
end


    
 %%% END EDIT

