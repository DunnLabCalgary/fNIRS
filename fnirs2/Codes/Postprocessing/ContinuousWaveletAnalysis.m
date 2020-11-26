function [X] = ContinuousWaveletAnalysis(FZband,Proc, RegionalData, Hemo, stimRects, procFile, sub, WaveletFig)
% profile on

% This function is used to calculate pairwise phase coherence and power
% time-series. It also calculates pairwise frequency varying vectors for
% each parameter.

% Input: 
%   Data = info from UI
%   t = time vector
%   s = stim vector
%   SD = source detector configuration struct
%   dc = preprocessed concentration data (concentration data x hemo parameter x channel)

% Output:
%   X = struct containing all the wavelet parameters calculated

%% Initialize
    tic
    t = Proc.t;
    fs = 1/(t(2)-t(1));
    SD = Proc.SD;
    f1 = FZband(1);
    f2 = FZband(2);
    nCh = 1:size(Proc.dc,3); % channel vector
    GoodCh = Proc.QC.dqChans(Proc.QC.dqChans<=max(nCh)); % Good channels after preprocessing... need to make it more universal, code it in the SD.MeasListAct format
    %Proc.QC.goodChannels;%find(find(SD.MeasListAct)<=max(nCh));
    
    %Construct stim times 
    si = {};
    for idx = 1:length(stimRects)
        r = stimRects(idx);
        lb = r{1}(1);
        ub = lb + r{1}(3);
        s = t>lb & t<ub;
        
        si{idx} = s;
    end
    
    rc = RegionalData.rc;
    
    if rc == 1
        wlabel = 'r1 intra coh';
       
        ChannelsX = intersect(GoodCh, RegionalData.R1Channels);
        ChannelsY = ChannelsX;
    elseif rc == 2
        wlabel = 'r1 intra coh';
       
        ChannelsX = intersect(GoodCh, RegionalData.R1Channels);
        ChannelsY = intersect(GoodCh, RegionalData.R2Channels);
    elseif rc == 3
        wlabel = 'r2 intra coh';
     
        ChannelsX = intersect(GoodCh, RegionalData.R2Channels);
        ChannelsY = ChannelsX;
    end
   
    X = regionalComparisonWaveletAnalysis(ChannelsX, ChannelsY, Proc,si, f1, f2, fs, Hemo, 0, procFile, sub, wlabel, WaveletFig);

% profile viewer
% profile off
end

function [X] = regionalComparisonWaveletAnalysis(ChannelsX, ChannelsY, Proc,si, f1, f2, fs, Hemo, detailLog, procFile, sub, wlabel, WaveletFig)

%% Pairwise Wavelet Analysis
for hb = 1:numel(Hemo)
    Hb = Hemo{hb};
    
    w = waitbar(0 ,['Processing ' Hb ' Wavelet Analysis' '-' sub '-' wlabel]);
    if Hb == "StO2"
        d = squeeze(Proc.dc(:,1,:)./Proc.dc(:,3,:)); % sto2 == hbo./thb
        %save the new proc var (preprocessed .mat file)  with sto2 in it.
        Proc.dc(:,4,:) = d;
        
        save([procFile filesep 'Preprocessing' filesep 'Processed.mat'],'Proc');
        
     
    
    else
        d = squeeze(Proc.dc(:,hb,:));
    end
    
    %% 3 for testing
    for idx = 1:length(ChannelsX) %Loop for Input Ch x
        ChX = ChannelsX(idx);
        
            waitbar(idx/length(ChannelsX),w,['Processing ' Hb ' Wavelet Analysis']);
         
            x = d(:,ChX); % Input data from channel (x)
            
            % wt = continuous wavelet transform power for channel (x),  
            % fzt = frequency in Hz
            % coiwt = cone of influence 
            [wt,fzt,~] = cwt(x,fs); % NumScalesToSmooth..?
            wt = abs(wt);
            
            %fi = binary vector with one in the places where the cwt freq
            %(fzt) is between f1 and f2 (specified in the gui).
            fi = fzt>f1 & fzt<f2;
            wt = squeeze(mean(wt(fi,:))); % collapses values across fz band of interest, fz-averaged time-varying vector of power
            
            if idx == 1 
                tp.wt = zeros(max(max(ChannelsX),max(ChannelsY)), length(wt));
            end 
            
            %store power vector values for later 
            tp.wt(ChX,:) = wt';  
            
            X.fzwt = fzt(fi);  %save the fz vector of interest 
            X.fzwtf = fzt; %save the entire fz vector returned from the cwt
           
            for ii = 1:length(ChannelsY) % Loop for input channel (y) to do cross-spectrum and coherence with (X)
               
                ChY = ChannelsY(ii);

                if ChX ~= ChY  % dont do coherence when (x) = (y)
                    
                    if detailLog
                        disp('comparing channels')
                        disp(ChX);
                        disp(ChY);
                    end
                
                    y = d(:,ChY); % Input data from channel (y)
                    
                     if tp.wt(ChY,1) == 0 %dont redo calc if the regions have overlapping channels
                        %% Power y
                        [wt,fzt,~] = cwt(y,fs); % NumScalesToSmooth..?
                        wt = abs(wt);
                        fi = fzt>f1 & fzt<f2;
                        wt = squeeze(mean(wt(fi,:))); % collapses values across fz band of interest, fz-averaged time-varying vector of power
                        %store power vector values for later 
                        tp.wt(ChY,:) = wt';  
                     end 
                    
                    %split up time for fz varying analysis;
%                     h = figure;
%                     b = axes;
                    for stim = 1:length(si)
                        
                        %data points during the stim
                        xs = x(si{stim});
                        ys = y(si{stim});
                        
                        if isempty(xs) || isempty(ys)
                            disp('not data points during stim')
                            return
                        end
                        
                        [wcoh,wcs,fz, coi] = wcoherence(xs,ys,fs);
                        
                        fvCoh = squeeze(nanmean(wcoh,2)); %fz varying coherence vector, collapse across time
                        
                        %insert the matrix into ch wise cell array
                        tp.wcohf(stim).stim(ChX,ChY,:) = fvCoh;
                       
                        if ~isfield(tp, 'wcohfzRange')
                            tp.wcohfzRange(stim).fz(:) = fz;
                        elseif stim > length(tp.wcohfzRange)
                            tp.wcohfzRange(stim).fz(:) = fz;
                        end
                        
                    end

%                     hold off;
                    
                    %over entire series for time varying coherence
                    [wcoh,wcs,fz,coi] = wcoherence(x,y,fs,'NumScalesToSmooth',6); %~, ~ -> wcs, coi, 
                    %% plot wavelet? 
                    if WaveletFig
                        helperPlotCoherence(wcoh,Proc.t,fz,coi,'Seconds','Hz', ChX, ChY, sub ,Hb);
                    end
                    
                    fi = fz>f1 & fz<f2; % same as other fi take fz band of interest
                    wcoh = squeeze(nanmean(wcoh(fi,:))); % collapses values across fz band of interest, fz-averaged time-varying vector
                    tp.wcoh(ChX,ChY,:) = wcoh'; % store phase vector in ch wise matrix
                    
                    % phase angle
                    wcs = angle(wcs); % wcs is a complex grid of size fz * #data points 
                    wcs(wcoh<0.5) = nan;

                    phase = squeeze(nanmean(wcs(fi,:))); % collapses values across fz band of interest, fz-averaged time-varying vector 
                    tp.phase(ChX,ChY,:) = phase'; % store phase vector in ch wise matrix 
%                      
                    X.fzwcoh = fz(fi); % save fz band vector of interest
                    
                end
            end
        
    end
    close
    %clean up tp values and save them in X for returning them
    %phase
    tp.phase(tp.phase==0)=NaN;
    X.phase.(Hb) = single(tp.phase);
    %tv coh
    tp.wcoh(tp.wcoh==0)=NaN;
    X.wcoh.(Hb) = single(tp.wcoh);
    %power
    tp.wt(tp.wt==0)=NaN;
    X.wt.(Hb) = single(tp.wt);
    %fv coh
    X.wcohf.(Hb) = tp.wcohf;
    X.wcohfzRange = tp.wcohfzRange;
  
%     tp.wtf(tp.wtf==0)=NaN;
%     X.wtf.(Hb) = single(tp.wtf);
%     tp.wcohf(tp.wcohf==0)=NaN;
%     X.wcohf.(Hb) = single(tp.wcohf);
%     tp.phasef(tp.phasef==0)=NaN;
%     X.phasef.(Hb) = single(tp.phasef);
    clear tp;
end
    X.t = Proc.t;
    X.s = Proc.s;
    X.ChX = ChannelsX;
    X.ChY= ChannelsY;
    X.RunTime = toc;
    X.DateTime = datetime;
end

