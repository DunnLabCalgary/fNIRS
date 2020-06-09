function [chanList detailedInfo]=hmrDataQualityCheck_CD(X,modeFlags,dBrange,...
    sdDistRange,snrMax,tRange,poxChan)
% chanList = hmrDataQualityCheck(filenm,modeFlags,dBrange,sdDistRange,snrMax,tRange)
%   OR
% chanList = hmrDataQualityCheck
%   OR
% [chanList,detailedInfo] = hmrDataQualityCheck(filenm,modeFlags,dBrange,sdDistRange,snrMax,tRange)
%   OR
% chanList = hmrDataQualityCheck(filenm,modeFlags,dbRange,sdDistRange,snrMax,tRange,poxChannel) 
%   OR
% [chanList,detailedInfo] = hmrDataQualityCheck(filenm,modeFlags,dbRange,sdDistRange,snrMax,tRange,poxChannel) 
%
% This function will prune away any channels that do not satisfy a number
% of requirements.
%    1. dBrange - the power of the raw data on the channel must be within
%                 a range.
%    2. snrMax - the maximum SNR is below a value.
%    3. sdRange - the distance between source and detector must be within
%                 a range.
%    4. HbMatch - corresponding oxy and deoxy channels are pruned if the
%                corresponding channel is bad.
%    5. cardiacPeak - (optional) a cardiac peak must be present on good channels
%   
%INPUTS:
% IF YOU CALL THIS COMMAND WITH NO INPUT ARGUEMENTS, EVERYTHING WILL 
% REVERT TO DEFAULT. THE WHOLE DIRECTORY WILL BE READ.
%
% filenm - either a filename to test
%       default (if empty string): entire directory
% modeFlags - [dBrangeTF snrMaxTF sdRangeTF HbMatchTF cardiacPeakMODE motionTF, cardiacFreq, DataQuality]
%       default (if empty array): [1 1 1 1 0 0]
%     dBrangeTF = use dB range? (0 or 1)
%     snrMaxTF = use SNR max? (0 or 1)
%     sdRangeTF = use source-detector range? (0 or 1)
%     HbMatchTF = use only pairs of oxy/deoxy that are both good? (0 or 1)
%     cardiacPeakMODE = use only channels with cardiac peak?
%               0 - no
%               1 - if no poxChan is provided, use the removeChannels script,
%                based on sorted FFT method. If poxChan is given, then
%                use the hmrTestNIRSforHR script, based on peak FFT method 
%     motionTF = obtain motion artifact data (must use 1 file, not
%                directory)
%     cardiacFreq = plot cardiac Freq? (0 or 1)
%     DataQuality = plot Data quality? (0 or 1)
% dBrange - the min and max allowable dB value for a raw signal
%       default (if empty array): [80 125]
% sdDistRange - the min and max allowable separation between
%       source and detector. Use 0 and/or inf as needed to unbound one side.
%       default (if empty array): [0 3.5]
% snrMax - the maximum value permitted for snr = std.dev.(d)/mean(d)
%       default (if empty array): 5
% tRange - range of values to examine for artifacts (must set motionTF 1)
% poxChan - the aux channel number for pulseox
%
%OUTPUTS:
% chanList - a set of channels indices that correspond to d from the nirs
%   files which remain after the pruning.
% detailedInfo - if you only test 1 file, this optional output will give
%   you more detail about why channels were pruned
%
% CALLS:
% hmrPruneChannels.m, removeStimuli.m, removeChannels.m, drawProbeGood.m
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if ~exist('filenm'),        filenm='';      end;
if ~exist('modeFlags'),     modeFlags=[];   end;
if ~exist('dBrange'),       dBrange=[];     end;
if ~exist('sdDistRange'),   sdDistRange=[]; end;
if ~exist('snrMax'),        snrMax=[];      end;
if ~exist('tRange'),        tRange=[];      end;
if ~exist('poxChan'),       poxChan=[];     end;
runFULLdirectory=0;

%% ESTABLISH DEFAULTS
% if isempty(filenm),         runFULLdirectory=1;     end;
if isempty(modeFlags),      modeFlags=[1 1 1 1 0 0 1 1]; end; % Chris added flags 7 and 8 for figure choice
if isempty(dBrange),        dBrange=[80 125];       end;
if isempty(sdDistRange),    sdDistRange=[0 3.5];    end;
if isempty(snrMax),         snrMax=5;               end;
if (modeFlags(1)==0),       dBrange = [-inf inf];   end;
if (modeFlags(2)==0),       snrMax= inf;            end;
if (modeFlags(3)==0),       sdDistRange=[0 inf];    end;

detailedInfo = struct;

%% Compute allchans and min,max,median stats
    if ~exist('X.ml')
        X.ml = X.SD.MeasList;
    end
        
    if sum(size(X.SD.MeasList)==size(X.ml))<2 
        X.SD.MeasList = X.ml;
    end
    [chans,SNR,dists] = hmrPruneChannels_CD(X.d,X.SD,sdDistRange,modeFlags(4),dBrange,snrMax);
    fprintf('%d ',chans); fprintf('\n');
    
    % Take min/max/median of log transform for dB amplitude criteria
    allchans = chans;
    minD = 20*log10(min(X.d));
    maxD = 20*log10(max(X.d));
    medianD = 20*log10(median(X.d));

%% If needed, remove non cardiac artifact containing channels
if isempty(poxChan)
    % if no pulseox channel is given, use the Bernhard method
    % which uses sorted FFT to pick channels
    fCard = [0.7 1.8]; % Cardiac frequency range
    threshCard = 4; % Cardiac power threshold
    fNoise = [2.5 3.0000]; % Define noise frequencies
    threshNoise = 0.5; % Noise threshold
    [detailedInfo.chinv detailedInfo.cardiacR detailedInfo.NoiseFloor] = ...
        removeChannels_CD(X.d,X.t,X.SD,fCard,threshCard,fNoise,threshNoise,allchans,modeFlags(7)); %Chris added allchans and figs Feb.4 2019
    goodHRchans = find(detailedInfo.chinv);
else
    % if a pulseox channel is given, use the Daniel method
    % which matches peak FFT values between channels and pulseox
    [goodHR,goodHRpairs,suspectHR,HR] = hmrTestNIRSforHR(X.d,X.t(1:2),X.SD,X.aux(:,poxChan),0);
    goodHRchans = find(goodHR);
end
detailedInfo.cardiacGoodChans = goodHRchans;
if modeFlags(5)==1 % only use cardiac criteria if it is selected
    allchans = intersect(allchans,goodHRchans);
end


%% additional calculations
sd = size(X.d,2);
sd2 = sd/2;
chansTFeach = zeros(1,size(X.d,2));
chansTFeach(allchans) = 1;
temp = chansTFeach(1:sd2) & chansTFeach((sd2+1):end);
chansTFall = [temp temp];
listlb1 = getNirsList(X.SD,X.ml,X.SD.Lambda(1));
listlb2 = getNirsList(X.SD,X.ml,X.SD.Lambda(2));
try
    gains = X.systemInfo.gain;
catch
    gains = zeros(1,200);
end
try
    DarkNoise = X.systemInfo.DarkNoise;
catch
    DarkNoise = zeros(2,200);
end
m = length(listlb1);


%% Require channels to match or not?
if modeFlags(4)==0
    chanList = allchans;
else
    chanList = find(chansTFall);
end

%% Text output
fprintf('Individual good channels:\n');
fprintf('%d ',allchans); fprintf('\n');
fprintf('Sets of good...\n');
fprintf('%d ',find(chansTFall)); fprintf('\n');
fprintf('Halfset...\n');
fprintf('%d ',find(temp)); fprintf('\n');

%% Plots
if modeFlags(8)==1 % 1 = YES for plot
Rect = [0.05, 0.05, 0.93, 0.9];
h1=figure('name',['Check Quality'],...
    'NumberTitle','off',...
    'Color','White',...
    'units','normalized',...
    'outerposition',[0 0 1 1],...
    'DefaultaxesPosition',Rect);
%Median min max plot
subplot(5,2,1);
hold on
for a=1:m
    plot((a)*[1 1],[minD(listlb1(a)) maxD(listlb1(a))],'b');
    plot((a),medianD(listlb1(a)),['xb']);
    plot((a)*[1 1],[minD(listlb2(a)) maxD(listlb2(a))],'r');
    plot((a),medianD(listlb2(a)),'xr');
end
plot([1 m],dBrange(1)*[1 1],'--k');
plot([1 m],dBrange(2)*[1 1],'--k');
xlim([1 m]);
ylim([dBrange(1)-15 dBrange(2)+15]);
hold off
title('Median,min,max')
xlabel('Channel #');
ylabel('Amplitude (dB)');
% legend(num2str(X.SD.Lambda))

%Gains plot
subplot(5,2,9);
if size(gains) == [X.SD.nSrcs X.SD.nDets] % indicates NIRx gains matrix rather than TechEn array
    idx = X.ml(listlb1,[1 2]);
    for ct = 1:size(idx,1)
        gplot(ct) = gains(idx(ct,1),idx(ct,2));
    end
%     thresh = ?
else
    gplot = gains(X.ml(listlb1,2));
%     thresh = ?
end
bar(gplot)
% plot thresh..?
xlim([1 m]);
title('Gains')
xlabel('Channel #')

%DarkNoise plot - from NIRX noise output with no light
subplot(5,2,7);
hold on
for a=1:X.SD.nDets
plot((a),DarkNoise(1,a),['xb']);
plot((a),DarkNoise(2,a),['xr']);
end
% threshold...?
legend(num2str(X.SD.Lambda))
title('Dark Noise')
xlabel('Detector #');
ylabel('Amplitude (mV)');
xlim([1 X.SD.nDets])

%SNR/CV plot
subplot(5,2,2);
hold on
for a=1:m
plot((a),SNR(a),['xb']);
plot((a),SNR(a+m),['xr']);
end
% plot([1 m],snrMax*[1 1],'-k'); threshold
% ylim([0 snrMax]);
legend(num2str(X.SD.Lambda))
title('Coefficient of Variation')
xlabel('Channel #');
xlim([1 m])

% Cardiac Pulse Power plot
subplot(5,2,3);
hold on
for a=1:m
        plot((a),detailedInfo.cardiacR(a),['xr']);
        plot((a),detailedInfo.cardiacR(a+m),['xb']);
end
plot([1 m],threshCard(1)*[1 1],'--k');
ylim([0 max(detailedInfo.cardiacR)+1]);
title('Cardiac Pulse Power Criteria')
xlabel('Channel #');
ylabel('Pulse Power index');
xlim([1 m])

% Cardiac Noise plot
subplot(5,2,5);
hold on
for a=1:m
        plot((a),detailedInfo.NoiseFloor(a),['xr']);
        plot((a),detailedInfo.NoiseFloor(a+m),['xb']);
end
plot([1 m],threshNoise(1)*[1 1],'--k');
mNF = max(detailedInfo.NoiseFloor);
ylim([0 (mNF+(mNF*0.5))]);
title('Cardiac Noise Criteria')
xlabel('Channel #');
ylabel('Noise Floor Index');
xlim([1 m])
% Channel Distance plot
subplot(5,2,4);
hold on
for a=1:m
    if dists(a) >= sdDistRange(1) & dists(a) <= sdDistRange(2)
        plot((a),dists(a),['xk']);
    else
        plot((a),dists(a),['xr']);
    end
end
plot([1 m],sdDistRange(1)*[1 1],'--k');
plot([1 m],sdDistRange(2)*[1 1],'--k');
ylim([0 dists(2)+10]);
title('Source-Detector Distances')
xlabel('Channel #');
ylabel('Distance (mm)');
xlim([1 m])
%Corresponding channel plot
subplot(2,2,4);
drawProbeGood_CD(X.SD,chansTFall);
title('Good & Bad Channel Display')
end
return



function thelist = getNirsList(SD,ml,freq)
% will tell you which channels are in that frequency
theInd = find(SD.Lambda == freq);
thelist = find(ml(:,4) == theInd);
return
