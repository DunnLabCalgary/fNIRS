function [chans,SNR,dists] = hmrPruneChannels_CD(d,SD,DISTminmax,correspondTF,ODrange,SNRmax)
%chans = nirsPruneChannels(d,SD,[DISTminmax,correspondTF,ODrange,SNRmax])
%
% This function will do everything it can to determine if a channel is
% valid or not. The following criteria are tested:
%  - Power range: 80-125 dB (by default, unless changed)
%  - stddev/mean=SNR < 5 (by default, unless changed)
%  - optional: oxy and deoxy corresponding channels are "good"
%  - optional: source-detector separation distance exclusion
%
%INPUTS:
% d = data by channel raw O.D.
% SD = data structure describing the probe
% [DISTminmax ]= the min/maximum distance permitted for source-detector pairs
%   if not specified, no distance criteria will be used
%   if only one value is specified, that is taken as max
%   if a vector of size 2 is specified, it is of the form [min max]
% correspondTF = 1 if yes, 0 if no. Default: yes.
% ODrange = lower and upper dB range for power (optional, default 80-125)
% SNRmax = max SNR value (optional, default 5)
%
%OUTPUTS:
% chans = indices of channels that are ok
%
%CALLS:
% (none)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preset values
if ~exist('correspondTF'); correspondTF=1;  end;
if ~exist('ODrange'),   ODrange=[];         end;
if ~exist('SNRmax'),    SNRmax=[];          end;
if isempty(ODrange),    ODrange = [80 125]; end;
if isempty(SNRmax),     SNRmax=5;           end;

%% check for o.d. range and SNR max
ddb = 20*log10(d);
SNR = std(ddb)./mean(ddb);
mins = min(ddb);
maxs = max(ddb);
chanList = (SNR<SNRmax) & (mins>ODrange(1)) & (maxs<ODrange(2));

%% check for valid distances
ml = SD.MeasList;
spos = SD.SrcPos;
dpos = SD.DetPos;
m2 = size(ml,1)/2; % number of channel pairs
sSet = spos(ml(:,1),1:2);
dSet = dpos(ml(:,2),1:2);
dists = sqrt((sSet(:,1)-dSet(:,1)).^2 + (sSet(:,2)-dSet(:,2)).^2);
if strcmp(SD.SpatialUnit,'cm') % coded in cm not mm so multiply by 10
    dists = dists*10;
end
if ~exist('DISTminmax'); DISTminmax=[]; end;
if length(DISTminmax)==1
    chanList = chanList & (dists<DISTminmax)';
elseif length(DISTminmax)==2
    chanList = chanList & (dists>DISTminmax(1))' & (dists<DISTminmax(2))';
end

%% Check that oxy and deoxy corresponding channels are good
if correspondTF==1
    % generate colorlist 2x(# of corresponding chans) - row1 is 690 row is 830
    colorList = [];
    wavelengths = SD.Lambda'; %Chris add to grab wavelengths from ml. Compatible with NIRx
    for w=1:2
        theInd = find(SD.Lambda == wavelengths(w));
        colorList(w,:) = find(ml(:,4) == theInd);
    end
    temp = [chanList(colorList(1,:)) & chanList(colorList(2,:))];
    chanList = [temp temp];
end

%% Get indices
chans = find(chanList);
return
