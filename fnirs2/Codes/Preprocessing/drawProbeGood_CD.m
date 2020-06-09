function drawProbeGood_CD(SD,goodInds)
%
% drawProbeGood(SD,goodInds)
% given an SD structure and a set of "good" indices from measlist
% a probe will be drawn with good or bad optodes.
%

ml = SD.MeasList;

if length(goodInds)~= size(ml,1)
    disp('Size mismatch error.')
    return
end

goodSrcs = zeros(1,SD.nSrcs);
for srcCount=1:SD.nSrcs
    sInds = find(ml(:,1) == srcCount);
    goodSrcs(srcCount) = sum(goodInds(sInds))>0;
end

% Now check the same deal per detector
goodDets = zeros(1,SD.nDets);
for detCount=1:SD.nDets
    dInds = find(ml(:,2) == detCount);
    goodDets(detCount) = sum(goodInds(dInds))>0;
end

for detCount=1:SD.nDets
    x=SD.DetPos(detCount,1);
    y=SD.DetPos(detCount,2);
    tf = goodDets(detCount);
    textColor = 0.7 * [(1-tf) tf 0];
    text(x,y,['o ' num2str(detCount)],'fontsize',12,'fontWeight','bold','Color',textColor);
end
for srcCount=1:SD.nSrcs
    x=SD.SrcPos(srcCount,1);
    y=SD.SrcPos(srcCount,2);
    tf = goodSrcs(srcCount);
    textColor = 0.7 * [(1-tf) tf 0];
    text(x,y,['x ' num2str(srcCount)],'fontsize',12,'fontWeight','bold','Color',textColor);
end
% now put in the connections
hold on
for PP=1:2
    for chan=1:size(ml,1)
        chanSrc = ml(chan,1);
        chanDet = ml(chan,2);
        chanType = ml(chan,4);
        chanX = mean([SD.SrcPos(chanSrc,1) SD.DetPos(chanDet,1)]);
        chanY = min(SD.SrcPos(chanSrc,2),SD.DetPos(chanDet,2));
        Yrange = abs(SD.SrcPos(chanSrc,2)-SD.DetPos(chanDet,2));
        chanY = chanY + ((chanType-1)*0.2 + 0.4)*0.8;
        tf = goodInds(chan);
        if PP==1
            plot([SD.SrcPos(chanSrc,1) SD.DetPos(chanDet,1)],...
                [SD.SrcPos(chanSrc,2) SD.DetPos(chanDet,2)]+0.1*(chanType-1),...
                'lineWidth',2,'color',[.5+(chanType*.2) .7 .5+((2-chanType)*.2)]);
        else
            if chan <= size(ml,1)/2 % Chris added so to not stack second wavelength number
            textColor = 0.7*[(1-tf) tf 0];
            text(chanX,chanY,num2str(chan),'Color',textColor);
            plot([SD.SrcPos(chanSrc,1) SD.DetPos(chanDet,1)],...
                [SD.SrcPos(chanSrc,2) SD.DetPos(chanDet,2)]+0.1*(chanType-1),...
                'lineWidth',1,'color',textColor);
            end
        end
    end
end

x=[SD.DetPos(:,1); SD.SrcPos(:,1)];
y=[SD.DetPos(:,2); SD.SrcPos(:,2)];
xl = [min(x) max(x)];
yl = [min(y) max(y)];
dx = diff(xl)*0.1;
dy = diff(yl)*0.1;

xlim([xl(1)-dx xl(2)+dy]);
ylim([yl(1)-dy yl(2)+dy]);
return