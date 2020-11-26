function helperPlotCoherence(wcoh,t,F,coi,xlab,ylab, ChX, ChY, sub, HB)
%   This helper function is provided solely in support of the featured
%   example wcoherdemo.m
%   This function may be changed or removed in a future release.

Yticks = 2.^(round(log2(min(F))):round(log2(max(F))));
figure('Name',['WCoherence' num2str(ChX) num2str(ChY) sub HB]);
imagesc(t,log2(F),wcoh);
set(gca,'YLim',log2([min(F),max(F)]), ...
    'layer','top', ...
    'YTick',log2(Yticks(:)), ...
    'YTickLabel',num2str(sprintf('%.2f\n',Yticks)), ...
    'layer','top','YDir','normal');
hold on;
plot(t,log2(coi),'w--');
xlabel(xlab);
ylabel(ylab);
title([ChX ChY sub HB]);
hcol = colorbar;
hcol.Label.String = 'Magnitude-Squared Coherence';
title('Wavelet Coherence');

if ~exist([pwd filesep 'WaveletCoherence'])
    mkdir([pwd filesep 'WaveletCoherence'])
end

f = [pwd filesep 'WaveletCoherence' filesep 'WCoherence' num2str(ChX) num2str(ChY) sub HB '.png'];
saveas(gcf,f)
close
