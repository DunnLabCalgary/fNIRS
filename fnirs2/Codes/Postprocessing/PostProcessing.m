%Postprocessing: 
%Takes in the preprocessed Data set and user inputs from the UI and
%calculates the continuous wavelet transform parameters (phase coherence,
%and power)

function PostProcessing(Data, PostUI)

status = 0; 

for S_count = 1:length(Data.SOI) % loop counter for all subjects selected
    sub = [Data.SOI{S_count}]; % subject variable string
    procFile = [Data.Path.study_path '/Processed/' sub];
    load([procFile filesep 'Processed.mat']); % load subject's preprocessed data struct
    cd(procFile);
    
    
    %load stims
    load([Data.Path.RawData filesep sub '.nirs'], '-mat');
    Rects = Snew.markerPos.StimRects;
    
    
    display(['Postprocessing Subject ' sub]); %console message
    
    % Calculate the pairwise continuous wavelet transform parameters (phase coherence,
    % and power).
    % input: Processed data & UI inputs
    % output: Wavelet struct with all wavelet parameters for the scan
    regionComparison = [Data.RegionalData.R1Intra, Data.RegionalData.R1R2, Data.RegionalData.R2Intra];
    regionNames = [convertCharsToStrings([Data.RegionalData.R1Name ' ' 'Intra']) , convertCharsToStrings([Data.RegionalData.R1Name ' ' Data.RegionalData.R2Name ' ' 'Inter']), convertCharsToStrings([Data.RegionalData.R2Name ' ' 'Intra'])];
    regionLabels = ["R1Intra", "R1R2", "R2Intra"];
    for rc = 1:length(regionComparison)
        if regionComparison(rc)
            
            Data.RegionalData.rc = rc;
            [Wavelet] = ContinuousWaveletAnalysis(PostUI.FZband,Proc, Data.RegionalData, PostUI.Hemo, Rects, [procFile filesep 'Processed.mat'], sub);
            
            Wavelet.HB = PostUI.Hemo;
            Data.RegionalData.regionNames = regionNames;
            Wavelet.RegionalData = Data.RegionalData;
            Wavelet.FZband = PostUI.FZband;
            
            F1 = strrep(num2str(min(PostUI.FZband)),'.','p'); % min fz band for string (replace periods with p)
            F2 = strrep(num2str(max(PostUI.FZband)),'.','p'); % max fz band for string
            
            if ~exist([procFile filesep 'PostProcessed'], 'dir')
                mkdir([procFile filesep 'PostProcessed'])
            end
            
            cd([procFile filesep 'PostProcessed'])
            
            save(['Wavelet' char(regionLabels(rc))], 'Wavelet') % save Wavelet data in processed folder
            
            if PostUI.Fig
                plotStaticCoherenceMatrix(Wavelet.wcoh, sub, Wavelet.ChX, Wavelet.ChY, PostUI.Hemo, char(regionNames(rc))) % plot & save pairwise coherence matrix
                WFile = ['WCohMatrix_' F1 '_' F2 '_Hz_' char(regionNames(rc)) '.jpg'];
                saveas(gcf,WFile)
                close
            end
            
            
        end
    end
    
end

disp('Done postprocessing')
end

function plotStaticCoherenceMatrix(X, sub, ChannelsX, ChannelsY, HB, regionLabel)

    h2 = figure('Name',['Static Matrices ' sub ],...
                    'Color','w',...
                    'NumberTitle','off');
    
    ha2 = tight_subplot(2,2,...
                    [.08 .05],[.05 .04],[.04 .01]);
                  
        for idx = 1:numel(HB)
            axes(ha2(idx))
         
            %before: X.HB{idx} = (numberChannels * numberChannels * numberDatapoints)
            %after applying nanmean: (numberChannels * numberChannels)
            %where the value in each place is the mean coh value,
            %numberChannels is ALL channels not just good ones
            im = nanmean(X.(HB{idx}),3);
            %take only the good channels to plot
            im = im(ChannelsX, ChannelsY);
            imagesc(im,[0 1]);
            yticks(1:length(ChannelsX));yticklabels(ChannelsX);%fit this wierdness CHX and CHY should be switched here 
            xticks(1:length(ChannelsY));xticklabels(ChannelsY);
            ax = gca; ax.FontSize = 6;
            disp(regionLabel);
            title([sub ' ' HB{idx} ' Time-averaged Matrix ' regionLabel], 'FontSize', 10,'FontWeight','bold');  
        end 
        
end
