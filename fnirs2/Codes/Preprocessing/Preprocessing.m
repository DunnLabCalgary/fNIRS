function [Data] = Preprocessing(Data, PreUI)

% This function is part of an automated pipeline for processing
% fNIRS time-series. Input from the user interface is used for assessing signal
% quality and rejecting poor data. Figures allow user to visualize data quality 
% and ensure that the algorithms and inputs are adequately cleaning the data.
%
% Motion artifact can be removed manually (MARA), or automatically (Spline
% nterpolation). The user can also choose artifact removal + smoothing
% (Sovitzky-Golay filter) which is ideal for analyzing the hrf, but may not
% be ideal for coherence analysis.
%
% After cleaning and filtering, the processed time-series is plotted alongside
% unprocessed data for visual confirmation that accurate preprocessing
% was performed.
%
% Input data should be raw data in .nirs format.
%
% (1) Initialize 
%       - Load raw data and setup
%       
% (2) Check Data Quality
%       - Assess signal quality and discard channels that don't meet
%       criteria. Methods and code are adapted from Homer2.8 toolbox.
%       - Key code: hmrDataQualityCheck_CD()
%
% (3) Convert Raw Data to Optical Density
%       - Code from Homer2.8 toolbox.
%       - Key code: hmrIntensity2OD()
%
% (4) Motion Artifact Detection/Correction
%       - Spline Interpolation, MARA, or SG filter
%       - Methods and code adapted from Homer2.8 and from Scholkmann et al. 2010. See citation
%       below.
%       - Key code: MotionCorrection_MARA_V1()
%                   hmrMotionCorrectSplineSG_CD()
%
% (5) Convert Optical Density (DoD) to ?Hemoglobin Concentration (dc)
%       - The differential pathlength factor (DPF) is a key factor in the Modified Beer Lambert
%           equation which is used to calculate hemoglobin concentration.
%           DPF has been shown to vary with age so CalcDPF adjusts the DPF
%           by age (input).If age is available in demographic table, it is
%           pulled to use in CalcDPF, before conversion to concentration.
%           See below for citation. hmrOD2Conc() is from Homer2.8
%       - Key code: CalcDPF()
%                   hmrOD2Conc()
%
% (6) Detrend, z-score normalize, and Downsample (5Hz)
%
% INPUT
%   Data:       Struct that contains all user input info
%
% OUTPUT
%   Data:       Updated Data struct
%
% SAVED
%   Each subjects preprocessed data is saved to Preprocessed folder in the study folder 
%
% -------------------------------------------------------------------- %
% Citations
%
% 1) Motion correction with MARA
%       "How to detect and reduce movement artifacts in near-infrared
%       imaging using moving standard deviation and spline interpolation."
%       Physiological Measurement, 31, 649-662. Scholkmann et al. 2010.
%
%       Requires toolbox 'MATLAB_ToolsForUCL_Felix_2014'.
%
% Sovitzky-Golay and Prune info.....
%
% Calculate age-dependent DPF (CalcDPF).
%       "General equation for the differential pathlength factor of the
%       frontal human head depending on wavelength and age".
%       Scholkmann and Wolf 2013.
%
% Homer2 toolbox required for conversion algorithms.
%       hmrIntensity2OD
%       hmrOD2Conc
%
%
%
%%% Author: Chris Duszynski
%%% Dunn Lab - Experimental Imaging Centre
%%% January 2019
%%% Contributers: Lia Hocke, Felix Scholkmann, Jeff Dunn, Carter Randall
%
% --------------------------------------------------------------------- %

status = 0; %used for loading bar in GUI
% Data.Path.stim = 'C:\Users\Chris\Google Drive\Matlab\Test2019\TimingNirx'; % need to make User Input

for S_count = 1:numel(Data.SOI) % SOI = subjects of interest (subjects selected in gui)
    
    %% Initialize

    sub = (Data.SOI{S_count}); % Subject of interest string
    
    if exist('w','var')
        waitbar(status,w,['Initializing Subject ' strrep(sub, '_', ' ')]); %strrep is to remove underscores and replace with spaces cause other wise the letter is subscripted (for niceness purposes only)
    else
        w = waitbar(0 ,['Initializing Subject ' strrep(sub, '_', ' ')]);
    end
    
    filename = [Data.Path.RawData filesep sub '.nirs'];  
    % future... check for nirx or nirs format and then convert automatically to nirs if nirx files
    
    X = load(filename,'-mat'); % load file 
    X.s = X.s(:,1);%to fix wierd data from sabrain where s had 2 columns 
    fs = 1/(X.t(2)-X.t(1)); % sampling rate
    
    if PreUI.NIRx == 1
        X.d = X.d*10^7; % Multiply to be on same scale as TechEn so taking the log does not make negative values. Need to check and confirm this has no negative effects
    end
    
    %% Quality Check And update MeasList Act
    
    QualityCheck = performQualityCheck(Data.Path.study_path, PreUI, X, sub);
    %ismember: logical where the data in MeasListAct is contained in the dqChans so
    %measListAct will become a logical array with 1 in places with good
    %channels and 0 else where.
    X.SD.MeasListAct = double(ismember(1:numel(X.SD.MeasListAct),QualityCheck.dqChans)'); %update SDMeasListAct with 0's for bad channels
    
    %% Convert raw data to optical density (DoD)
    
    d = inf_resolution_d(X.d); % take care of inf's or nans in the data
    dod = hmrIntensity2OD(d); % for motion correction
    dodnomc = dod; %was: hmrIntensity2OD(d); % for non-corrected visual 
    dodmc = inf_resolution(d,dod); % take care of inf's or nans in the data
    
    %% Assess Signal Quality with Figures
    
    if ~exist('w','var')
        w = waitbar(0 ,['Initializing Subject ' sub]);
    else
        waitbar((.50*S_count)/length(Data.SOI),w,['Assessing signal quality - ' sub]);
    end
    
    %figs
    Ch_w1 = find(X.SD.MeasList(:,4)==1); % channel list for wavelength 1
    Ch_w2 = find(X.SD.MeasList(:,4)==2); % channel list for wavelength 2
    
    if PreUI.Figs.sigQuality == 1
        % Define axes stretched past the default positions
        m = 10; % 9 rows
        n = ceil(length(Ch_w1)/m); % # of columns determined by # of channels to plot
        l = m*n; % total # of plots (SNR counts as entire row of plots)
        Rect = [0.05, 0.05, 0.93, 0.9];
        AxisPos = moPlotPos(n,m, Rect);
        % Figure handle
        h1=figure('name',[sub ' Signal Quality'],...
            'NumberTitle','off',...
            'Color','White',...
            'units','normalized',...
            'outerposition',[0 0 1 1],...
            'DefaultaxesPosition',Rect);
        
        % plot OD raw 
        for ii = 1:size(Ch_w1)
            subplot(m,n,ii)
            set(gca,'Position',AxisPos(ii,:))
            plot(X.t,dodnomc(:,Ch_w1(ii)),'r');hold on;
            plot(X.t,dodnomc(:,Ch_w2(ii)),'b')
            xlim([0 X.t(end)]);
            ylabel('OD');
            if ~ismember(Ch_w1(ii),QualityCheck.dqChans) 
                yl=ylim;plot(X.t,linspace(yl(1),yl(2),length(X.t)),'--r','LineWidth',2)
            end
            grid on
            title(num2str(Ch_w1(ii)));
        end
        
        figure(h1);
        saveas(gca,[Data.Path.study_path filesep 'Processed' filesep sub '_SignalQuality'],'jpg');
        close
    end
    clearvars -except Data PreUI d dod dodmc dodnomc filenm QualityCheck S_count sub w fs status X
            
        if ~exist('w','var')
            w = waitbar(0 ,['Performing Motion Correction - ' sub]);
        else
            waitbar((.50*S_count)/length(Data.SOI),w,['Performing Motion Correction - ' sub]);
        end
           
    %% Motion correction
       disp('performing automated (SplineSG) motion correction')
        %% Automated motion correction
    % Altered hmrMotionCorrectSplineSG to accept input of whether or not to
    % smooth with Sovitzky-Golay (not recommended for time-frequency
    % analysis). Motion is detected with hmrtInc_baselineshift_Ch which
    % detects motion based on std variations, gradient outliers, baseline
    % shift and spikes, in consideration of heart rate amplitude, and then
    % corrected using Spline Interpolation. 
        p = 0.99;
        dodmc = hmrMotionCorrectSplineSG_CD(dodmc, d, X.t, X.SD, p, PreUI.SDseconds, PreUI.SGfilter);
    
    
    if exist('h1','var')
    close (h1)
    end  
    
    %% Convert DoD to dc [Hb Conc]
    % Extract age from demographic table for differential pathlength factor
    % in the modified beer lambert equation (conversion to concentration)
    
    if ~exist('w','var')
        w = waitbar(0 ,['Initializing Subject ' sub]);
    else
        waitbar((.75*S_count)/length(Data.SOI),w,['Processing Time-Series - ' sub]);
    end
    
%     if size(Data.Demographics,1) == length(Data.PreUI.SOI)
%         Age = Data.Demographics.Age(S_count);
%     else %need to fix bug here to deal with demographic list doesnt exactly match Subject list
        dm = table2array(Data.Demographics(:,1));
%         if isnumeric(dm{S_count,1})
%             ct = find(table2array(dm(:,1)) == Data.PreUI.SOI(S_count));
%         else
%             ct = find(strcmp(table2array(dm),Data.PreUI.SOI{S_count}));
%         end

        % find age that corresponds to subject name
        for idx = 1:numel(dm)
            q = dm{idx};
            if strcmp(q,Data.SOI{S_count})
                ct = idx;
            end
        end
        
        if exist('ct', 'var')
            Age = Data.Demographics.Age(ct);
        else
            Age = 30;
        end
        
        
%     end
    
    if isempty(Age) || isnan(Age)
        Age = 30;    % If no age available, use default age of 30
        warning('DPF was not age corrected. Default age of 30 was used')
    end
    
    for lm = 1:length(X.SD.Lambda)
        DPF(lm) = CalcDPF (Age,X.SD.Lambda(lm));
    end
    
    clearvars Age ct
    
    % Convert Optical Density (DoD) to dConcentration (dc)
    dc = hmrOD2Conc(dodmc,X.SD,DPF); %corrected
    dcnomc = hmrOD2Conc(dodnomc,X.SD,DPF); %uncorrected
    
    %% Detrend, z-score normalize, and Downsample (5Hz)
    
    TFs = PreUI.TFs;    % sampling rate target after downsampling (Hz)
    Ch = 1:length(X.SD.MeasList)/2; %channels of interest
    GC = QualityCheck.dqChans; % need to implement with chanList
    
    % organize OD data based on channels of interest for plot
    dd(:,:,1) = dodnomc(:,Ch);
    dd(:,:,2) = dodnomc(:,Ch+length(dodnomc(1,:))/2);
    
    if isempty(GC) 
        warndlg('All channels where discarded for subject:')
        clearvars -except Data PreUI S_count status w
        continue
    end
    
    % Detrend, decimate (downsample), and z score normalize
    t1 = X.t;
    for ii = 1:length(Ch)
        if ismember(Ch(ii),GC)
            for HB = 1:3
                
                if PreUI.Detrend %detrend
                    xd(:,HB,ii) = detrend(dc(:,HB,ii));
                else
                    xd(:,HB,ii) = dc(:, HB, ii);
                end
                
                if PreUI.Downsample && fs > TFs %% downsample
                    n = round(fs/TFs);
                    dcf(:,HB,ii) = decimate(xd(:,HB,ii),n); % Fs -> TF 
                    t = 0+1/TFs:1/TFs:length(dcf)/TFs;t=t';
                    times = t1(find(X.s));
                    smod = zeros(length(t), 1);
                    
                    for k = 1:length(times)
                        [~,tIndex] = min(abs(t-times(k)));
                        smod(tIndex) = 1;
                    end
                    
                else 
                    dcf(:,HB,ii) = xd(:,HB,ii);
                    t = X.t;
                    
                    smod = X.s;
                end
            end
            
            % zscore normalization
            if PreUI.Normalize
                dcfz = zscore(dcf,0,1); 
            end

            % Figures for assessing preprocessing, plot each channel
            if PreUI.Figs.corrSig == 1
                
                if ~exist([Data.Path.RawData filesep 'Processed' filesep sub])
                    mkdir([Data.Path.RawData filesep 'Processed' filesep sub]);
                end
                
                % Figure handle
                
                h=figure('name',[sub ' - Ch ' num2str(Ch(ii))],...
                    'Color','White',...
                    'units','normalized',...
                    'outerposition',[0 0 1 1]);
                % Make/adjust subplot axes handle
                ha1 = tight_subplot(2,2,...
                    [.08 .05],[.1 .04],[.04 .01]);
                % Plot 1: Raw DOD
                axes(ha1(1));
                y1=dd(:,ii,1);y2=dd(:,ii,2); % DOD 690 and 830
                plot(t1,y1,'k','LineWidth',1); hold on;
                plot(t1,y2,'Color',[0 0 0]+0.5,'LineWidth',1);
                title (['Raw Data - Ch ' num2str(Ch(ii))])
                ylabel('Delta Optical Density');xlabel('time (s)')
                legend([num2str(X.SD.Lambda(1)) ' nm'],...
                    [num2str(X.SD.Lambda(2)) ' nm'],'Location','southwest')
                % Plot 2: Raw Concentration
                axes(ha1(2));
                plot(t1,dcnomc(:,2,ii),'b','LineWidth',1);
                hold on
                plot(t1,dcnomc(:,1,ii),'r','LineWidth',1);
                title(['Uncorrected & Unfiltered Concentration'])
                ylabel(['Delta Concentration (M)']);xlabel('time (s)')
                legend(['Hb'],...
                    ['HbO'],...
                    'Location','southwest');
                % Plot 3: Cleaned/corrected Concentration
                axes(ha1(3));
                plot(t1,dc(:,2,ii),'b','LineWidth',1);
                hold on
                plot(t1,dc(:,1,ii),'r','LineWidth',1);
                title(['Motion Corrected Concentration'])
                ylabel(['Delta Concentration (M)']);xlabel('time (s)')
                legend(['Hb'],...
                    ['HbO'],...
                    'Location','southwest');
                % Plot 4: Motion Corrected, Smoothed, and Detrended
                % Concentration
                axes(ha1(4));
                plot(t,dcfz(:,2,ii),'b','LineWidth',1);hold on
                plot(t,dcfz(:,1,ii),'r','LineWidth',1);
                title('Detrended, normalized, & downsampled')
                ylabel('Delta Concentration (M)');xlabel('time (s)')
                legend(['Hb'],['Hb0'],'Location','southwest')
                if ~exist([Data.Path.study_path filesep 'Processed' filesep sub])
                    mkdir([Data.Path.study_path filesep 'Processed' filesep sub]);
                end
                saveas(gca,[Data.Path.study_path filesep 'Processed' filesep sub filesep ...
                    'Ch' num2str(Ch(ii))],'jpg');close
            end
        end
    end
    
    if ~exist('w','var')
        w = waitbar(0 ,['Initializing Subject ' sub]);
    else
        status = S_count/length(Data.SOI);
        waitbar(status,w, ['Saving data - ' sub]);
    end
    
    %% Save preprocessed data struct
    if PreUI.Normalize && exist('dcfz', 'var')
        Proc.dc = dcfz;
    elseif PreUI.Downsample && exist('dcf', 'var')
        Proc.dc = dcf;
    elseif PreUI.Detrend && exist('xd', 'var')
        Proc.dc = xd;
    else
        Proc.dc = dc;
    end
    Proc.smod = smod; %%s vector to match length of t if downsampled
    Proc.t = t;
    Proc.SD = X.SD;
    Proc.aux = X.aux;
    Proc.s = X.s;
    Proc.QC = QualityCheck;
    Proc.DateTime = datetime('today');  % Date and time stamp of processing
    
%     % Preprocessed data goes in subject processed folder
    if ~exist([Data.Path.study_path filesep 'Processed' filesep sub])
        mkdir([Data.Path.study_path filesep 'Processed' filesep sub]);
    end
    save([Data.Path.study_path filesep 'Processed' filesep sub filesep 'Processed.mat'],'Proc');
    
    % Save XL sheet with details **** need to make this output with all
    % Proc info and write to csv
     %xlswrite([Data.Path.study_path '\Preprocessed\' sub '\DiscardedChannels.xlsx'],Data.Subjects.(sub).Info.Discards)
    
    clearvars -except Data PreUI S_count status w
    
end
close(w)
disp('done preprocessing');
end

function [QualityCheck] = performQualityCheck(study_path, PreUI, X, sub)
    %% Check Data Quality
    
    %%%get values from UI struct%%%
    motionTF = 0;
    cardiacPeak = 0;
    modeFlags = [PreUI.dbRangeBox, PreUI.SNRbox, PreUI.sdRangebox, PreUI.hbMatchbox, cardiacPeak, motionTF, PreUI.Figs.cardiacFreq, PreUI.Figs.DataQuality];                      
    dBrange = [PreUI.dBmin, PreUI.dBmax]; 
    sdDistRange = [PreUI.SDmin,PreUI.SDmax]; % mm coded in SD, not cm
    snrMax = PreUI.SNRmax; 
    
    %%%channel list is the indeces of remaining channels (good channels) corresponding to d in the MeasListAct%%%
    [channelList, detailedInfo] = hmrDataQualityCheck_CD(X,modeFlags,dBrange,sdDistRange,snrMax); % remove channels that don't meet criteria
    
    if ~exist([study_path filesep 'Processed'])
        mkdir([study_path filesep 'Processed']);
    end
  
    %% clean up the folder saving structure here a bit later
    %%%save figures%%%
    if PreUI.Figs.DataQuality == 1
        if ~exist([study_path filesep 'QualityCheck'])
            mkdir([study_path filesep 'QualityCheck']);
        end
        
    saveas(gcf,[study_path filesep 'Processed' filesep sub 'QualityCheck'],'jpg');close;
    end
    
    if PreUI.Figs.cardiacFreq == 1
        if ~exist([study_path filesep 'QualityCheck'])
            mkdir([study_path filesep 'QualityCheck']);
        end
    saveas(gcf,[study_path filesep 'Processed' filesep sub 'CardiacPlot'],'jpg');close;
    end 
    
    %%%Collect output values from quality check%%%
    QualityCheck.modeFlags = modeFlags; clear modeFlags
    QualityCheck.dBrange = dBrange; clear dBrange
    QualityCheck.SDrange = sdDistRange; clear sdDistRange
    QualityCheck.SNRmax = snrMax; clear snrMax
    QualityCheck.dqChans = channelList; clear channelList
    QualityCheck.detailedInfo = detailedInfo; clear detailedInfo

end

