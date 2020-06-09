function PostAnalysis_CD(Data)
    
    regionLabels = ["R1Intra", "R1R2", "R2Intra"];
    roiCombinations = {};
    wcoh = {};
    Cohfft = {};
    phase = {};
    stims = {};
    stimNames = {};
    for rc = 1:length(regionLabels)
        
        status = 0;
        flag = 0;
        
        for S_count = 1:length(Data.SOI)
            
            sub = [Data.SOI{S_count}]; % sub string

             if exist('w', 'var')
                 waitbar(status, w, ['Analyzing ' strrep(sub, '_', ' ') ' for ' char(regionLabels(rc))]);
             else
                 w = waitbar(status, ['Analyzing ' strrep(sub, '_', ' ') ' for ' char(regionLabels(rc))]);
             end
            
            status = (S_count / length(Data.SOI));
            
            %load processed file
            
            procFile = ([Data.Path.study_path '/Processed' filesep sub]);
            load([procFile filesep 'Processed.mat']);
            
            %get wavelet file
            WFile = [procFile filesep 'PostProcessed' filesep 'Wavelet' char(regionLabels(rc)) '.mat'];
            
            
            if exist(WFile)
                flag = 1;
                load(WFile);
                display(['Analyzing Subject ' sub ' ' char(regionLabels(rc))]);
                
                %construct stim times
                load([Data.Path.RawData filesep sub '.nirs'], '-mat'); %important to load this first so t doesint over ride the Proc.t
                
                t = Proc.t;
                smod = Proc.smod;
                
                Rects = Snew.markerPos.StimRects;
            
                si = {};
                for idx = 1:length(Rects)
                    r = Rects(idx);
                    lb = r{1}(1);
                    ub = lb + r{1}(3);
                    stimLogic = t>lb & t<ub;
                    si{idx} = stimLogic; 
                end
                
                %get good ch
                nCh = 1:size(Proc.dc,3);
                GoodCh = Proc.QC.dqChans(Proc.QC.dqChans<=max(nCh));
                
                FZband = Wavelet.FZband;
                Hemo = Wavelet.HB;
                t = Wavelet.t;
                
                if size(t,1) == 1 %ensure t is a column array
                    t = t';
                end
                
                fs = 1/(t(2)-t(1)); % sampling rate
                
                RegionalData = Wavelet.RegionalData;
                regionalNames = RegionalData.regionNames;
                
                if isfield(RegionalData, 'R1Channels')
                    RegionalData.R1Channels = intersect(GoodCh, RegionalData.R1Channels);
                end

                if isfield(RegionalData, 'R2Channels')
                    RegionalData.R2Channels = intersect(GoodCh, RegionalData.R2Channels);
                end
                
                R1Channels = [];
                R2Channels = [];

                if rc == 1
                    if S_count == 1
                        roiCombinations{end + 1} = [1,1];
                    end
                    R1Channels = RegionalData.R1Channels;
                    R2Channels = R1Channels;
                    
                elseif rc == 2
                    if S_count == 1
                        roiCombinations{end + 1} = [1,2];
                    end
                    R1Channels = RegionalData.R1Channels;
                    R2Channels = RegionalData.R2Channels;
                elseif rc == 3
                    if S_count == 1
                        roiCombinations{end + 1} = [2,2];
                    end
                    R1Channels = RegionalData.R2Channels;
                    R2Channels = R1Channels;
                end
                
                filePath = [procFile filesep char(regionLabels(rc)) '_'];
                regionCombLabel = char(regionalNames(rc));
                X = analysisPlot(si, Proc.dc, Wavelet, fs, R1Channels, R2Channels, Hemo, FZband, t, t(find(smod)), filePath, [char(regionCombLabel) ' ' strrep(sub, '_', ' ')], {RegionalData.R1Name, RegionalData.R2Name}, Data.postUI.analysisFigures, Snew.StimNames);
                
                fftcoh{S_count} = X.fftcoh;
                tvcoh{S_count} = X.tvcoh;
                
                tvwt = X.tvwt;
                
                if rc == 1
                    wt{S_count, 1} = squeeze(tvwt);
                elseif rc == 2
                    wt{S_count, 1} = squeeze(tvwt(1,:,:));
                    wt{S_count, 2} = squeeze(tvwt(2,:,:));
                else 
                    wt{S_count, 2} = squeeze(tvwt);
                end
                
                if ~isequal(R1Channels, R2Channels)
                    tvphase{S_count} = X.tvphase;
                end
                
                if ~exist('si')
                    stims{S_count} = {};
                    stimNames{S_count} = {};
                else
                    stims{S_count} = si;
                    stimNames{S_count} = Snew.StimNames;
                end 
                
            end
            
        end
        %end S_count
        
        if flag
            Cohfft{end + 1} = fftcoh;
            wcoh{end + 1} = tvcoh;
            
            if rc == 2
                phase{end + 1} = tvphase;
            else
                phase{end + 1} = {};
            end 
        end
        
        clear fftcoh tvcoh tvwt tvphase
        
    end
    
    if exist('w', 'var')
        waitbar(0.99, w, ['Outputting Data Tables']);
    end
    
    output_DataTables(Data, FZband(1), FZband(2), Hemo, fs, wcoh{1}, wt, phase, Cohfft, stims, sub, roiCombinations, stimNames) %for each region comparison we output one data table

    if exist('w', 'var')
        close(w)
    end


 
        