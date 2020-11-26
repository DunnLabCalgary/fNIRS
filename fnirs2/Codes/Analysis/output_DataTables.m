

function [] = output_DataTables(Data,F1,F2,hemo,fs,wcoh,wt,phase,Cohfft, si, sub, roiCombinations, tk)

%% make name strings: pass in name data from gui, chechk intra and inter specifications.

rNames = ["Region1", "Region2"];
for r = 1:length(roiCombinations)
    ROIs{r} = strjoin([rNames(roiCombinations{r}(1)) '->' rNames(roiCombinations{r}(2))]);
end 


%% Coherence
% for each subject for each combination of hb, task, roi take wcoh mean and
% write to csv
id = 1;

for roi = 1:length(ROIs)
    for sb = 1:length(Data.SOI)
        for hb = 1:numel(hemo)
            for stim = 1:length(si{sb})
                
    %               sn = Data.AllSubs(sb);
                    taskVector = si{sb}{stim}'; 

            %coh stuff
                    w = wcoh{roi, sb}(hb, taskVector); %wcoh vector for region, sub, hb 
                    x = nanmean(w); %nanmean coh during stim
                    y = nanstd(w); % std coh value for roi x hb x task
                    
                    %phase stuff
                    if ~isequal(phase{roi}, {}) 
                        p = phase{roi}{sb}(hb,taskVector); %phase from region, sub, hb
                        xph = nanmean(p); % mean coh value for roi x hb x task
                        xabs = nanmean(abs(p)); % abs phase value for roi x hb x task
                        yph = nanstd(p); % std coh value for roi x hb x task
                    else
                        xph = 0;
                        xabs = 0;
                        yph = 0; % std coh value for roi x hb x task
                    end

                    yfft = nanmean(Cohfft{roi}{sb}(hb,:)); % fft coh value for roi x hb

                    %save for later table output
                    m(id,:) = [sub,ROIs{roi},tk{stim, sb},hemo{hb},yfft,x,y,xph,xabs,yph]; 

                    id = id + 1;
                
               
            end
        end
    end
end

header_string = {'ID','Connection','Task','Hb','Cohfft','CohMean','CohSTD','PhaseMean','AbsPhaseMean','PhaseSTD'};
m = array2table(m,'VariableNames',header_string);
%% make sure dir exists here (output)
filename = [Data.Path.study_path filesep 'Processed' filesep sub filesep 'Analysis' filesep 'WCOHresults_' num2str(F1) '_' num2str(F2) 'Hz.csv'];
writetable(m,filename); % write to csv
clear m header_string


%% Wavelet Power
% for each subject for each combination of hb, task, roi take wcoh mean and
% write to csv
id = 1;
ROISplit{1} = ['Region1']; %%make these the region names 
ROISplit{2} = ['Region2'];

for region = 1:length(wt)
    for sb = 1:length(Data.SOI)
        for hb = 1:numel(hemo)
            for stim = 1:length(si{sb})
            
                taskVector = si{sb}{stim};
                w = wt{sb, region}(hb, taskVector);
                x = nanmean(w);
                y = nanstd(w);
                
                m(id,:) = [sb,ROISplit{region},tk{stim, sb},hemo{hb},x,y];
                
                id = id + 1;
            end
        end
    end
end


header_string = {'ID','ROI','Task','Hb','PowerMean','PowerSTD'};
m = cell2table(m,'VariableNames',header_string);
filename = [Data.Path.study_path filesep 'Processed' filesep sub filesep 'Analysis' filesep 'Powerresults_' num2str(F1) '_' num2str(F2) 'Hz.csv'];
writetable(m,filename); % write to csv
clear m

%% Phase
% for each subject for each combination of hb, task, roi take wcoh mean and
% write to csv
 id = 1;
 for hb = 1:numel(hemo)
     for stim = 1:length(si{sb})
         for sb = 1:length(Data.SOI)
             
             n = 2;
             if length(phase) > 1
                 if isempty(phase{2}) %phase only is r1r2 ==> (2)
                    continue
                 end
             else
                 % if the one region we did is r1r2 
                 if ROIs{1} == strjoin([rNames(1) '->' rNames(2)])
                    if isempty(phase{1})
                        continue
                    else
                        n = 1;
                    end
                 end
             end
                 
             taskVector = si{sb}{stim};
             x = nanmean(phase{n}{sb}(hb,taskVector));
             xabs = nanmean(abs(phase{n}{sb}(hb,taskVector)));
             y = nanstd(phase{n}{sb}(hb,taskVector));
             m(id, :) = [sb,"Region 1 -> Region 2",tk{stim, sb},hb,x,xabs,y];
             
         end
     end
 end
 if exist('m')
    header_string = {'ID','Connection','Task','Hb','PhaseMean','AbsPhase','PhaseSTD'};
    m = array2table(m,'VariableNames',header_string);
    filename = [Data.Path.study_path filesep 'Processed' filesep sub filesep 'Analysis' filesep 'Phaseresults_' num2str(F1) '_' num2str(F2) 'Hz.csv'];
    writetable(m,filename); % write to csv
 end
   
 disp('Done Analysis')
 
