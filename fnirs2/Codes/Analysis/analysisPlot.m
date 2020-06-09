function [X] = analysisPlot(si, dc, Wavelet,fs, R1Ch, R2Ch, Hemo, FZband, t, xt, filePath, label, rNames, analysisFigures, stimNames)
    
    
    %selected channels during the stim)
    X.fftcoh = fftcoherence(si, dc, fs, R1Ch, R2Ch, Hemo, FZband); %calculate fourier transform coherence
    
    X.tvcoh = plotTVcoh(R1Ch, R2Ch, Wavelet.wcoh, t, xt, Hemo, label, analysisFigures(1));
    if analysisFigures(1)
        saveas(gcf,[filePath 'tvCoh.png']);
        close
    end
    
    if ~isequal(R1Ch, R2Ch)
        %phase makes no sense for intra regional.
        X.tvphase = plotTVphase(R1Ch, R2Ch, Wavelet.phase, t, xt, Hemo, label, analysisFigures(2));
        if analysisFigures(2)
            saveas(gcf,[filePath 'tvPhase.png']);
            close
        end
    end
    
    X.tvwt = plotTVwt(R1Ch, R2Ch, Wavelet.wt, t, xt, Hemo, label, analysisFigures(3), rNames, FZband);
    if analysisFigures(3)
        saveas(gcf,[filePath 'wt.png']);
        close
    end
    
    X.fzcoh = plotFvCohByTask(Wavelet.wcohf, Wavelet.wcohfzRange, Hemo, R1Ch, R2Ch, label, analysisFigures(4), stimNames);
    if analysisFigures(4)
        saveas(gcf,[filePath 'fzByTask.png']);
        close
    end
    
    %% Move to new function for study specific analysis
    % [wcohX,wtX,phaseX] = tvGroupPlot(Data,F1,F2,hemo,fs,s_all,wcoh,wt,phase); 
    % clearvars wt wcoh phase
    % [wt,wcoh,phase] = tvMean(Data,F1,F2,wcohX,wtX,phaseX);
    % [wt,wcoh,phase] = tvSTD(Data,F1,F2,wcohX,wtX,phaseX,wt,wcoh,phase);
    % corrplots(Data,F1,F2,wt,wcoh,phase)
    
end
    
        function ym = fftcoherence(si, d, fs, R1Ch, R2Ch, hemo, FZband)
            % do fft coherence to compare to wavelet coherence
            % Cohfft output is a channel by channel matrix of static coherence
            % for each task (3rd dimension)

    %         stimIndeces = find(s);
            for hb = 1:numel(hemo) %for each hemo param

                 for stim = 1:length(si)
                     %d = (data points*hemo*channels) so this is selecting all
                     %the data points between stim start stimIndeces(stim) and stim end
                     %stimIndeces(stim+1), resulting size is (#data points * 1 * #channels)
                     data = d(si{stim},hb,:);
                     %get rid of singleton dim so we have matrix of (data pts * channels)
                     reducedData = squeeze(data);

                     %compute fftCoh
                     [Coh(:,:,stim)] = fftCoh2(reducedData,fs,min(FZband),max(FZband));

                 end

                 y = Coh(R1Ch, R2Ch,:);
                 %ym1 = (1 * numR2Channels * nStims) where 1 is the mean along
                 %the R1Channels
                 ym1 = nanmean(y);
                 % ym2 = (1*1*8) so its the mean of coh during each sitm
                 ym2 = nanmean(ym1);
                 %get rid of the singleton dimensions
                 ym(hb,:) = squeeze(ym2);
                 
            end

        end   
      
        function all = plotTVcoh(R1Ch, R2Ch, d, t, xt, hemo, label, shouldPlot)
        
            if shouldPlot
                Rect = [0.05, 0.05, 0.93, 0.9];
                h2 = figure('Name',['Time-varying Coherence ' label],...
                'Color','w',...
                'NumberTitle','off',...
                'units','normalized',...
                'outerposition',[0 0 1 1],...
                'DefaultaxesPosition',Rect);
                ha2 = tight_subplot(2,2,...
                [.08 .05],[.05 .04],[.04 .01]);
                cl = get(gca,'colororder');
                
            end

            for hb = 1:numel(hemo)
                %get hemo string
                HBString = hemo{hb};
                %select axes

                y = d.(HBString);
                
                ych = y(R1Ch, R2Ch, :);
                y = reshape(ych,[],numel(t));
                if length(R1Ch) > 1 && length(R2Ch) > 1
                    ym = nanmean(y)';
                else
                    ym = y;
                end
                
                
                if shouldPlot
                    axes(ha2(hb));
                    hold on;

                    plot(t,ym, 'Color',cl(hb,:),'LineWidth',2);

                    xticklabels('auto')
                    yticklabels('auto')
                    title([HBString ' Coherence ' label])
                    ylabel('Coherence')
                    xlabel('time (t)')
                    set(gca,'FontSize',15,'FontWeight','bold')
                    for nn = 1:numel(xt)
                        plot([xt(nn) xt(nn)],ylim,'--k')
                    end
                end

                all(hb,:) = ym;

            end
        end
        
        function all = plotTVphase(R1Ch, R2Ch ,d,t,xt,hemo, label, shouldPlot)
                
            if shouldPlot
                Rect = [0.05, 0.05, 0.93, 0.9];
                h2 = figure('Name',['Time-varying Coherence ' label],...
                    'Color','w',...
                    'NumberTitle','off',...
                    'units','normalized',...
                    'outerposition',[0 0 1 1],...
                    'DefaultaxesPosition',Rect);
                ha2 = tight_subplot(2,2,...
                    [.08 .05],[.05 .04],[.04 .01]);
                cl = get(gca,'colororder');
            end
                for hb = 1:numel(hemo)
                    HBString = hemo{hb};
                    
                    y = d.(HBString);
                    y = reshape(y(R1Ch, R2Ch, :),[],numel(t));
                    
                    if length(R1Ch) > 1 && length(R2Ch) > 1
                        ym = nanmean(y)';
                    else
                        ym = y;
                    end
                   
                    ys = nanstd(y)';
                    lo = ym - ys;
                    hi = ym + ys;
                    
                    if shouldPlot
                        axes(ha2(hb))
                        hold on;
                        h(hb) = plot(t,ym,'Color',cl(hb,:),'LineWidth',2);

                        xticklabels('auto')
                        yticklabels('auto')
                        title([HBString ' Phase ' label])
                        ylabel('Radians')
                        xlabel('time (t)')
                        set(gca,'FontSize',15,'FontWeight','bold')
                        for nn = 1:numel(xt)
                            plot([xt(nn) xt(nn)],ylim,'--k')
                        end
                    
                    end
                    
                    all(hb,:) = ym;

                end
        end 
        
        function all = plotTVwt(R1Ch, R2Ch, d, t, xt, hemo, label, shouldPlot, rNames, FZband)
            if shouldPlot 
                Rect = [0.05, 0.05, 0.93, 0.9];
                h2 = figure('Name',['Time-varying Power ' label num2str(FZband(1)) '-' num2str(FZband(2))],...
                    'Color','w',...
                    'NumberTitle','off',...
                    'units','normalized',...
                    'outerposition',[0 0 1 1],...
                    'DefaultaxesPosition',Rect);
                ha2 = tight_subplot(2,2,...
                    [.08 .05],[.05 .04],[.04 .01]);

                cl = get(gca,'colororder');
            end

            channels = {};
        
            if ~isequal(R1Ch, R2Ch) %then we have an inter
                channels{1} = R1Ch;
                channels{2} = R2Ch;
            else
                channels{1} = R1Ch;
            end
            
            for hb = 1:numel(hemo)
                HBString = hemo{hb};
                
                if shouldPlot
                    axes(ha2(hb))
                end

                for j = 1:length(channels) 

                    y = d.(HBString);

                    y = y(channels{j},:);
                    
                    if length(R1Ch) > 1 && length(R2Ch) > 1
                        ym = nanmean(y)';
                    else 
                        ym = y;
                    end
                    
                    ys = nanstd(y)';
                    lo = ym - ys;
                    hi = ym + ys;
                    
                    if shouldPlot
                        hold on;
                        h(j) = plot(t,ym,'Color',cl(hb,:),'LineWidth',2);
                        
                        xticklabels('auto')
                        yticklabels('auto')
                        title([HBString ' Power ' label])
                        ylabel('Power')
                        xlabel('time (t)')
                        set(gca,'FontSize',15,'FontWeight','bold');
                       
                    end

                    all(j, hb,:) = ym;

                end
                
                if shouldPlot
                    
                    for nn = 1:numel(xt)
                        plot([xt(nn) xt(nn)],ylim,'--k')
                    end
                    
                    legend(h, rNames)
                end
                
                
            end
        end
        
        function all = plotFvCohByTask(d, fzRange, hemo, R1Ch, R2Ch, label, shouldPlot, stimNames)

            if shouldPlot
                Rect = [0.05, 0.05, 0.93, 0.9];
                h2 = figure('Name',['Frequency-varying Coherence ' label],...
                    'Color','w',...
                    'NumberTitle','off',...
                    'units','normalized',...
                    'outerposition',[0 0 1 1],...
                    'DefaultaxesPosition',Rect);
                ha2 = tight_subplot(2,2,...
                    [0.12 .05],[.05 .04],[.04 .01]);
            end

            for hb = 1:length(hemo)

                HBString = hemo{hb};
                
                if shouldPlot
                    axes(ha2(hb));
                end

                data = d.(HBString);

                for i = 1:length(data)
                     y = data(i).stim;
                     fz = fzRange(i).fz;

                     ych = y(R1Ch, R2Ch, :);
                     sz = size(ych);
                     y = reshape(ych, [], sz(3));
                     y(y==0)=NaN;
                     
                     if length(R1Ch) > 1 && length(R2Ch) > 1 
                         ym = nanmean(y)';
                     else
                         ym = y;
                     end
                     
                     ys = nanstd(y)';

                    fzcut = fz < 0.3; 
                    fz = fz(fzcut);
                    ym = ym(fzcut);
                    
                    if shouldPlot
                        plot(log2(fz),ym,'LineWidth',2) 
                        legend(stimNames)
                        yticklabels('auto')
                        title([HBString ' Coherence ' label])
                        xlabel('frequency (Hz)')
                        Xticks = [0.01,0.02,0.03,0.06,0.12,0.25,0.50,1.00];
                        set(ha2(hb),'FontSize',15,'FontWeight','bold','XTick',log2(Xticks(:)), ...
                           'XTickLabel',num2str(sprintf('%.2f\n',Xticks)));

                        hold on;
                    end
                    
                    all{hb, i} = ym;
                end

            end

            hold off;
        end
           
%          function [all] = plotFVwt(ROI,d,hemo)
% 
%                  for hb = 1:numel(hemo)
%                      HB = hemo{hb};
%                      ct = 1;
%                      ctt = 1;
%                      for r1 = 1:numel(ROI.Names)
%                          roi = str2num(ROI.Channels{r1});
%                          y = d.(HB);
%                          ym = squeeze(nanmean(y));
%                          all(ctt,hb,:,:) = ym;
%                          ctt = ctt+1;
%                      end
%                      ct = ct+1;
% 
%                  end
%          end
           
    %      function [all] = plotFVphase(ROI,d,hemo)
    %     
    %          for hb = 1:numel(hemo)
    %              HB = hemo{hb};
    %              ct = 1;
    %              ctt = 1;
    %              ac =[1,1];
    %              for r1 = 1:numel(ROI.Names)
    %                  roi1 = str2num(ROI.Channels{r1});
    %                  for r2 = 1:numel(ROI.Names)
    %                      if r1 ~= r2 & ~ismember(sort([r1,r2]),sort(ac,2),'rows')
    %                          roi2 = str2num(ROI.Channels{r2});
    %                          y = d.(HB);
    %                          y = reshape(y(roi1,roi2,:),[],size(y,3),size(y,4));
    %                          ym = squeeze(nanmean(y));
    %                          all(ctt,hb,:,:) = ym;
    %                          ctt = ctt+1;
    %                      end
    %                      ac(ct,:) = [r1,r2];
    %                      ct = ct+1;
    %                  end
    %              end    
    
    %          end
    
    

% % % % % % 

