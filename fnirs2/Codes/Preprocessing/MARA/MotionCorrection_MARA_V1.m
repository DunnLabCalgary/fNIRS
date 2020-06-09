function [dodmc] = MotionCorrection_MARA_V1(dod,fs,mldq)

% This function is allows for easy user input to correct channels of
% interest for motion artifact using the MARA_1_1 function. User can
% choose only those channels suspected of motion artifact, based on plot
% view of all channels



% hb = Data.ui.hemo{HB};
% d = Data.Subjects.(sub).Data.dc;
% fs = 1/(t(2)-t(1));
% Ch_List = ;



ch = inputdlg('Which channels would you like to apply motion correction to?                        ie... 1 2 3 4 or 1:4');
ch = str2num(ch{:});
close all;

% INPUT:
%   x:              Signal (one dimensional vector)
%   T:              Threshold for artifact detection
%   L:              Length of the moving-window for moving std
%                   calculation
%   p:              Window length for artifact correction (i.e. length of the LOESS smoothing window)

%%% Here we start loop through selected channels
for ii = 1:length(ch)
    k = ch(ii);kk = ['Channel ' num2str(k)];
    %     idx = find(Ch_List == k);
    idx = k;
    % Defaults
    
%     T = 1*10^-6; % default set high for typical concentration scale
    L = round(fs*(1/2)); % default set to 1 seconds
    p = round(fs*(1/2)); % default set to 0.10 seconds
    yay = movingstd(dod(:,idx),L);
    T = mean(yay)+(3*std(yay));
%     T = max(yay)-(0.1*mean(yay))
    
    pp = {num2str(T),num2str(L),num2str(p)};
    
    myflag = true;
    
    while myflag % while loop allows user to adjust parameters before moving on to next channel
        pp = inputdlg({'Enter T (threshold):',...
            'Enter L (window length of std detection):',...
            'Enter p (window length for artifact correction):'},...
            'MARA Input',...
            1,pp);
        
        try
            [dodmc(:,idx)] = MARA_1_1(dod(:,idx),fs,...
                str2double(pp{1}),str2double(pp{2}),str2double(pp{3}),kk);
            
        catch
            display('Error: Re-adjust the threshold level') % if error occurs, readjust parameters, typically threshold is too far off
            pause(2)
            
        end
        pause(2)
        answer = input(['1: Continue' newline ...
        '2: Adjust Parameters' newline ...
        '3: No Correction' newline ...
        '4: Correct, then re-adjust parameters' newline newline ...
        'Press a key____']) % answer here determines if you exit while loop and move on to next channel
        
        if answer == 1
            myflag = false;     
        elseif answer == 3
            dodmc(:,idx) = dod(:,idx);
            myflag = false;  
            
        elseif answer == 4
            dod(:,idx) = dodmc(:,idx);
        end

        
    end
    
    
end

% save([sub 'processed'],'Pre')


end
