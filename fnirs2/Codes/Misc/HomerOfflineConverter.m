function HomerOfflineConverter(files)
%This function converts recorded NIRx datasets into Homer2 format (*.nirs)
%Please provide as input the root path, within all datasets are located
%If no probeInfo file can be found, it is necessary to choose it manually
%
%Note: to account for the original inter-optode distances, a few functions
%are added to properly import the 2D coordinates and convert them to 3D
%
%By: NIRx Medical Technologies
%Contact: support@nirx.net
%
%Last Update: 2019-02-12 - JWP

%INPUT: files = ['filePathForSubject1', 'filePathForSubject2'...]
% will be better for me to just write my own script to get the file names
% since we only need a few at a time in most cases. Also add something to
% check if the converted one exists before converting.
for subj=1:length(files)
    disp('converting file for subj:')
    disp(subj)
    wb = waitbar(0,'Converting dataset into Homer2 format... Please wait.');
    set(wb,'Name','Converter');

    dataset = files(subj);
    
    [pathstr, name] = fileparts(dataset);

    config_filename = [pathstr filesep name '_config.txt'];
    fp =  fopen(config_filename); 

    while (1)

        sline = fgets(fp);
        if sline == -1
            break;
        end

        try
            if contains(sline, 'SamplingRate')
                sr = sscanf(sline, 'SamplingRate=%f\n');
            elseif contains(sline, 'waveLength_N')
                WLs = sscanf(sline, 'waveLength_N=%d\n');
            elseif contains(sline, 'Wavelengths')
                if WLs == 2
                    lambda = sscanf(sline, 'Wavelengths=[%d %d]', [1 2]);
                else
                    lambda = sscanf(sline, 'Wavelengths=[%d	%d %d %d]', [1 WLs]);
                end
            elseif contains(sline, 'source_N')
                nsrc_cfg = sscanf(sline, 'source_N=%d\n');
            elseif contains(sline, 'detector_N')
                ndet_cfg = sscanf(sline, 'detector_N=%d\n');
            end

        catch
            disp(['Error loading ' files(subj).name]);
            break;
        end

    end

    fclose(fp);

    crtFile1=[pathstr filesep name '_probeInfo.mat'];
    if exist(crtFile1, 'file') ~= 2
        disp(['Error: probeInfo file in folder ' pathstr ' not found.']);
        [filename, pathname] = uigetfile('*_probeInfo.mat', 'Pick probeInfo file');
        crtFile1 = [pathname filename];
    end    
    coms = ['load ' crtFile1 '  probeInfo'];
    eval(coms);

    nsrc = probeInfo.probes.nSource0;
    ndet = probeInfo.probes.nDetector0;

    if (nsrc > nsrc_cfg) || (ndet > ndet_cfg)
        flag = 0;
    else
        flag = 1;
    end

    while flag == 0    

        uiwait(msgbox('The selected probeInfo file does not correspond to the saved dataset. Please choose a new one.','modal'));
        [filename1, pathname1] = uigetfile('*_probeInfo.mat', 'Pick new _probeInfo file');

        %in case user cancels the probeInfo loading, end loop
        if isequal(filename1,0) || isequal(pathname1,0)
            uiwait(msgbox('The process of formatting the data for Homer2 has been cancelled.'));
            return;
        end

        probe_path = [pathname1, filename1];
        comms = ['load ' probe_path '  probeInfo'];
        eval(comms); %loading new probeInfo file   

        nsrc = probeInfo.probes.nSource0;
        ndet = probeInfo.probes.nDetector0;

        if (nsrc <= nsrc_cfg) || (ndet <= ndet_cfg)
            flag = 1;
        end

    end        

    %SD.MeasList field
    nchan = probeInfo.probes.nChannel0;
    demo.SD.MeasList = zeros([WLs*nchan 4]); %The number of columns is always 4, as explained inthe Homer user's guide.
    demo.SD.MeasList(:,3) = 1; %The Homer user's guide explains that this column is presently unused.
    for i=1:WLs
        demo.SD.MeasList(1+nchan*(i-1):i*nchan,4) = i; %Initializing rows to contain data for wavelength i
    end
    demo.SD.MeasList(:,1:2) = repmat(probeInfo.probes.index_c,[WLs 1]); % Modified so that MeasList is not initialized to contain 4 wavelength data

    %SD.Lambda field
    demo.SD.Lambda = lambda; %wavelengths information

    newcoords = rotate_clusters(probeInfo); %cluster-based rotation of 3D coords

    %SD.SrcPos and SD.DetPos fields
    demo.SD.SrcPos = zeros([nsrc 3]); 
    demo.SD.DetPos = zeros([ndet 3]);

    %check how the probeInfo was created (NIRSite or nirsLAB)
    if isfield(probeInfo, 'headmodel') && strcmp(probeInfo.headmodel, 'ICBM152')
        fprintf('ProbeInfo file is using the ICBM-152 head model\n');
        demo.SD.SrcPos(:,1:2) = newcoords(1:nsrc,1:2);    
        demo.SD.DetPos(:,1:2) = newcoords(nsrc+1:end,1:2);
    else
        fprintf('ProbeInfo file is using the nirsLAB head model\n');
        demo.SD.SrcPos(:,1:2) = -newcoords(1:nsrc,1:2);     %additional 180deg rotation
        demo.SD.DetPos(:,1:2) = -newcoords(nsrc+1:end,1:2); %additional 180deg rotation
    end    

    demo.SD.SrcPos(:,3) = newcoords(1:nsrc,3); %src coords
    demo.SD.DetPos(:,3) = newcoords(nsrc+1:end,3); %det coords

    %SD.nSrcs and SD.nDets fields
    demo.SD.nSrcs = nsrc; demo.SD.nDets = ndet; demo.SD.SpatialUnit = 'cm';
    demo.SD.xmin = 0; %xmin, xmax, ymin, ymax are initialized here
    demo.SD.xmax = 0; %only so the field MeasListAct will be added
    demo.SD.ymin = 0; %on the expected position -> this is required
    demo.SD.ymax = 0; %to avoid conflict when loading data on Homer2
    demo.SD.MeasListAct = ones(size(demo.SD.MeasList,1),1);
    demo.SD.MeasListVis = ones(size(demo.SD.MeasList,1),1);
    %SD sructure array completed.

    %initialize struct by loading the last wavelength
    wl_ini{WLs} = load([pathstr filesep name '.wl' num2str(WLs)]); 
    for i=1:(WLs-1)
        wl_ini{i} = load([pathstr filesep name '.wl' num2str(i)]);
    end

    tot = nsrc_cfg*ndet_cfg;
    channels = probeInfo.probes.index_c;

    %if masked saved data -> adapt dataset
    if size(wl_ini{WLs},2) < tot 
        wl{WLs} = zeros(size(wl_ini{WLs},1),tot);
        for i=1:(WLs-1)
            wl{i} = wl{WLs};
        end
        for i=1:WLs
            for j=1:size(channels,1)
                wl{i}(:,(ndet_cfg*(channels(j,1)-1))+channels(j,2)) = wl_ini{i}(:,j);
            end
        end
    %if data saved as expected -> copy to 'wl'
    else 
        wl = wl_ini;   
    end

    %'t' (time) field of the *.nirs structure array
    nt = size(wl{WLs},1); %number of frames recorded
    demo.t = (1:nt)'/sr;
    %'d' (data) field of the *.nirs structure array
    demo.d = zeros([nt WLs*nchan]);

    %check for saturated data and mark corresponding channel as 'bad'
    for i = 1:WLs
        for k1 = 1:nchan
            demo.d(:,k1+((i-1)*nchan)) = wl{i}(:,ndet_cfg*(demo.SD.MeasList(k1,1) - 1) + demo.SD.MeasList(k1,2));
            if nonzeros(isnan(demo.d(:,k1+((i-1)*nchan)))) %look for NaN values
                demo.SD.MeasListAct(k1+((i-1)*nchan)) = 0; %mark channel as bad
            end
        end
    end

    %'s' (event triggers) field of the *.nirs structure array
    evt = load([pathstr filesep name '.evt']);
    if size(evt,1) ~= 0    
        trgaux = repmat([1 2 4 8 16 32 64 128],size(evt,1),1); %binary entries
        markers = sum(evt(:,2:end).*trgaux,2); %triggers in decimal unit
        evt = [evt(:,1) markers];
        uc = 1:max(markers);
        uc = uc';
    else
        uc = [];
        markers = 1;
    end
%     uc = unique(evt(:,2));
    if size(uc,1) ~= 0
        nc = size(uc,1); %This is the number of distinct conditions.
    else
        nc = 1;
    end
    demo.s = zeros([nt min(max(markers),1)]);
%     demo.s = zeros([nt nc]); 
    for k2 = 1:nc
        for k3 = 1:size(evt,1)
            if evt(k3,1) ~= 0
                if evt(k3,2:end) == uc(k2,:) %Find the time frames that have start triggers for the k2-th condition.
                    demo.s(evt(k3,1),k2) = 1; %Flag the correponding rows in the k2-th column of the demo.s field.
                end
            end
        end
    end

    %'aux' (auxiliary) field of the *.nirs structure array. 
    %This is optional in terms of using it during data processing with Homer. 
    %But Homer will not load the *.nirs file if this field doesn't exist.
    demo.aux = ones(size(demo.t));
    
    %Save current dataset into *.nirs file
    save([pathstr filesep name '.nirs'],'-struct','demo');
    disp(['Dataset converted: ' files(subj).name]);
    close(wb)
    clear dataset nsrc ndet nchan probeInfo flag demo newcoords wl_ini wl

end

function newcoords = rotate_clusters(probeInfo)

    src = probeInfo.probes.coords_s3;
    det = probeInfo.probes.coords_d3;
    channels = probeInfo.probes.index_c;
    
    pi.nsources = length(src);
    pi.ndetectors = length(det);
    
    pi.optode_coords = [src; det]; %stack src and det 3D coords
    
    pi.channel_indices = zeros(length(channels),2);
    pi.channel_distances = zeros(length(channels),1);
    
    for i=1:length(channels)
        src_i = channels(i,1);
        det_i = channels(i,2);
        pi.channel_indices(i,:) = [src_i (length(src)+det_i)];
        pi.channel_distances(i) = norm(src(src_i,:) - det(det_i,:));
    end
    
    origin = find_origin(pi);
    for i=1:length(pi.optode_coords)
        pi.optode_coords(i,:) = pi.optode_coords(i,:) - origin;
    end
    
    clusters = cluster_search_mat(pi);
    
    newcoords = zeros(size(pi.optode_coords));
    
    for i = 1:size(clusters,2)
        
        idx = clusters{i};
        center = mean(pi.optode_coords(idx',:));
        
        %center in spherical coordinates
        center_r = norm(center);
        center_phi = atan2(center(2), center(1));
        center_theta = acos(center(3)/center_r);
        
        %phi tangent vector
        tangent = [-sin(center_phi) cos(center_phi) 0];
        
        %angle between r and z unit vectors
        angle = acos((center/norm(center))*[0 0 1]');
        
        mat = rotmat(center, tangent, -angle);
        
        coords = ( mat * [pi.optode_coords(idx',:) ones(length(idx),1)]' )';
        
        newcoords(idx',:) = coords(:,1:3);
        
    end

    
function origin = find_origin(pi)
    
    origin = fminsearch(@(x) fun(x, pi), [0 0 0]);
    
    function f = fun(x, pi)
        
        k = ones(length(pi.optode_coords),1);
        for i=1:length(k)
            k(i) = norm(pi.optode_coords(i,:) - x);
        end
        f = std(k);
     
        
function found = cluster_search_mat(pi)

    index = pi.channel_indices;
    k = index(1,1);
    source = 1;
    cluster = 1;
    found{cluster} = k;
    found_src(1) = k;
    found_det = [];
    
    while size(index,1) > 0
        if source
            chn = find(index(:,1) == k); %channels with source = k
            if ~isempty(chn)
                for i=1:size(chn)
                    if ~any(found{cluster} == index(chn(i),2))
                        found{cluster} = [found{cluster} index(chn(i),2)];
                        found_det = [found_det index(chn(i),2)'];
                    end
                end
                index(chn,:) = []; %remove channels from index list
            end
            if isempty(index)
                break;
            end
            found_src(1) = []; %remove current source index
            if isempty(found_src) %change to detector indexes
                source = 0;
                if isempty(found_det)
                    k = index(1,2); %if both are empty, re-initialize
                    found_det(1) = k;
                    cluster = cluster + 1; %go to next cluster
                    found{cluster} = k;
                else
                    k = found_det(1);
                end
            else
                k = found_src(1);
            end
        else
            chn = find(index(:,2) == k); %channels with detector = k
            if ~isempty(chn)
                for i=1:size(chn)
                    if ~any(found{cluster} == index(chn(i),1))
                        found{cluster} = [found{cluster} index(chn(i),1)];
                        found_src = [found_src index(chn(i),1)];
                    end
                end
                index(chn,:) = []; %remove channels from index list
            end
            if isempty(index)
                break;
            end
            found_det(1) = []; %remove current source index
            if isempty(found_det) %change to detector indexes
                source = 1;
                if isempty(found_src)
                    k = index(1,1); %if both are empty, re-initialize
                    found_src(1) = k;
                    cluster = cluster + 1; %go to next cluster
                    found{cluster} = k;
                else
                    k = found_src(1);
                end
            else
                k = found_det(1);
            end
        end
    end


function mat = rotmat(point, direction, theta)

    a = point(1); b = point(2); c = point(3);
    
    t = direction/norm(direction);
    u = t(1); v = t(2); w = t(3);

    si = sin(theta);
    co = cos(theta);

    mat = zeros(4);

    % rotational part    
    mat(1:3, 1:3) = [ (u*u + (v*v + w*w) * co) (u*v*(1-co) - w*si)     (u*w*(1-co) + v*si);
                      (u*v*(1-co) + w*si)      (v*v + (u*u + w*w)*co)  (v*w*(1-co) - u*si);
                      (u*w*(1-co) - v*si)      (v*w*(1-co) + u*si)     (w*w + (u*u + v*v)*co) ];

    % translational part
    mat(1,4) = (a*(v*v+w*w)-u*(b*v+c*w)) * (1-co) + (b*w-c*v)*si;
    mat(2,4) = (b*(u*u+w*w)-v*(a*u+c*w)) * (1-co) + (c*u-a*w)*si;
    mat(3,4) = (c*(u*u+v*v)-w*(a*u+b*v)) * (1-co) + (a*v-b*u)*si;
    mat(4,4) = 1;