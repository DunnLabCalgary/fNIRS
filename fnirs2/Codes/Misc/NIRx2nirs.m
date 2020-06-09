function NIRx2nirs(NIRx_foldername, SD_filename)

% NIRx2nirs.m version 1.0
% #####################################################
% This script takes a folder (NIRx_foldername) containing the NIRx output
% data files (.hdr, .wl1, .wl2) and a pre-defined SD file (SD_filename) 
% (built using the Homer2 SDgui), which matches the source-detector layout 
% used in the NIRx acquisition and creates a .nirs file for use in Homer2

% To use this script, the user must first create an SD file which matches
% their NIRx probe layour using the SDgui function of Homer2.  It is 
% essential that the SD file loaded matches the NIRx acquisition layout as 
% this is assumed to be correct by this script.  This includes maintaining 
% the real-world NIRx source and detector numbers in the SD file, which may
% necessitate padding the SD file if consecutively numbered sources and
% detectors, starting from 1, were not used.

% This code was written and debugged using data from the NIRx NIRSCOUT, it
% may not be applicable to other models.

% Rob J Cooper, University College London, August 2013
% robert.cooper@ucl.ac.uk

% #########################################################################

% Select NIRx folder containing .wl1, .wl2 and .hr files
if ~exist('NIRx_foldername','var');
NIRx_foldername = uigetdir(pwd,'Select NIRx Data Folder...');
end

% Load SD_file
if ~exist('SD_filename','var');
    SD_filename = uigetfile('.SD','Select associated SD file...');
end
load(SD_filename,'-mat');

% Load wavelength d
% #######################################################################
wl1_dir = dir([NIRx_foldername '/*.wl1']);
if length(wl1_dir) == 0; 
    error('ERROR: Cannot find NIRx .wl1 file in selected directory...'); 
    disp('Please ensure there are NIRx output data files (.hdr, .wl1, .wl2) in the Raw Data directory'); end;
wl1 = load([NIRx_foldername '/' wl1_dir(1).name]);
wl2_dir = dir([NIRx_foldername '/*.wl2']);
if length(wl2_dir) == 0; 
    error('ERROR: Cannot find NIRx .wl2 file in selected directory...');
    disp('Please ensure there are NIRx output data files (.hdr, .wl1, .wl2) in the Raw Data directory'); end;
wl2 = load([NIRx_foldername '/' wl2_dir(1).name]);

d=[wl1 wl2]; % d matrix from .wl1 and .wl2 files

% Read and interpret .hdr d ############################################
% #########################################################################
hdr_dir = dir([NIRx_foldername '/*.hdr']);
if length(hdr_dir) == 0; error('ERROR: Cannot find NIRx header file in selected directory...'); end;
fid = fopen([NIRx_foldername '/' hdr_dir(1).name]);
tmp = textscan(fid,'%s','delimiter','\n');%This just reads every line
hdr_str = tmp{1};
fclose(fid);

%Find number of sources
keyword = 'Sources=';
tmp = strfind(hdr_str,keyword);
ind = find(~cellfun(@isempty,tmp)); %This gives cell of hdr_str with keyword
tmp = hdr_str{ind};
NIRx_Sources = str2num(tmp(length(keyword)+1:end));

%Find number of sources
keyword = 'Detectors=';
tmp = strfind(hdr_str,keyword);
ind = find(~cellfun(@isempty,tmp)); %This gives cell of hdr_str with keyword
tmp = hdr_str{ind};
NIRx_Detectors = str2num(tmp(length(keyword)+1:end));

%Compare to SD file for checking...
if NIRx_Sources ~= SD.nSrcs || NIRx_Detectors ~= SD.nDets;
   error('The number or sources and detectors in the NIRx files does not match your SD file...');
end

%Find Sample rate
keyword = 'SamplingRate=';
tmp = strfind(hdr_str,keyword);
ind = find(~cellfun(@isempty,tmp)); %This gives cell of hdr_str with keyword
tmp = hdr_str{ind};
fs = str2num(tmp(length(keyword)+1:end));

%Find Active Source-Detector pairs (these will just be ordered by source,
%then detector (so, for example d(:,1) = source 1, det 1 and d(:,2) =
%source 1 det 2 etc.
keyword = 'S-D-Mask="#';
tmp = strfind(hdr_str,keyword);
ind = find(~cellfun(@isempty,tmp)) + 1; %This gives cell of hdr_str with keyword
tmp = strfind(hdr_str(ind+1:end),'#');
ind2 = find(~cellfun(@isempty,tmp)) - 1;
ind2 = ind + ind2(1);
sd_ind = cell2mat(cellfun(@str2num,hdr_str(ind:ind2),'UniformOutput',0));
sd_ind = sd_ind';
sd_ind = find([sd_ind(:);sd_ind(:)]);
d = d(:,sd_ind);

%Find Event Markers and build S vector
keyword = 'Events="#';
tmp = strfind(hdr_str,keyword);
ind = find(~cellfun(@isempty,tmp)) + 1; %This gives cell of hdr_str with keyword
tmp = strfind(hdr_str(ind+1:end),'#');
ind2 = find(~cellfun(@isempty,tmp)) - 1;
ind2 = ind + ind2(1);
events = cell2mat(cellfun(@str2num,hdr_str(ind:ind2),'UniformOutput',0));
events = events(:,2:3);
markertypes = unique(events(:,1));
s = zeros(length(d),length(markertypes));
for i = 1:length(markertypes);
    s(events(find(events(:,1)==markertypes(i)),2),i) = 1;
end

%Create t, aux varibles
aux = zeros(length(d),8);
t = 0:1/fs:length(d)/fs - 1/fs;

outname = uiputfile('*.nirs','Save .nirs file ...');

fprintf('Saving as %s ...\n',outname);
save(outname,'d','s','t','aux','SD');
