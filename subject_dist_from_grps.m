function subject_dist_from_grps(fps, TMin, TMax, tbin, y_thresh)

%% Find average distance of subject when fish from the groups are closer to the subject tank
%% TMin and TMax are optional inputs, if left empty, the entire time is considered

%% Input :
%% Optional Inputs :
% You can use [] to use default: Default values can be changed below as
% indicated in the script
%
%   fps - frames per second. Default - 30
%   TMin - Minimum time (seconds). Default - First Frame in video
%   TMax - Maximum time (seconds). Default - Last Frame in video
%   tbin - bin data over specified time bin (seconds). Default - one frame
%   y_thresh - y below which groups of fish are considered closer to
%   subject. Default = 600

%% Check for inputs else use default values

% If user specified inputs, else default
if exist('fps') && ~isempty(fps)
    Frames_per_sec = fps;
else
    Frames_per_sec = 30;
end

if exist('tbin') && ~isempty(tbin)
    frame_bin = tbin*Frames_per_sec;
else
    frame_bin = 1;
end

if exist('TMin') && ~isempty(TMin)
    FirstFrame = round(TMin*fps);
else
    FirstFrame = 1;
end

if exist('TMax') && ~isempty(TMax)
    LastFrame = round(TMax*fps);
else
    LastFrame = round(TMax*fps);
end



if exist('TMax') && ~isempty(TMax)
    LastFrame = round(TMax*Frames_per_sec);
else
    LastFrame = NumFrames;
end





%% Main Script
PathName = uigetdir(pwd, 'Select modified trajectories file');
FileName = dir([PathName, filesep,'*modified*.mat']);


%Approximate area that should be occupied by groups of fish.
%[xleft,xright,ytop,ybottom];
Group1_area =[150, 650, 950, 450];
Group2_area =[650, 1150, 950, 450];
Subject_area =[150, 1150, 450, 200];

for ii = 1:length(FileName)
    %Save Figures in This Folder
    SaveName = FileName(ii).name(1:strfind(FileName(1).name, 'modified')-2);
    Result_Folder = [PathName, filesep, 'ExcelFiles', filesep, SaveName];
    mkdir(Result_Folder);
    
    disp(['Processing Folder...', SaveName]);
    
    
    % Load trajectories
    traj = load([PathName, filesep, FileName(ii).name]);
    NumFrames = size(traj.grp1_XY_mod, 1);
    
    
    
    
    for jj = FirstFrame:LastFrame
        
    end
    
end
