function group_distance_from_subject(fps, TMin, TMax, tbin, y_thresh_min, y_thresh_max)

%% Find which fish were close to the subject when they are in specific quadrants

close all
warning off

%% Input :
%% Optional Inputs :
% You can use [] to use default: Default values can be changed below as
% indicated in the script.
%
%   fps - frames per second. Default - 30
%   TMin - Minimum time (seconds). Default - First Frame in video
%   TMax - Maximum time (seconds). Default - Last Frame in video
%   tbin - bin data over specified time bin (seconds). Default - one frame
%   y_thresh_min and y_thresh_max - ROI over which to analyze group
%   behavior (in percent). Default - 0% (border closest to subject fish) to 25%

%% Check for inputs else use default values
pixel_to_mm_change = 3.05; % Using approx. 3.05 pixels/mm

% If user specified inputs, else default
if exist('fps','var') && ~isempty(fps)
    Frames_per_sec = fps;
else
    Frames_per_sec = 30;
end

if exist('tbin','var') && ~isempty(tbin)
    frame_bin = tbin*Frames_per_sec;
else
    frame_bin = 1;
end

if exist('TMin','var') && ~isempty(TMin)
    FirstFrame = round(TMin*Frames_per_sec);
else
    FirstFrame = 1;
end

if exist('y_thresh_min','var') && ~isempty(y_thresh_min)
    Minimum_ythresh = y_thresh_min;
else
    Minimum_ythresh = 0;
end

if exist('y_thresh_max','var') && ~isempty(y_thresh_max)
    Maximum_ythresh = y_thresh_max;
else
    Maximum_ythresh = 25;
end

Inputs_provided.Frames_per_sec = Frames_per_sec;
Inputs_provided.frame_bin = frame_bin;
Inputs_provided.FirstFrame = FirstFrame;
Inputs_provided.Minimum_ythresh = Minimum_ythresh;
Inputs_provided.Maximum_ythresh = Maximum_ythresh;

%% Main Script
PathName = uigetdir(pwd, 'Select modified trajectories file');
FileName = dir([PathName, filesep,'*modified*.mat']);

if isempty(FileName)
    return;
end

% Find y threshold using data of all fish to find the most accurate boundary.
grp_traj_X = [];
grp_traj_Y = [];
sub_traj_X = [];
sub_traj_Y = [];
grp1_traj_X = [];
grp2_traj_X = [];

for ii = 1:length(FileName)
    % Load trajectories
    traj = load([PathName, filesep, FileName(ii).name]);
    
    %Check which group is present
    if size(traj.grp2_XY_mod,2) == 0
        disp('Only Group 1 fish exist');
        grp_XY_mod = traj.grp1_XY_mod;
        flag = 1;
        % Collect all group X and Y coordinates
        grp_traj_X = [grp_traj_X; reshape(grp_XY_mod(:,:,1),size(grp_XY_mod,1)*size(grp_XY_mod,2),1)];
        grp_traj_Y = [grp_traj_Y; reshape(grp_XY_mod(:,:,2),size(grp_XY_mod,1)*size(grp_XY_mod,2),1)];
        
        sub_traj_X = [sub_traj_X; traj.subject_XY_mod(:,1,1)];
        sub_traj_Y = [sub_traj_Y; traj.subject_XY_mod(:,1,2)];
        
    elseif size(traj.grp1_XY_mod,2) == 0
        disp('Only Group 2 fish exist');
        grp_XY_mod = traj.grp2_XY_mod;
        flag = 1;
        % Collect all group X and Y coordinates
        grp_traj_X = [grp_traj_X; reshape(grp_XY_mod(:,:,1),size(grp_XY_mod,1)*size(grp_XY_mod,2),1)];
        grp_traj_Y = [grp_traj_Y; reshape(grp_XY_mod(:,:,2),size(grp_XY_mod,1)*size(grp_XY_mod,2),1)];
        
        sub_traj_X = [sub_traj_X; traj.subject_XY_mod(:,1,1)];
        sub_traj_Y = [sub_traj_Y; traj.subject_XY_mod(:,1,2)];
        
    else
        % Collect all group X and Y coordinates
        flag = 0;
        grp_traj_X = [grp_traj_X; reshape(traj.grp1_XY_mod(:,:,1),size(traj.grp1_XY_mod,1)*size(traj.grp1_XY_mod,2),1); ...
            reshape(traj.grp2_XY_mod(:,:,1),size(traj.grp2_XY_mod,1)*size(traj.grp2_XY_mod,2),1)];
        grp_traj_Y = [grp_traj_Y; reshape(traj.grp1_XY_mod(:,:,2),size(traj.grp1_XY_mod,1)*size(traj.grp1_XY_mod,2),1); ...
            reshape(traj.grp2_XY_mod(:,:,2),size(traj.grp2_XY_mod,1)*size(traj.grp2_XY_mod,2),1)];
        
        grp1_traj_X = [grp1_traj_X; reshape(traj.grp1_XY_mod(:,:,1),size(traj.grp1_XY_mod,1)*size(traj.grp1_XY_mod,2),1)];
        grp2_traj_X = [grp2_traj_X; reshape(traj.grp2_XY_mod(:,:,1),size(traj.grp2_XY_mod,1)*size(traj.grp2_XY_mod,2),1)];
        
        sub_traj_X = [sub_traj_X; traj.subject_XY_mod(:,1,1)];
        sub_traj_Y = [sub_traj_Y; traj.subject_XY_mod(:,1,2)];
        
    end
end

if flag == 0 %Both groups are present
    
    % Find y min and y max using the grp trajectories
    coordinates_y_grp(1) = min(grp_traj_Y) + (max(grp_traj_Y)-min(grp_traj_Y))*(Minimum_ythresh/100);
    coordinates_y_grp(2) = min(grp_traj_Y) + (max(grp_traj_Y)-min(grp_traj_Y))*(Maximum_ythresh/100);
    
    coordinates_X_sub = (max(grp1_traj_X)+min(grp2_traj_X))/2;
    coordinates_Y_sub = (max(sub_traj_Y)+min(sub_traj_Y))/2;
else %Only a single group is present
    
    % Find y min and y max using the grp trajectories
    coordinates_y_grp(1) = min(grp_traj_Y) + (max(grp_traj_Y)-min(grp_traj_Y))*(Minimum_ythresh/100);
    coordinates_y_grp(2) = min(grp_traj_Y) + (max(grp_traj_Y)-min(grp_traj_Y))*(Maximum_ythresh/100);
    
    coordinates_X_sub = (max(sub_traj_X)+min(sub_traj_X))/2;
    coordinates_Y_sub = (max(sub_traj_Y)+min(sub_traj_Y))/2;
end

clear grp_traj_X grp_traj_Y sub_traj_X sub_traj_Y traj grp_XY_mod traj grp1_traj_X grp2_traj_X


for ii = 1:length(FileName)
    
    %Save Figures and excel files in these Folder
    SaveName = FileName(ii).name(1:strfind(FileName(1).name, 'modified')-2);
    Result_Folder_excel = [PathName, filesep, 'Excel', filesep, SaveName];
    mkdir(Result_Folder_excel);
    Result_Folder_matfiles = [PathName, filesep, 'Matfiles', filesep, SaveName];
    mkdir(Result_Folder_matfiles);
    
    disp(['Processing Folder...', SaveName]);
    
    % Load trajectories
    traj = load([PathName, filesep, FileName(ii).name]);
    NumFrames = size(traj.grp1_XY_mod, 1);
    
    %Get Fish Numbers
    FishNum = load([PathName, filesep,SaveName, '_FishNumber.mat']);
    
    % Find Last frame (TMax)
    if exist('TMax','var') && ~isempty(TMax)
        LastFrame = FirstFrame + frame_bin*fix((TMax*Frames_per_sec-FirstFrame)/frame_bin);
    else
        LastFrame = FirstFrame + frame_bin*fix((NumFrames-FirstFrame)/frame_bin);
    end
    
    
    %Check which group is present
    if size(traj.grp2_XY_mod,2) == 0 %Only Group1 is present
        Inputs_provided.grp_string = 'Group1';
        Quadrant = get_frames_within_ROI(traj.grp1_XY_mod, traj.subject_XY_mod, Inputs_provided.grp_string,Frames_per_sec,...
            FirstFrame, LastFrame,coordinates_y_grp,coordinates_X_sub,coordinates_Y_sub,...
            Inputs_provided.Minimum_ythresh, Inputs_provided.Maximum_ythresh,SaveName, ...
            Result_Folder_matfiles, Result_Folder_excel,FishNum.FishNumber.grp1);
        
    elseif size(traj.grp1_XY_mod,2) == 0 %Only Group2 is present
        Inputs_provided.grp_string = 'Group2';
        Quadrant = get_frames_within_ROI(traj.grp2_XY_mod, traj.subject_XY_mod, Inputs_provided.grp_string,Frames_per_sec,...
            FirstFrame, LastFrame,coordinates_y_grp,coordinates_X_sub,coordinates_Y_sub,...
            Inputs_provided.Minimum_ythresh, Inputs_provided.Maximum_ythresh,SaveName, ...
            Result_Folder_matfiles, Result_Folder_excel,FishNum.FishNumber.grp2);
        
    else %Both Group1 and Group2 are present
        grp_string = 'Group1';
        Quadrant = get_frames_within_ROI(traj.grp1_XY_mod, traj.subject_XY_mod, grp_string,Frames_per_sec,...
            FirstFrame, LastFrame,coordinates_y_grp,coordinates_X_sub,coordinates_Y_sub,...
            Inputs_provided.Minimum_ythresh, Inputs_provided.Maximum_ythresh,SaveName, ...
            Result_Folder_matfiles, Result_Folder_excel,FishNum.FishNumber.grp1);
        
        grp_string = 'Group2';
        Quadrant = get_frames_within_ROI(traj.grp2_XY_mod, traj.subject_XY_mod, grp_string,Frames_per_sec,...
            FirstFrame, LastFrame,coordinates_y_grp,coordinates_X_sub,coordinates_Y_sub,...
            Inputs_provided.Minimum_ythresh, Inputs_provided.Maximum_ythresh,SaveName, ...
            Result_Folder_matfiles, Result_Folder_excel,FishNum.FishNumber.grp2);
    end
    
end

end

%% Function to find frames where groups were when subject was ina  specific quadrant

function Quadrant = get_frames_within_ROI(grp_XY_mod, subject_XY_mod,grp_string,Frames_per_sec,...
    FirstFrame, LastFrame,coordinates_y_grp,coordinates_X_sub,coordinates_Y_sub,...
    Minimum_ythresh, Maximum_ythresh, SaveName, Result_Folder_matfiles,Result_Folder_excel,FishNum)

% 1. Get distances from corresponding frames in the subject
Data_subj = subject_XY_mod(FirstFrame:LastFrame,:,:);

% 2. Find quadrants where the subjects are located. There are 4 quadrants.
Quadrant_subject = zeros(size(Data_subj,1),1);

Quadrant_subject(Data_subj(:,1,1)<coordinates_X_sub & Data_subj(:,1,2)>coordinates_Y_sub) = 1; %Group1 close
Quadrant_subject(Data_subj(:,1,1)>coordinates_X_sub & Data_subj(:,1,2)>coordinates_Y_sub) = 2; %Group2 close

Quadrant_subject(Data_subj(:,1,1)<coordinates_X_sub & Data_subj(:,1,2)<coordinates_Y_sub) = 3; %Group1 far
Quadrant_subject(Data_subj(:,1,1)>coordinates_X_sub & Data_subj(:,1,2)<coordinates_Y_sub) = 4; %Group2 far


% 3. Find frames where the subject were in a particular quadrant for 1
% second. Also Find group details based on that.
temp_quad = unique(fix(Quadrant_subject));
Quadrant.QuadNum = temp_quad;

for ii = 1:length(temp_quad)
    % Find frames where the subject is in a quadrant for one second
    ind_quad = fix(find(Quadrant_subject==temp_quad(ii))./Frames_per_sec);
    ind_quad_unique = unique(ind_quad);
    ind_quad_within_timethresh = ind_quad_unique((histc(ind_quad, ind_quad_unique)>=Frames_per_sec-4));
    
    % Get trajectories of groups fish and find if they are near the subject
    % in these frames
    for jj = 1:length(ind_quad_within_timethresh)
        Data_grpy = median(squeeze(grp_XY_mod(ismember(ind_quad, ind_quad_within_timethresh(jj)),:,2)));
        Frames_grp(jj,:) = Data_grpy>coordinates_y_grp(1)&Data_grpy<coordinates_y_grp(2);
    end
    
    Quadrant.(['Quad',int2str(temp_quad(ii)), '_Time']) = ind_quad_within_timethresh;
    
    for kk = 1:length(FishNum)
        Quadrant.(['Quad',int2str(temp_quad(ii)), '_Fishid', int2str(FishNum(kk))]) = Frames_grp(:,kk);
        Quadrant.(['Fishid', int2str(FishNum(kk))])(temp_quad(ii),1) = sum(Frames_grp(:,kk));
    end
    clear Frames_grp
end


% Save files
%Get a name file using all inputs to create unique files for each input
name_file = [SaveName, '_',grp_string,'_from_S_T=',int2str(round(FirstFrame./Frames_per_sec)), 'to', int2str(round(LastFrame./Frames_per_sec)),'secs',...
    '_ythresh_', int2str(Minimum_ythresh), '%to' , int2str(Maximum_ythresh), '%'];

% 1. As a matfile
save([Result_Folder_matfiles, filesep, name_file], 'Quadrant');

% 2. As a excel file
save_as_excel(Quadrant,name_file, Result_Folder_excel, FishNum)

end

%% Save as Excel
function save_as_excel(Quadrant,name_file, Result_Folder_excel, FishNum)

%Arrange structure to save
for ii = 1:length(Quadrant.QuadNum)
    for jj = 1:length(FishNum)
        Quadrant1.(['Quad', int2str(Quadrant.QuadNum(ii)),'Fishid', int2str(FishNum(jj))]) = Quadrant.(['Quad',int2str(Quadrant.QuadNum(ii)), '_Fishid', int2str(FishNum(jj))]);
    end
    count_quad(ii) = (length(FishNum))*ii;
end

Quadrant1.Quadrant = Quadrant.QuadNum;

for jj = 1:length(FishNum)
    Quadrant1.(['Fishid', int2str(FishNum(jj))]) = Quadrant.(['Fishid', int2str(FishNum(jj))]);
end

Temp_Dat = fieldnames(Quadrant1);
filename = [ Result_Folder_excel,filesep,name_file,'.xls'];
fid = fopen(filename, 'w+');

% Go through and save as a cell in a format suitable for excel files
count = 0;
count_quad = count_quad+1;
for kk = 1:length(Temp_Dat)
    if find((ismember(count_quad, kk) == 1))
        count = count+3;
    else
        count = count+1;
    end
    
    Xls_Dat{1,count} = Temp_Dat{kk};
    for ii = 1:size(Quadrant1.(Temp_Dat{kk}),1)
        temp1 = Quadrant1.(Temp_Dat{kk})(ii);
        if temp1 ~= 0
            Xls_Dat{ii+1,count} = temp1;
        else
            Xls_Dat{ii+1,count} = 0;
        end
    end
    
end

%Save the cell as excel
fid = fopen(filename, 'a');
[nrows,ncols]= size(Xls_Dat);

for row = 1:nrows
    for col = 1:ncols
        if row == 1
            fprintf(fid, '%s\t', Xls_Dat{row,col});
        else
            fprintf(fid, '%4.2f\t', Xls_Dat{row,col});
        end
    end
    fprintf(fid, '\n');
end

end


