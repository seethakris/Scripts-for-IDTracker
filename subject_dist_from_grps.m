
function subject_dist_from_grps(fps, TMin, TMax, tbin, y_thresh_min, y_thresh_max,Num_fish_close_to_subject, time_near)

%% Find average distance of subject when fish from the groups are closer to the subject tank
%% TMin and TMax are optional inputs, if left empty, the entire time is considered

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
%   Num_fish_close_to_subject - the least number of fish that need to be near the subject
%   to include that frame. Default - 3
%   time_near - time (seconds) to find stats on those frames where groups of fish were
%   near the subject for the time specified. Default - 1 second


%% Check for inputs else use default values
pixel_to_mm_change = 3.05; % Using approx. 3.05 pixels/mm


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
    FirstFrame = round(TMin*Frames_per_sec);
else
    FirstFrame = 1;
end

if exist('y_thresh_min') && ~isempty(y_thresh_min)
    Minimum_ythresh = y_thresh_min;
else
    Minimum_ythresh = 0;
end

if exist('y_thresh_max') && ~isempty(y_thresh_max)
    Maximum_ythresh = y_thresh_max;
else
    Maximum_ythresh = 25;
end

if exist('Num_fish_close_to_subject') && ~isempty(Num_fish_close_to_subject)
    Num_fish_close_to_subject;
else
    Num_fish_close_to_subject = 3;
end

if exist('time_near') && ~isempty(time_near)
    Time_threshold = time_near;
else
    Time_threshold = 1;
end

Inputs_provided.Frames_per_sec = Frames_per_sec;
Inputs_provided.frame_bin = frame_bin;
Inputs_provided.FirstFrame = FirstFrame;
Inputs_provided.Minimum_ythresh = Minimum_ythresh;
Inputs_provided.Maximum_ythresh = Maximum_ythresh;
Inputs_provided.Num_fish_close_to_subject = Num_fish_close_to_subject;
Inputs_provided.Time_threshold = Time_threshold;


%% Main Script
PathName = uigetdir(pwd, 'Select modified trajectories file');
FileName = dir([PathName, filesep,'*modified*.mat']);

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
    
    % Collect all group X and Y coordinates
    grp_traj_X = [grp_traj_X; reshape(traj.grp1_XY_mod(:,:,1),size(traj.grp1_XY_mod,1)*size(traj.grp1_XY_mod,2),1); ...
        reshape(traj.grp2_XY_mod(:,:,1),size(traj.grp2_XY_mod,1)*size(traj.grp2_XY_mod,2),1)];
    grp_traj_Y = [grp_traj_Y; reshape(traj.grp1_XY_mod(:,:,2),size(traj.grp1_XY_mod,1)*size(traj.grp1_XY_mod,2),1); ...
        reshape(traj.grp2_XY_mod(:,:,2),size(traj.grp2_XY_mod,1)*size(traj.grp2_XY_mod,2),1)];
    
    grp1_traj_X = [grp1_traj_X; reshape(traj.grp1_XY_mod(:,:,1),size(traj.grp1_XY_mod,1)*size(traj.grp1_XY_mod,2),1)];
    grp2_traj_X = [grp2_traj_X; reshape(traj.grp2_XY_mod(:,:,1),size(traj.grp2_XY_mod,1)*size(traj.grp2_XY_mod,2),1)];
    
    sub_traj_X = [sub_traj_X; traj.subject_XY_mod(:,1,1)];
    sub_traj_Y = [sub_traj_Y; traj.subject_XY_mod(:,1,2)];
end

% Find y min and y max using the grp trajectories
coordinates_y_grp(1) = min(grp_traj_Y) + (max(grp_traj_Y)-min(grp_traj_Y))*(Minimum_ythresh/100);
coordinates_y_grp(2) = min(grp_traj_Y) + (max(grp_traj_Y)-min(grp_traj_Y))*(Maximum_ythresh/100);

coordinates_X_sub = (max(grp1_traj_X)+min(grp2_traj_X))/2;
coordinates_Y_sub = (max(sub_traj_Y)+min(sub_traj_Y))/2;

clear grp_traj_X grp_traj_Y traj grp1_traj_X grp2_traj_X sub_traj_X sub_traj_Y


for ii = 1:length(FileName)
    
    %Save Figures and excel files in these Folder
    SaveName = FileName(ii).name(1:strfind(FileName(1).name, 'modified')-2);
    Result_Folder_excel = [PathName, filesep, 'Excel', filesep, SaveName];
    mkdir(Result_Folder_excel);
    Result_Folder_figures = [PathName, filesep, 'Figures', filesep, SaveName];
    mkdir(Result_Folder_figures);
    Result_Folder_matfiles = [PathName, filesep, 'Matfiles', filesep, SaveName];
    mkdir(Result_Folder_matfiles);
    
    
    disp(['Processing Folder...', SaveName]);
    
    % Load trajectories
    traj = load([PathName, filesep, FileName(ii).name]);
    NumFrames = size(traj.grp1_XY_mod, 1);
    
    % Find Last frame (TMax)
    if exist('TMax') && ~isempty(TMax)
        LastFrame = FirstFrame + frame_bin*fix((TMax*Frames_per_sec-FirstFrame)/frame_bin);
    else
        LastFrame = FirstFrame + frame_bin*fix((NumFrames-FirstFrame)/frame_bin);
    end
    
    Inputs_provided.LastFrame = LastFrame;
    
    % Find frames where groups of fish (2 or more) are within the ROI
    if frame_bin == 1
        
        [Distance, Quadrant_Stats, Time_Threshold] = get_frames_within_ROI(traj, Frames_per_sec, FirstFrame, LastFrame,coordinates_y_grp,...
            coordinates_X_sub,coordinates_Y_sub,Num_fish_close_to_subject,pixel_to_mm_change,Time_threshold);
        
        
        % Save Files
        % Get a name file using all inputs to create unique files for each input
        name_file = ['D_sub_grp_input', '_T=',int2str(round(FirstFrame./Frames_per_sec)), 'to', int2str(round((LastFrame)./Frames_per_sec)),'secs',...
            '_ythresh_', int2str(Inputs_provided.Minimum_ythresh), '%to' , int2str(Inputs_provided.Maximum_ythresh), '%',...
            '_leastfish_', int2str(Inputs_provided.Num_fish_close_to_subject), '_timethreshold_', int2str(Inputs_provided.Time_threshold),'secs'];
        
        % 1. As a matfile
        save([Result_Folder_matfiles, filesep, name_file], 'Distance', 'Quadrant_Stats', 'Inputs_provided');
        
        % 2. As a figure
        %Save trajectories figure
        [fs1] = plot_trajectories(traj,Distance,pixel_to_mm_change);
        name_file = ['Thresholded_Positions_Median', '_T=',int2str(round(FirstFrame./Frames_per_sec)), 'to', int2str(round((LastFrame)./Frames_per_sec)),'secs',...
            '_ythresh_', int2str(Inputs_provided.Minimum_ythresh), '%to' , int2str(Inputs_provided.Maximum_ythresh), '%', ...
            '_leastfish_', int2str(Inputs_provided.Num_fish_close_to_subject), '_timethreshold_', int2str(Inputs_provided.Time_threshold),'secs'];
        
        set(gcf, 'PaperPositionMode','auto','InvertHardCopy', 'off')
        saveas(fs1, [Result_Folder_figures, filesep, name_file], 'tif');
        
        %Save distance we had between subject and groups
        [fs2] = plot_distance_sub_grp(Distance,pixel_to_mm_change);
        name_file = ['Distance_sub_grp_', '_T=',int2str(round(FirstFrame./Frames_per_sec)), 'to', int2str(round((LastFrame)./Frames_per_sec)),'secs',...
            '_ythresh_', int2str(Inputs_provided.Minimum_ythresh), '%to' , int2str(Inputs_provided.Maximum_ythresh), '%', ...
            '_leastfish_', int2str(Inputs_provided.Num_fish_close_to_subject), '_timethreshold_', int2str(Inputs_provided.Time_threshold),'secs'];
        
        set(gcf, 'PaperPositionMode','auto','InvertHardCopy', 'off')
        saveas(fs2, [Result_Folder_figures, filesep, name_file], 'tif');
        
        % 3. As a excel file
        name_file = ['D_sub_grp_input', '_T=',int2str(round(FirstFrame./Frames_per_sec)), 'to', int2str(round((LastFrame)./Frames_per_sec)),'secs',...
            '_ythresh_', int2str(Inputs_provided.Minimum_ythresh), '%to' , int2str(Inputs_provided.Maximum_ythresh), '%', ...
            '_leastfish_', int2str(Inputs_provided.Num_fish_close_to_subject), '_timethreshold_', int2str(Inputs_provided.Time_threshold),'secs'];
        
        save_as_excel(Distance,Quadrant_Stats,Time_Threshold, name_file, Result_Folder_excel);
        
    else
        
        
        for jj = FirstFrame:frame_bin:LastFrame %If the user requires binning of time
            
            [Distance, Quadrant_Stats,Time_Threshold] = get_frames_within_ROI(traj, Frames_per_sec, jj, jj+frame_bin,coordinates_y_grp,...
                coordinates_X_sub, coordinates_Y_sub, Num_fish_close_to_subject,pixel_to_mm_change, Time_threshold);
            
            % Save files
            % 1. As a matfile
            %Get a name file using all inputs to create unique files for each input
            name_file = ['D_sub_grp_input', '_T=',int2str(round(FirstFrame./Frames_per_sec)), 'to', int2str(round((LastFrame)./Frames_per_sec)),'secs',...
                '_ythresh_', int2str(Inputs_provided.Minimum_ythresh), '%to' , int2str(Inputs_provided.Maximum_ythresh), '%', ...
                '_leastfish_', int2str(Inputs_provided.Num_fish_close_to_subject), '_timethreshold_', int2str(Inputs_provided.Time_threshold),'secs'];
            
            save([Result_Folder_matfiles, filesep, name_file], 'Distance', 'Inputs_provided');
            
            % 2. As a figure
            %Save trajectories figure
            [fs1] = plot_trajectories(traj,Distance,pixel_to_mm_change);
            name_file = ['Thresholded_Positions_Median', '_T=',int2str(round(FirstFrame./Frames_per_sec)), 'to', int2str(round((LastFrame)./Frames_per_sec)),'secs',...
                '_ythresh_', int2str(Inputs_provided.Minimum_ythresh), '%to' , int2str(Inputs_provided.Maximum_ythresh), '%', ...
                '_leastfish_', int2str(Inputs_provided.Num_fish_close_to_subject), '_timethreshold_', int2str(Inputs_provided.Time_threshold),'secs'];
            
            set(gcf, 'PaperPositionMode','auto','InvertHardCopy', 'off')
            saveas(fs1, [Result_Folder_figures, filesep, name_file], 'tif');
            
            
            %Save distance we had between subject and groups
            [fs2] = plot_distance_sub_grp(Distance,pixel_to_mm_change);
            name_file = ['Distance_sub_grp_', '_T=',int2str(round(FirstFrame./Frames_per_sec)), 'to', int2str(round((LastFrame)./Frames_per_sec)),'secs',...
                '_ythresh_', int2str(Inputs_provided.Minimum_ythresh), '%to' , int2str(Inputs_provided.Maximum_ythresh), '%', ...
                '_leastfish_', int2str(Inputs_provided.Num_fish_close_to_subject), '_timethreshold_', int2str(Inputs_provided.Time_threshold),'secs'];
            
            set(gcf, 'PaperPositionMode','auto','InvertHardCopy', 'off')
            saveas(fs2, [Result_Folder_figures, filesep, name_file], 'tif');
            
            % 3. As a excel file
            name_file = ['D_sub_grp_input', '_T=',int2str(round(FirstFrame./Frames_per_sec)), 'to', int2str(round((LastFrame)./Frames_per_sec)),'secs',...
                '_ythresh_', int2str(Inputs_provided.Minimum_ythresh), '%to' , int2str(Inputs_provided.Maximum_ythresh), '%', ...
                '_leastfish_', int2str(Inputs_provided.Num_fish_close_to_subject), '_timethreshold_', int2str(Inputs_provided.Time_threshold),'secs'];
            
            save_as_excel(Distance,Quadrant_Stats,Time_Threshold, name_file, Result_Folder_excel);
            
        end
    end
    
end
end




%% Function to find frames where groups of fish are within user given ROI
function [Distance, Quadrant_Stats, Time_Threshold] = get_frames_within_ROI(traj, Frames_per_sec, FirstFrame, LastFrame,coordinates_y_grp, coordinates_X_sub, coordinates_Y_sub, Num_fish,pixel_to_mm_change, Time_threshold)

close all

%% Group 1
% 1. Find frames and x,y position where groups of fish were within specified ROI
Data_grpy = traj.grp1_XY_mod(FirstFrame:LastFrame,:,2);
Data_grp1 = traj.grp1_XY_mod(FirstFrame:LastFrame,:,:);
Frames_grp1 = find((sum(Data_grpy>coordinates_y_grp(1)&Data_grpy<coordinates_y_grp(2),2)>=Num_fish)==1);
Frames_grp1_pos = Data_grp1(Frames_grp1,:,:);

% 2. Get distances from corresponding frames in the subject
Data_subj = traj.subject_XY_mod(FirstFrame:LastFrame,:,:);
Frames_subject_grp1_pos = Data_subj(Frames_grp1,:,:);

% 3. Get center of mass by taking median of x,y positions of group
Centre_mass_grp1 = squeeze(median(Frames_grp1_pos,2));

% 4.Take euclidean distances between centre of mass of group and subject
Dist_subject_from_grp1_pos = sqrt((squeeze(Frames_subject_grp1_pos(:,1,1))-Centre_mass_grp1(:,1)).^2 + ...
    (squeeze(Frames_subject_grp1_pos(:,1,2))-Centre_mass_grp1(:,2)).^2 );

% 5. Convert pixels to mm
Dist_subject_from_grp1_pos_mm = Dist_subject_from_grp1_pos/pixel_to_mm_change;
Centre_mass_grp1_mm = Centre_mass_grp1/pixel_to_mm_change;

% 6. Find quadrants where the subjects are located. There are 4 quadrants.
Quadrant_subject_grp1_pos = zeros(size(Frames_subject_grp1_pos,1),1);

Quadrant_subject_grp1_pos(Frames_subject_grp1_pos(:,1)<coordinates_X_sub & Frames_subject_grp1_pos(:,2)>coordinates_Y_sub) = 1; %Group1 close
Quadrant_subject_grp1_pos(Frames_subject_grp1_pos(:,1)>coordinates_X_sub & Frames_subject_grp1_pos(:,2)>coordinates_Y_sub) = 2; %Group2 close

Quadrant_subject_grp1_pos(Frames_subject_grp1_pos(:,1)<coordinates_X_sub & Frames_subject_grp1_pos(:,2)<coordinates_Y_sub) = 3; %Group1 far
Quadrant_subject_grp1_pos(Frames_subject_grp1_pos(:,1)>coordinates_X_sub & Frames_subject_grp1_pos(:,2)<coordinates_Y_sub) = 4; %Group2 far

% 7. Find frames where groups of fish were near subject for the amount
% of time specified by user (in seconds)
Frames_sec_grp1 = fix(Frames_grp1./Frames_per_sec);
count = 1;
temp_seconds = unique(Frames_sec_grp1);

Time_found_grp1 = 0;
Time_centremass_grp1 = 0;
Time_subdistance_grp1 = 0;
Time_quadrant_grp1 = 0;

for ii = 1:length(temp_seconds)
    if size(find(ismember(Frames_sec_grp1, temp_seconds(ii):temp_seconds(ii)+Time_threshold-1)==1),1)>=Frames_per_sec*Time_threshold
        Time_found_grp1(count,1) = temp_seconds(ii);
        Time_centremass_grp1(count,1) = mean(Centre_mass_grp1_mm(ismember(Frames_sec_grp1, temp_seconds(ii):temp_seconds(ii)+Time_threshold-1)==1));
        Time_subdistance_grp1(count,1) = mean(Dist_subject_from_grp1_pos_mm(ismember(Frames_sec_grp1, temp_seconds(ii):temp_seconds(ii)+Time_threshold-1)==1));
        Time_quadrant_grp1(count,1) = fix(median(Quadrant_subject_grp1_pos(ismember(Frames_sec_grp1, temp_seconds(ii):temp_seconds(ii)+Time_threshold-1)==1)));
        count = count+1;
    end
end

%% Group2
% 1. Find frames and x,y position where groups of fish were within specified ROI
Data_grpy = traj.grp2_XY_mod(FirstFrame:LastFrame,:,2);
Data_grp2 = traj.grp2_XY_mod(FirstFrame:LastFrame,:,:);
Frames_grp2 = find((sum(Data_grpy>coordinates_y_grp(1)&Data_grpy<coordinates_y_grp(2),2)>=Num_fish)==1);
Frames_grp2_pos = Data_grp2(Frames_grp2,:,:);

% 2. Get distances from corresponding frames in the subject
Frames_subject_grp2_pos = Data_subj(Frames_grp2,:,:);

% 3. Get center of mass by taking median of x,y positions of group
Centre_mass_grp2 = squeeze(median(Frames_grp2_pos,2));

% 4.Take euclidean distances between centre of mass of group and subject
Dist_subject_from_grp2_pos = sqrt((squeeze(Frames_subject_grp2_pos(:,1,1))-Centre_mass_grp2(:,1)).^2 + ...
    (squeeze(Frames_subject_grp2_pos(:,1,2))-Centre_mass_grp2(:,2)).^2 );

% 5. Convert pixels to mm
Dist_subject_from_grp2_pos_mm = Dist_subject_from_grp2_pos/pixel_to_mm_change;
Centre_mass_grp2_mm = Centre_mass_grp2/pixel_to_mm_change;


% 6. Find quadrants where the subjects are located. There are 4 quadrants.
Quadrant_subject_grp2_pos = zeros(size(Frames_subject_grp2_pos,1),1);

Quadrant_subject_grp2_pos(Frames_subject_grp2_pos(:,1)<coordinates_X_sub & Frames_subject_grp2_pos(:,2)>coordinates_Y_sub) = 1; %Group1 close
Quadrant_subject_grp2_pos(Frames_subject_grp2_pos(:,1)>coordinates_X_sub & Frames_subject_grp2_pos(:,2)>coordinates_Y_sub) = 2; %Group2 close

Quadrant_subject_grp2_pos(Frames_subject_grp2_pos(:,1)<coordinates_X_sub & Frames_subject_grp2_pos(:,2)<coordinates_Y_sub) = 3; %Group1 far
Quadrant_subject_grp2_pos(Frames_subject_grp2_pos(:,1)>coordinates_X_sub & Frames_subject_grp2_pos(:,2)<coordinates_Y_sub) = 4; %Group2 far


% 7. Find frames where groups of fish were near subject for the amount
% of time specified by user (in seconds)
Frames_sec_grp2 = fix(Frames_grp2./Frames_per_sec);
count = 1;
temp_seconds = unique(Frames_sec_grp2);

Time_found_grp2 = 0;
Time_centremass_grp2 = 0;
Time_subdistance_grp2 = 0;
Time_quadrant_grp2 = 0;

for ii = 1:length(temp_seconds)
    if size(find(ismember(Frames_sec_grp2, temp_seconds(ii):temp_seconds(ii)+Time_threshold-1)==1),1)>=Frames_per_sec*Time_threshold
        Time_found_grp2(count,1) = temp_seconds(ii);
        Time_centremass_grp2(count,1) = mean(Centre_mass_grp2_mm(ismember(Frames_sec_grp2, temp_seconds(ii):temp_seconds(ii)+Time_threshold-1)==1));
        Time_subdistance_grp2(count,1) = mean(Dist_subject_from_grp2_pos_mm(ismember(Frames_sec_grp2, temp_seconds(ii):temp_seconds(ii)+Time_threshold-1)==1));
        Time_quadrant_grp2(count,1) = fix(median(Quadrant_subject_grp2_pos(ismember(Frames_sec_grp2, temp_seconds(ii):temp_seconds(ii)+Time_threshold-1)==1)));
        count = count+1;
    end
end

% 8. Find frames where both groups are at the ROI
Frames_both_grp1 = ismember(Frames_grp1, Frames_grp2); %0-only group1, 1-both group1 and 2
Frames_both_grp2 = ismember(Frames_grp2, Frames_grp1); %0-only group2, 1-both group1 and 2


% 9. Seperate stats based on Quadrant
for ii = 1:4
    Quadrant_num(ii,1) = ii;
    Quadrant_count_grp1(ii,1) = size(find(Quadrant_subject_grp1_pos==ii),1);
    Quadrant_Centremass_grp1(ii,1) = mean(Centre_mass_grp1(Quadrant_subject_grp1_pos==ii)/pixel_to_mm_change);
    Quadrant_subdist_grp1(ii,1) = mean(Dist_subject_from_grp1_pos_mm(Quadrant_subject_grp1_pos==ii));
    
    Quadrant_count_grp2(ii,1) = size(find(Quadrant_subject_grp2_pos==ii),1);
    Quadrant_Centremass_grp2(ii,1) = mean(Centre_mass_grp2(Quadrant_subject_grp2_pos==ii)/pixel_to_mm_change);
    Quadrant_subdist_grp2(ii,1) = mean(Dist_subject_from_grp2_pos_mm(Quadrant_subject_grp2_pos==ii));
end



%% Save as Matfile
Distance.Frames_grp1 = Frames_grp1;
Distance.Frames_grp1_pos = Frames_grp1_pos;
Distance.Frames_subject_grp1_pos = Frames_subject_grp1_pos;
Distance.Centre_mass_grp1 = Centre_mass_grp1;
Distance.Centre_mass_grp1_mm = Centre_mass_grp1_mm;
Distance.Dist_subject_from_grp1_pos = Dist_subject_from_grp1_pos;
Distance.Dist_subject_from_grp1_pos_mm = Dist_subject_from_grp1_pos_mm;
Distance.Frames_both_grp1 = Frames_both_grp1;
Distance.Quadrant_subject_grp1_pos = Quadrant_subject_grp1_pos;

Distance.Frames_grp2 = Frames_grp2;
Distance.Frames_grp2_pos = Frames_grp2_pos;
Distance.Frames_subject_grp2_pos = Frames_subject_grp2_pos;
Distance.Centre_mass_grp2 = Centre_mass_grp2;
Distance.Centre_mass_grp2_mm = Centre_mass_grp2_mm;
Distance.Dist_subject_from_grp2_pos = Dist_subject_from_grp2_pos;
Distance.Dist_subject_from_grp2_pos_mm = Dist_subject_from_grp2_pos_mm;
Distance.Frames_both_grp2 = Frames_both_grp2;
Distance.Quadrant_subject_grp2_pos = Quadrant_subject_grp2_pos;

Quadrant_Stats.Quadrant_num = Quadrant_num;
Quadrant_Stats.Quadrant_count_grp1 = Quadrant_count_grp1;
Quadrant_Stats.Quadrant_Centremass_grp1 = Quadrant_Centremass_grp1;
Quadrant_Stats.Quadrant_subdist_grp1 = Quadrant_subdist_grp1;
Quadrant_Stats.Quadrant_count_grp2 = Quadrant_count_grp2;
Quadrant_Stats.Quadrant_Centremass_grp2 = Quadrant_Centremass_grp2;
Quadrant_Stats.Quadrant_subdist_grp2 = Quadrant_subdist_grp2;

Time_Threshold.Time_found_grp1 = Time_found_grp1;
Time_Threshold.Time_centremass_grp1 = Time_centremass_grp1;
Time_Threshold.Time_subdistance_grp1 = Time_subdistance_grp1;
Time_Threshold.Time_quadrant_grp1 = Time_quadrant_grp1;
Time_Threshold.Time_found_grp2 = Time_found_grp2;
Time_Threshold.Time_centremass_grp2 = Time_centremass_grp2;
Time_Threshold.Time_subdistance_grp2 = Time_subdistance_grp2;
Time_Threshold.Time_quadrant_grp2 = Time_quadrant_grp2;

end



%% Save as Excel
function save_as_excel(Distance, Quadrant_Stats, Time_Threshold, name_file, Result_Folder_excel)

% 1. First save distances with names suitable for excel files

Distance1.Group1_frames = Distance.Frames_grp1;
Distance1.Group1_centre_mass_mm = Distance.Centre_mass_grp1_mm;
Distance1.Group1_subjectdist_mm = Distance.Dist_subject_from_grp1_pos_mm;
Distance1.Group1_and_group2_frames = Distance.Frames_both_grp1;
Distance1.Group1_subject_quadrant = Distance.Quadrant_subject_grp1_pos;

group1_matrices = 5; % to know how many Group2 vs Group1 quantifications - provide a gap while saving in excel

Distance1.Group2_frames = Distance.Frames_grp2;
Distance1.Group2_centre_mass_mm = Distance.Centre_mass_grp2_mm;
Distance1.Group2_subjectdist_mm = Distance.Dist_subject_from_grp2_pos_mm;
Distance1.Group2_and_group1_frames = Distance.Frames_both_grp2;
Distance1.Group2_subject_quadrant = Distance.Quadrant_subject_grp2_pos;

group_matrices = 10; % to know how many group quantifications - provide a gap while saving in excel

Distance1.Quadrant = Quadrant_Stats.Quadrant_num;
Distance1.Count_grp1 = Quadrant_Stats.Quadrant_count_grp1;
Distance1.Centremass_grp1 = Quadrant_Stats.Quadrant_Centremass_grp1;
Distance1.Subdist_grp1 = Quadrant_Stats.Quadrant_subdist_grp1;
Distance1.Count_grp2 = Quadrant_Stats.Quadrant_count_grp2;
Distance1.Centremass_grp2 = Quadrant_Stats.Quadrant_Centremass_grp2;
Distance1.Subdist_grp2 = Quadrant_Stats.Quadrant_subdist_grp2;

quad_matrices = 17;

Distance1.Group1_time_secs = Time_Threshold.Time_found_grp1;
Distance1.Group1_time_centre_mass_mm = Time_Threshold.Time_centremass_grp1;
Distance1.Group1_time_subjectdist_mm = Time_Threshold.Time_subdistance_grp1;
Distance1.Group1_time_subject_quadrant = Time_Threshold.Time_quadrant_grp1;

group_time_matrices = 21;

Distance1.Group2_time_secs = Time_Threshold.Time_found_grp2;
Distance1.Group2_time_centre_mass_mm = Time_Threshold.Time_centremass_grp2;
Distance1.Group2_time_subjectdist_mm = Time_Threshold.Time_subdistance_grp2;
Distance1.Group2_time_subject_quadrant = Time_Threshold.Time_quadrant_grp2;

clear Distance Quadrant_Stats

Temp_Dat = fieldnames(Distance1);
filename = [ Result_Folder_excel,filesep,name_file,'.xls'];
fid = fopen(filename, 'w+');

% Go through and save as a cell in a format suitable for excel files
count = 0;
for kk = 1:length(Temp_Dat)
    if kk == group1_matrices+1 || kk == group_matrices+1 || kk == quad_matrices+1 || kk == group_time_matrices+1
        count = count+3;
    else
        count = count+1;
    end
    
    Xls_Dat{1,count} = Temp_Dat{kk};
    for ii = 1:size(Distance1.(Temp_Dat{kk}),1)
        temp1 = Distance1.(Temp_Dat{kk})(ii);
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

%% Plotting:
function [fs1] = plot_trajectories(traj, Distance,pixel_to_mm_change)

% Plot trajectories in chosen frames for user and save as
% Thresholded_Positions Distance_sub_grp and in Figures folder

fs1 = figure(1);
set(fs1,'color','white')
hold on

plot(squeeze(Distance.Frames_grp1_pos(:,:,1)), squeeze(Distance.Frames_grp1_pos(:,:,2)),'.')
plot(Distance.Centre_mass_grp1(:,1), Distance.Centre_mass_grp1(:,2), 'k+', 'MarkerSize', 10)
plot(squeeze(Distance.Frames_grp2_pos(:,:,1)), squeeze(Distance.Frames_grp2_pos(:,:,2)),'*')
plot(Distance.Centre_mass_grp2(:,1), Distance.Centre_mass_grp2(:,2), 'k+', 'MarkerSize', 10)
plot(squeeze(Distance.Frames_subject_grp1_pos(:,:,1)), squeeze(Distance.Frames_subject_grp1_pos(:,:,2)),'r.', 'MarkerSize', 10)
plot(squeeze(Distance.Frames_subject_grp2_pos(:,:,1)), squeeze(Distance.Frames_subject_grp2_pos(:,:,2)),'g*')

hold off
set(gca, 'TickDir','out', 'FontSize',12)
set(gca, 'YDir', 'reverse');
box off

%Convert to mm
x = get(gca, 'Xtick');
set(gca, 'xTickLabel',strread(int2str(x/pixel_to_mm_change),'%s'));
y = get(gca, 'Ytick');
set(gca, 'yTickLabel',strread(int2str(y/pixel_to_mm_change),'%s'));
xlabel(gca,'x distance (mm)', 'FontSize',12);
ylabel(gca,'y distance (mm)', 'FontSize',12);
title(gca, 'Positions of fish in selected ROI')
set(gcf,'position',get(0,'screensize'))

%Create legend string
for ii = 1:size(traj.grp1_XY_mod, 2)
    legend_str1{ii} = ['Group1 Fish',int2str(ii)];
end
legend_str1{size(traj.grp1_XY_mod, 2)+1} = 'Center of Mass Group1';

for ii = size(traj.grp1_XY_mod, 2)+2:size(traj.grp1_XY_mod, 2)+size(traj.grp2_XY_mod, 2)+1
    legend_str1{ii} = ['Group2 Fish',int2str(ii)];
end
legend_str1{size(traj.grp1_XY_mod, 2)+size(traj.grp2_XY_mod, 2)+2} = 'Center of Mass Group2';

legend_str1{size(traj.grp1_XY_mod, 2)+size(traj.grp2_XY_mod, 2)+3} = 'Subject and Group1';
legend_str1{size(traj.grp1_XY_mod, 2)+size(traj.grp2_XY_mod, 2)+4} = 'Subject and Group2';

legend(legend_str1{:});

end



function [fs2] = plot_distance_sub_grp(Distance,pixel_to_mm_change)
% 2. Plot distance from centre of mass of groups and subject for selected
% frames

fs2 = figure(2);
set(fs2,'color','white')
subplot(2,1,1)
plot(Distance.Frames_grp1, Distance.Dist_subject_from_grp1_pos_mm, 'r.')
hold on
plot(Distance.Frames_grp2, Distance.Dist_subject_from_grp2_pos_mm, 'g.')
hold off

set(gca, 'TickDir','out', 'FontSize',12)
box off
xlabel(gca,'Frame Number', 'FontSize',12);
ylabel(gca,'distance (mm)', 'FontSize',12);
title(gca, 'Distance of subject from groups of fish', 'FontSize',12)

legend('Subject and Group1', 'Subject and Group2');

subplot(2,1,2)
[mean1, ci1] = normfit(Distance.Dist_subject_from_grp1_pos_mm);
errorbar(1,mean1,ci1,'r')
hold on
[mean2, ci2] = normfit(Distance.Dist_subject_from_grp2_pos_mm);
errorbar(2,mean2,ci2,'g')
hold off

set(gca, 'TickDir','out', 'FontSize',12)
box off
xlabel(gca,'Frame', 'FontSize',12);
ylabel(gca,'distance (mm)', 'FontSize',12);
title(gca, 'Distance of subject from groups of fish', 'FontSize',12)

legend('Subject and Group1', 'Subject and Group2');

end


