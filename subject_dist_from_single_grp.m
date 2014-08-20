
function subject_dist_from_single_grp(fps, TMin, TMax, tbin, y_thresh_min, y_thresh_max,Num_fish_close_to_subject, time_near)

%% Find average distance of subject when fish from the groups are closer to the subject tank
%% In this case there is only one group of fish

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

if exist('Num_fish_close_to_subject','var') && ~isempty(Num_fish_close_to_subject)
    Num_fish_close_to_subject;
else
    Num_fish_close_to_subject = 3;
end

if exist('time_near','var') && ~isempty(time_near)
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

if isempty(FileName)
    return;
end

% Find y threshold using data of all fish to find the most accurate boundary.
grp_traj_X = [];
grp_traj_Y = [];
sub_traj_X = [];
sub_traj_Y = [];

%Check which group exists and only use that one group for gurther analysis
for ii = 1:length(FileName)    
    % Load trajectories
    traj = load([PathName, filesep, FileName(ii).name]);
    
    %Check which group is present
    if size(traj.grp2_XY_mod,2) == 0
        disp('Only Group 1 fish exist');
        grp_XY_mod = traj.grp1_XY_mod;
        
    elseif size(traj.grp1_XY_mod,2) == 0
        disp('Only Group 2 fish exist');
        grp_XY_mod = traj.grp2_XY_mod;
        
    end
    
    % Collect all group X and Y coordinates
    grp_traj_X = [grp_traj_X; reshape(grp_XY_mod(:,:,1),size(grp_XY_mod,1)*size(grp_XY_mod,2),1)];
    grp_traj_Y = [grp_traj_Y; reshape(grp_XY_mod(:,:,2),size(grp_XY_mod,1)*size(grp_XY_mod,2),1)];
    
    sub_traj_X = [sub_traj_X; traj.subject_XY_mod(:,1,1)];
    sub_traj_Y = [sub_traj_Y; traj.subject_XY_mod(:,1,2)];    
end

%Plot
fs = figure(1);
plot(grp_traj_X, grp_traj_Y, 'b')
hold on
plot(sub_traj_X, sub_traj_Y, 'r')
hold off
set(gca, 'TickDir','out', 'FontSize',12)
set(gca, 'YDir', 'reverse');
box off

%Check with the user if trajectories are OK
button = questdlg('Do you want to keep these experiments in the same folder?', 'Verify trajectories', 'OK', 'Cancel','OK');
if strcmp(button,'Cancel')
    disp('Canceled file operation')
    return;
end

% Find y min and y max using the grp trajectories
coordinates_y_grp(1) = min(grp_traj_Y) + (max(grp_traj_Y)-min(grp_traj_Y))*(Minimum_ythresh/100);
coordinates_y_grp(2) = min(grp_traj_Y) + (max(grp_traj_Y)-min(grp_traj_Y))*(Maximum_ythresh/100);

coordinates_X_sub = (max(sub_traj_X)+min(sub_traj_X))/2;
coordinates_Y_sub = (max(sub_traj_Y)+min(sub_traj_Y))/2;

clear grp_traj_X grp_traj_Y sub_traj_X sub_traj_Y traj grp_XY_mod traj

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
    if exist('TMax','var') && ~isempty(TMax)
        LastFrame = FirstFrame + frame_bin*fix((TMax*Frames_per_sec-FirstFrame)/frame_bin);
    else
        LastFrame = FirstFrame + frame_bin*fix((NumFrames-FirstFrame)/frame_bin);
    end
    
    Inputs_provided.LastFrame = LastFrame;
    
    % Find frames where groups of fish (2 or more) are within the ROI
    if frame_bin == 1     
        if size(traj.grp2_XY_mod,2) == 0
            Inputs_provided.grp_string = 'Group1';
            [Distance, Quadrant_Stats, Time_Threshold] = get_frames_within_ROI(traj.grp1_XY_mod, traj.subject_XY_mod, Frames_per_sec,...
                FirstFrame, LastFrame,coordinates_y_grp,coordinates_X_sub,coordinates_Y_sub,Num_fish_close_to_subject,pixel_to_mm_change,Time_threshold);            
        elseif size(traj.grp1_XY_mod,2) == 0
            Inputs_provided.grp_string = 'Group2';
            [Distance, Quadrant_Stats, Time_Threshold] = get_frames_within_ROI(traj.grp2_XY_mod, traj.subject_XY_mod, Frames_per_sec,...
                FirstFrame, LastFrame,coordinates_y_grp,coordinates_X_sub,coordinates_Y_sub,Num_fish_close_to_subject,pixel_to_mm_change,Time_threshold);
        end
               
        % Save Files
        % Get a name file using all inputs to create unique files for each input
        name_file = ['D_sub_grp_input', '_T=',int2str(round(FirstFrame./Frames_per_sec)), 'to', int2str(round((LastFrame)./Frames_per_sec)),'secs',...
            '_ythresh_', int2str(Inputs_provided.Minimum_ythresh), '%to' , int2str(Inputs_provided.Maximum_ythresh), '%',...
            '_leastfish_', int2str(Inputs_provided.Num_fish_close_to_subject), '_timethreshold_', int2str(Inputs_provided.Time_threshold),'secs'];
        
        % 1. As a matfile
        save([Result_Folder_matfiles, filesep, name_file], 'Distance', 'Quadrant_Stats', 'Inputs_provided');
        
        % 2. As a figure
        %Save trajectories figure
        [fs1] = plot_trajectories(Distance,pixel_to_mm_change, Inputs_provided.grp_string);
        name_file = [SaveName,'Thresholded_Positions_Median', '_T=',int2str(round(FirstFrame./Frames_per_sec)), 'to', int2str(round((LastFrame)./Frames_per_sec)),'secs',...
            '_ythresh_', int2str(Inputs_provided.Minimum_ythresh), '%to' , int2str(Inputs_provided.Maximum_ythresh), '%', ...
            '_leastfish_', int2str(Inputs_provided.Num_fish_close_to_subject), '_timethreshold_', int2str(Inputs_provided.Time_threshold),'secs'];
        
        set(gcf, 'PaperPositionMode','auto','InvertHardCopy', 'off')
        saveas(fs1, [Result_Folder_figures, filesep, name_file], 'tif');
        
        %Save distance we had between subject and groups
        [fs2] = plot_distance_sub_grp(Distance, Inputs_provided.grp_string);
        name_file = [SaveName, 'Distance_sub_grp_', '_T=',int2str(round(FirstFrame./Frames_per_sec)), 'to', int2str(round((LastFrame)./Frames_per_sec)),'secs',...
            '_ythresh_', int2str(Inputs_provided.Minimum_ythresh), '%to' , int2str(Inputs_provided.Maximum_ythresh), '%', ...
            '_leastfish_', int2str(Inputs_provided.Num_fish_close_to_subject), '_timethreshold_', int2str(Inputs_provided.Time_threshold),'secs'];
        
        set(gcf, 'PaperPositionMode','auto','InvertHardCopy', 'off')
        saveas(fs2, [Result_Folder_figures, filesep, name_file], 'tif');
        
        % 3. As a excel file
        name_file = [SaveName, '_T=',int2str(round(FirstFrame./Frames_per_sec)), 'to', int2str(round((LastFrame)./Frames_per_sec)),'secs',...
            '_ythresh_', int2str(Inputs_provided.Minimum_ythresh), '%to' , int2str(Inputs_provided.Maximum_ythresh), '%', ...
            '_leastfish_', int2str(Inputs_provided.Num_fish_close_to_subject), '_timethreshold_', int2str(Inputs_provided.Time_threshold),'secs'];
        
        save_as_excel(Distance,Quadrant_Stats,Time_Threshold, name_file, Result_Folder_excel, Inputs_provided.grp_string);
        
    else
        
        
        for jj = FirstFrame:frame_bin:LastFrame %If the user requires binning of time
            
            if size(traj.grp2_XY_mod,2) == 0
                disp('Only Group 1 fish exist');
                Inputs_provided.grp_string = 'Group1';
                [Distance, Quadrant_Stats, Time_Threshold] = get_frames_within_ROI(traj.grp1_XY_mod, traj.subject_XY_mod, Frames_per_sec,...
                    jj, jj+frame_bin,coordinates_y_grp,coordinates_X_sub, coordinates_Y_sub, Num_fish_close_to_subject,pixel_to_mm_change, Time_threshold);
                
            elseif size(traj.grp1_XY_mod,2) == 0
                disp('Only Group 2 fish exist');
                Inputs_provided.grp_string = 'Group2';
                [Distance, Quadrant_Stats, Time_Threshold] = get_frames_within_ROI(traj.grp2_XY_mod, traj.subject_XY_mod, Frames_per_sec,...
                    jj, jj+frame_bin,coordinates_y_grp,coordinates_X_sub, coordinates_Y_sub, Num_fish_close_to_subject,pixel_to_mm_change, Time_threshold);
            end
            
            % Save files
            % 1. As a matfile
            %Get a name file using all inputs to create unique files for each input
            name_file = ['D_sub_grp_input', '_T=',int2str(round(FirstFrame./Frames_per_sec)), 'to', int2str(round((LastFrame)./Frames_per_sec)),'secs',...
                '_ythresh_', int2str(Inputs_provided.Minimum_ythresh), '%to' , int2str(Inputs_provided.Maximum_ythresh), '%', ...
                '_leastfish_', int2str(Inputs_provided.Num_fish_close_to_subject), '_timethreshold_', int2str(Inputs_provided.Time_threshold),'secs'];
            
            save([Result_Folder_matfiles, filesep, name_file], 'Distance', 'Inputs_provided');
            
            % 2. As a figure
            %Save trajectories figure
            [fs1] = plot_trajectories(Distance,pixel_to_mm_change,Inputs_provided.grp_string);
            name_file = ['Thresholded_Positions_Median', '_T=',int2str(round(FirstFrame./Frames_per_sec)), 'to', int2str(round((LastFrame)./Frames_per_sec)),'secs',...
                '_ythresh_', int2str(Inputs_provided.Minimum_ythresh), '%to' , int2str(Inputs_provided.Maximum_ythresh), '%', ...
                '_leastfish_', int2str(Inputs_provided.Num_fish_close_to_subject), '_timethreshold_', int2str(Inputs_provided.Time_threshold),'secs'];
            
            set(gcf, 'PaperPositionMode','auto','InvertHardCopy', 'off')
            saveas(fs1, [Result_Folder_figures, filesep, name_file], 'tif');
            
            
            %Save distance we had between subject and groups
            [fs2] = plot_distance_sub_grp(Distance,Inputs_provided.grp_strings);
            name_file = ['Distance_sub_grp_', '_T=',int2str(round(FirstFrame./Frames_per_sec)), 'to', int2str(round((LastFrame)./Frames_per_sec)),'secs',...
                '_ythresh_', int2str(Inputs_provided.Minimum_ythresh), '%to' , int2str(Inputs_provided.Maximum_ythresh), '%', ...
                '_leastfish_', int2str(Inputs_provided.Num_fish_close_to_subject), '_timethreshold_', int2str(Inputs_provided.Time_threshold),'secs'];
            
            set(gcf, 'PaperPositionMode','auto','InvertHardCopy', 'off')
            saveas(fs2, [Result_Folder_figures, filesep, name_file], 'tif');
            
            % 3. As a excel file
            name_file = [SaveName, '_T=',int2str(round(FirstFrame./Frames_per_sec)), 'to', int2str(round((LastFrame)./Frames_per_sec)),'secs',...
                '_ythresh_', int2str(Inputs_provided.Minimum_ythresh), '%to' , int2str(Inputs_provided.Maximum_ythresh), '%', ...
                '_leastfish_', int2str(Inputs_provided.Num_fish_close_to_subject), '_timethreshold_', int2str(Inputs_provided.Time_threshold),'secs'];
            
            save_as_excel(Distance,Quadrant_Stats,Time_Threshold, name_file, Result_Folder_excel, Inputs_provided.grp_string);
            
        end
    end
end
end

%% Function to find frames where groups of fish are within user given ROI
function [Distance, Quadrant_Stats, Time_Threshold] = get_frames_within_ROI(grp_XY_mod, subject_XY_mod, Frames_per_sec,...
        FirstFrame, LastFrame,coordinates_y_grp, coordinates_X_sub, coordinates_Y_sub, Num_fish,pixel_to_mm_change, Time_threshold)
    

%% Group
% 1. Find frames and x,y position where groups of fish were within specified ROI
Data_grpy = grp_XY_mod(FirstFrame:LastFrame,:,2);
Data_grp = grp_XY_mod(FirstFrame:LastFrame,:,:);
Frames_grp = find((sum(Data_grpy>coordinates_y_grp(1)&Data_grpy<coordinates_y_grp(2),2)>=Num_fish)==1);
Frames_grp_pos = Data_grp(Frames_grp,:,:);

% 2. Get distances from corresponding frames in the subject
Data_subj = subject_XY_mod(FirstFrame:LastFrame,:,:);
Frames_subject_grp_pos = Data_subj(Frames_grp,:,:);

% 3. Get center of mass by taking median of x,y positions of group
Centre_mass_grp = squeeze(median(Frames_grp_pos,2));

% 4.Take euclidean distances between centre of mass of group and subject
Dist_subject_from_grp_pos = sqrt((squeeze(Frames_subject_grp_pos(:,1,1))-Centre_mass_grp(:,1)).^2 + ...
    (squeeze(Frames_subject_grp_pos(:,1,2))-Centre_mass_grp(:,2)).^2 );

% 5. Convert pixels to mm
Dist_subject_from_grp_pos_mm = Dist_subject_from_grp_pos/pixel_to_mm_change;
Centre_mass_grp_mm = Centre_mass_grp/pixel_to_mm_change;

% 6. Find quadrants where the subjects are located. There are 4 quadrants.
Quadrant_subject_grp_pos = zeros(size(Frames_subject_grp_pos,1),1);

Quadrant_subject_grp_pos(Frames_subject_grp_pos(:,1)<coordinates_X_sub & Frames_subject_grp_pos(:,2)>coordinates_Y_sub) = 1; %Group1 close
Quadrant_subject_grp_pos(Frames_subject_grp_pos(:,1)>coordinates_X_sub & Frames_subject_grp_pos(:,2)>coordinates_Y_sub) = 2; %Group2 close

Quadrant_subject_grp_pos(Frames_subject_grp_pos(:,1)<coordinates_X_sub & Frames_subject_grp_pos(:,2)<coordinates_Y_sub) = 3; %Group1 far
Quadrant_subject_grp_pos(Frames_subject_grp_pos(:,1)>coordinates_X_sub & Frames_subject_grp_pos(:,2)<coordinates_Y_sub) = 4; %Group2 far

% 7. Find frames where groups of fish were near subject for the amount
% of time specified by user (in seconds)
Frames_sec_grp = Frames_grp./Frames_per_sec;
count = 1;
temp_seconds = unique(fix(Frames_sec_grp));

Time_found_grp = 0;
Time_centremass_grp = 0;
Time_subdistance_grp = 0;
Time_quadrant_grp = 0;

for ii = 1:length(temp_seconds)
    if size(find(ismember(fix(Frames_sec_grp), temp_seconds(ii):temp_seconds(ii)+Time_threshold-1)==1),1)>=Frames_per_sec*Time_threshold
        Time_found_grp(count,1) = temp_seconds(ii);
        Time_centremass_grp(count,1) = mean(Centre_mass_grp_mm(ismember(fix(Frames_sec_grp), temp_seconds(ii):temp_seconds(ii)+Time_threshold-1)==1));
        Time_subdistance_grp(count,1) = mean(Dist_subject_from_grp_pos_mm(ismember(fix(Frames_sec_grp), temp_seconds(ii):temp_seconds(ii)+Time_threshold-1)==1));
        Time_quadrant_grp(count,1) = fix(median(Quadrant_subject_grp_pos(ismember(fix(Frames_sec_grp), temp_seconds(ii):temp_seconds(ii)+Time_threshold-1)==1)));
        count = count+1;
    end
end


% 9. Seperate stats based on Quadrant for all frames and those thresholded
% based on time

for ii = 1:4
    
    Quadrant_num(ii,1) = ii;
    Quadrant_count_grp(ii,1) = size(find(Quadrant_subject_grp_pos==ii),1);
    Quadrant_centremass_grp(ii,1) = mean(Centre_mass_grp_mm(Quadrant_subject_grp_pos==ii));
    Quadrant_subdist_grp(ii,1) = mean(Dist_subject_from_grp_pos_mm(Quadrant_subject_grp_pos==ii));  
    
    Time_quadrant_count_grp(ii,1) = size(find(Time_quadrant_grp==ii),1);
    Time_quadrant_centremass_grp(ii,1) = mean(Time_centremass_grp(Time_quadrant_grp==ii));
    Time_quadrant_subdist_grp(ii,1) = mean(Time_subdistance_grp(Time_quadrant_grp==ii));
end


%% Save as Matfile
Distance.Frames_grp = Frames_grp;
Distance.Frames_grp_sec = Frames_sec_grp;
Distance.Frames_grp_pos = Frames_grp_pos;
Distance.Frames_subject_grp_pos = Frames_subject_grp_pos;
Distance.Centre_mass_grp = Centre_mass_grp;
Distance.Centre_mass_grp_mm = Centre_mass_grp_mm;
Distance.Dist_subject_from_grp_pos = Dist_subject_from_grp_pos;
Distance.Dist_subject_from_grp_pos_mm = Dist_subject_from_grp_pos_mm;
Distance.Quadrant_subject_grp_pos = Quadrant_subject_grp_pos;

Quadrant_Stats.Quadrant_num = Quadrant_num;
Quadrant_Stats.Quadrant_count_grp = Quadrant_count_grp;
Quadrant_Stats.Quadrant_centremass_grp = Quadrant_centremass_grp;
Quadrant_Stats.Quadrant_subdist_grp = Quadrant_subdist_grp;

Time_Threshold.Time_found_grp = Time_found_grp;
Time_Threshold.Time_centremass_grp = Time_centremass_grp;
Time_Threshold.Time_subdistance_grp = Time_subdistance_grp;
Time_Threshold.Time_quadrant_grp = Time_quadrant_grp;

Time_Threshold.Time_quadrant_count_grp = Time_quadrant_count_grp;
Time_Threshold.Time_quadrant_centremass_grp =  Time_quadrant_centremass_grp;
Time_Threshold.Time_quadrant_subdist_grp = Time_quadrant_subdist_grp;

end


%% Save as Excel
function save_as_excel(Distance, Quadrant_Stats, Time_Threshold, name_file, Result_Folder_excel, grp_string)

% 1. First save distances with names suitable for excel files

Distance1.([grp_string,'_frames']) = Distance.Frames_grp;
Distance1.([grp_string,'_frames_sec']) = fix(Distance.Frames_grp_sec);
Distance1.([grp_string,'_centre_mass_mm']) = Distance.Centre_mass_grp_mm;
Distance1.([grp_string,'_subjectdist_mm']) = Distance.Dist_subject_from_grp_pos_mm;
Distance1.([grp_string,'_subject_quadrant']) = Distance.Quadrant_subject_grp_pos;

group_matrices = 5; % to know how many Group2 vs Group1 quantifications - provide a gap while saving in excel


Distance1.Quadrant = Quadrant_Stats.Quadrant_num;
Distance1.(['Count_',grp_string]) = Quadrant_Stats.Quadrant_count_grp;
Distance1.(['Centremass_',grp_string]) = Quadrant_Stats.Quadrant_centremass_grp;
Distance1.(['Subdist_',grp_string]) = Quadrant_Stats.Quadrant_subdist_grp;

quad_matrices = 9;

Distance1.([grp_string,'_time_secs']) = Time_Threshold.Time_found_grp;
Distance1.([grp_string,'_time_centre_mass_mm']) = Time_Threshold.Time_centremass_grp;
Distance1.([grp_string,'_time_subjectdist_mm']) = Time_Threshold.Time_subdistance_grp;
Distance1.([grp_string,'_time_subject_quadrant']) = Time_Threshold.Time_quadrant_grp;

group_time_matrices = 13;

Distance1.Time_Quadrant = Quadrant_Stats.Quadrant_num;
Distance1.(['Time_Count_',grp_string]) = Time_Threshold.Time_quadrant_count_grp;
Distance1.(['Time_Centremass_',grp_string]) = Time_Threshold.Time_quadrant_centremass_grp;
Distance1.(['Time_Subdist_',grp_string]) = Time_Threshold.Time_quadrant_subdist_grp;


clear Distance Quadrant_Stats

Temp_Dat = fieldnames(Distance1);
filename = [ Result_Folder_excel,filesep,name_file,'.xls'];
fid = fopen(filename, 'w+');

% Go through and save as a cell in a format suitable for excel files
count = 0;
for kk = 1:length(Temp_Dat)
    if kk == group_matrices+1 || kk == quad_matrices+1 || kk == group_time_matrices+1
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
function [fs1] = plot_trajectories(Distance,pixel_to_mm_change, grp_string)

% Plot trajectories in chosen frames for user and save as
% Thresholded_Positions Distance_sub_grp and in Figures folder

fs1 = figure(1);
clf;
set(fs1,'color','white')
hold on
plot(squeeze(Distance.Frames_grp_pos(:,:,1)), squeeze(Distance.Frames_grp_pos(:,:,2)),'.')
plot(Distance.Centre_mass_grp(:,1), Distance.Centre_mass_grp(:,2), 'k+', 'MarkerSize', 10)
plot(squeeze(Distance.Frames_subject_grp_pos(:,:,1)), squeeze(Distance.Frames_subject_grp_pos(:,:,2)),'r.', 'MarkerSize', 10)
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
for ii = 1:size(Distance.Frames_grp_pos, 2)
    legend_str1{ii} = [grp_string, ' Fish',int2str(ii)];
end
legend_str1{size(Distance.Frames_grp_pos, 2)+1} = ['Center of Mass ', grp_string];

legend_str1{size(Distance.Frames_grp_pos, 2)+2} = ['Subject and ', grp_string];

legend(legend_str1{:});

end



function [fs2] = plot_distance_sub_grp(Distance, grp_string)
% 2. Plot distance from centre of mass of groups and subject for selected
% frames

fs2 = figure(2);
set(fs2,'color','white')
subplot(2,1,1)
plot(Distance.Frames_grp, Distance.Dist_subject_from_grp_pos_mm, 'r.')

set(gca, 'TickDir','out', 'FontSize',12)
box off
xlabel(gca,'Frame Number', 'FontSize',12);
ylabel(gca,'distance (mm)', 'FontSize',12);
title(gca, 'Distance of subject from groups of fish', 'FontSize',12)

legend(['Subject and',grp_string]);

subplot(2,1,2)
[mean1, ci1] = normfit(Distance.Dist_subject_from_grp_pos_mm);
errorbar(1,mean1,ci1,'r')

set(gca, 'TickDir','out', 'FontSize',12)
box off
xlabel(gca,'Frame', 'FontSize',12);
ylabel(gca,'distance (mm)', 'FontSize',12);
title(gca, 'Distance of subject from groups of fish', 'FontSize',12)

legend(['Subject and ',grp_string]);

end


