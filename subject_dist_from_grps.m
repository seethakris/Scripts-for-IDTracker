function subject_dist_from_grps(fps, TMin, TMax, tbin, y_thresh_min, y_thresh_max,Num_fish_close_to_subject)

%% Find average distance of subject when fish from the groups are closer to the subject tank
%% TMin and TMax are optional inputs, if left empty, the entire time is considered

close all
warning off

%% Input :
%% Optional Inputs :
% You can use [] to use default: Default values can be changed below as
% indicated in the script. Using approx. 3.05 pixels/mm
%
%   fps - frames per second. Default - 30
%   TMin - Minimum time (seconds). Default - First Frame in video
%   TMax - Maximum time (seconds). Default - Last Frame in video
%   tbin - bin data over specified time bin (seconds). Default - one frame
%   y_thresh_min and y_thresh_max - ROI over which to analyze group
%   behavior (in percent). Default - 0% (border closest to subject fish) to 25%
%   Num_fish_close_to_subject - the least number of fish that need to be near the subject
%   to include that frame. Default - 3

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


Inputs_provided.Frames_per_sec = Frames_per_sec;
Inputs_provided.frame_bin = frame_bin;
Inputs_provided.FirstFrame = FirstFrame;
Inputs_provided.Minimum_ythresh = Minimum_ythresh;
Inputs_provided.Maximum_ythresh = Maximum_ythresh;
Inputs_provided.Num_fish_close_to_subject = Num_fish_close_to_subject;


%% Main Script
PathName = uigetdir(pwd, 'Select modified trajectories file');
FileName = dir([PathName, filesep,'*modified*.mat']);


% Find y threshold using data of all fish to find the most accurate boundary.
grp_traj_X = [];
grp_traj_Y = [];
subject_traj_X = [];
subject_traj_Y = [];


for ii = 1:length(FileName)
    % Load trajectories
    traj = load([PathName, filesep, FileName(ii).name]);
    
    % Collect all group X and Y coordinates
    grp_traj_X = [grp_traj_X; reshape(traj.grp1_XY_mod(:,:,1),size(traj.grp1_XY_mod,1)*size(traj.grp1_XY_mod,2),1); ...
        reshape(traj.grp2_XY_mod(:,:,1),size(traj.grp2_XY_mod,1)*size(traj.grp2_XY_mod,2),1)];
    grp_traj_Y = [grp_traj_Y; reshape(traj.grp1_XY_mod(:,:,2),size(traj.grp1_XY_mod,1)*size(traj.grp1_XY_mod,2),1); ...
        reshape(traj.grp2_XY_mod(:,:,2),size(traj.grp2_XY_mod,1)*size(traj.grp2_XY_mod,2),1)];
end

% Find y min and y max using the grp trajectories
coordinates_y_grp(1) = [min(grp_traj_Y) + (max(grp_traj_Y)-min(grp_traj_Y))*(Minimum_ythresh/100)];
coordinates_y_grp(2) = [min(grp_traj_Y) + (max(grp_traj_Y)-min(grp_traj_Y))*(Maximum_ythresh/100)];

clear grp_traj_X grp_traj_Y traj

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
        get_frames_within_ROI(traj, FirstFrame, LastFrame,Frames_per_sec,coordinates_y_grp,Num_fish_close_to_subject, ...
            Result_Folder_excel,Result_Folder_figures,Result_Folder_matfiles,Inputs_provided);
        
    else
        for jj = FirstFrame:frame_bin:LastFrame %If the user requires binning of time
            get_frames_within_ROI(traj, jj, jj+frame_bin,Frames_per_sec,coordinates_y_grp,Num_fish_close_to_subject,...
                Result_Folder_excel,Result_Folder_figures,Result_Folder_matfiles, Inputs_provided);
        end
    end
    
end
end

%% Function to find frames where groups of fish are within user given ROI
function get_frames_within_ROI(traj, FirstFrame, LastFrame,Frames_per_sec,coordinates_y_grp,Num_fish,...
    Result_Folder_excel,Result_Folder_figures,Result_Folder_matfiles,Inputs_provided)

close all

%Group 1
Data_grpy = traj.grp1_XY_mod(FirstFrame:LastFrame,:,2);
Data_grp1 = traj.grp1_XY_mod(FirstFrame:LastFrame,:,:);
Frames_grp1 = find((sum(Data_grpy>coordinates_y_grp(1)&Data_grpy<coordinates_y_grp(2),2)>=Num_fish)==1); %Find frames where atleast two fish are within ROI
Frames_grp1_pos = Data_grp1(Frames_grp1,:,:);

Data_subj = traj.subject_XY_mod(FirstFrame:LastFrame,:,:);
Frames_subject_grp1_pos = Data_subj(Frames_grp1,:,:); % Get distances from corresponding frames in the subject

Centre_mass_grp1 = squeeze(median(Frames_grp1_pos,2)); %Centre of mass of group1
Dist_subject_from_grp1_pos = sqrt((squeeze(Frames_subject_grp1_pos(:,1,1))-Centre_mass_grp1(:,1)).^2 + ...
    (squeeze(Frames_subject_grp1_pos(:,1,2))-Centre_mass_grp1(:,2)).^2 ); %Take equilidean distance of centre of mass of group and subject
Dist_subject_from_grp1_pos_mm = Dist_subject_from_grp1_pos/3.05;


% Group2
Data_grpy = traj.grp2_XY_mod(FirstFrame:LastFrame,:,2);
Data_grp2 = traj.grp2_XY_mod(FirstFrame:LastFrame,:,:);
Frames_grp2 = find((sum(Data_grpy>coordinates_y_grp(1)&Data_grpy<coordinates_y_grp(2),2)>=Num_fish)==1); %Find frames where atleast two fish are within ROI
Frames_grp2_pos = Data_grp2(Frames_grp2,:,:);

Frames_subject_grp2_pos = Data_subj(Frames_grp2,:,:); % Get distances from corresponding frames in the subject

Centre_mass_grp2 = squeeze(median(Frames_grp2_pos,2)); %Centre of mass of group2
Dist_subject_from_grp2_pos = sqrt((squeeze(Frames_subject_grp2_pos(:,1,1))-Centre_mass_grp2(:,1)).^2 + ...
    (squeeze(Frames_subject_grp2_pos(:,1,2))-Centre_mass_grp2(:,2)).^2 ); %Take equilidean distance of centre of mass of group and subject
Dist_subject_from_grp2_pos_mm = Dist_subject_from_grp2_pos/3.05;


%% Save as Matfile
Distance.Frames_grp1 = Frames_grp1;
Distance.Frames_grp1_pos = Frames_grp1_pos;
Distance.Frames_subject_grp1_pos = Frames_subject_grp1_pos;
Distance.Centre_mass_grp1 = Centre_mass_grp1;
Distance.Dist_subject_from_grp1_pos = Dist_subject_from_grp1_pos;
Distance.Dist_subject_from_grp1_pos_mm = Dist_subject_from_grp1_pos_mm;

Distance.Frames_grp2 = Frames_grp2;
Distance.Frames_grp2_pos = Frames_grp2_pos;
Distance.Frames_subject_grp2_pos = Frames_subject_grp2_pos;
Distance.Centre_mass_grp2 = Centre_mass_grp2;
Distance.Dist_subject_from_grp2_pos = Dist_subject_from_grp2_pos;
Distance.Dist_subject_from_grp2_pos_mm = Dist_subject_from_grp2_pos_mm;

%Get a name file using all inputs to create unique files for each input
name_file = ['D_sub_grp_input', '_T=',int2str(round(FirstFrame./Frames_per_sec)), 'to', int2str(round((LastFrame)./Frames_per_sec)),'secs',...
    '_ythresh_', int2str(Inputs_provided.Minimum_ythresh), '%to' , int2str(Inputs_provided.Maximum_ythresh), '%', '_leastfish_', int2str(Inputs_provided.Num_fish_close_to_subject)];

save([Result_Folder_matfiles, filesep, name_file], 'Distance', 'Inputs_provided');

clear Distance
%% Save as Excel



%% Plotting:
% 1. Plot trajectories in chosen frames for user and save as
% Plot_Thresholded_Positions in Figures folder
fs1 = figure(1);
set(fs1,'color','white')
hold on
plot(squeeze(Frames_grp1_pos(:,:,1)), squeeze(Frames_grp1_pos(:,:,2)),'.')
plot(Centre_mass_grp1(:,1), Centre_mass_grp1(:,2), 'k+', 'MarkerSize', 10)
plot(squeeze(Frames_grp2_pos(:,:,1)), squeeze(Frames_grp2_pos(:,:,2)),'*')
plot(Centre_mass_grp2(:,1), Centre_mass_grp2(:,2), 'k+', 'MarkerSize', 10)
plot(squeeze(Frames_subject_grp1_pos(:,:,1)), squeeze(Frames_subject_grp1_pos(:,:,2)),'r.', 'MarkerSize', 10)
plot(squeeze(Frames_subject_grp2_pos(:,:,1)), squeeze(Frames_subject_grp2_pos(:,:,2)),'g*')
hold off
set(gca, 'TickDir','out', 'FontSize',12)
set(gca, 'YDir', 'reverse');
box off
%Convert to mm
x = get(gca, 'Xtick');
set(gca, 'xTickLabel',strread(int2str(x/3.05),'%s'));
y = get(gca, 'Ytick');
set(gca, 'yTickLabel',strread(int2str(y/3.05),'%s'));
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

%Save figure
name_file = ['Thresholded_Positions_Median', '_T=',int2str(round(FirstFrame./Frames_per_sec)), 'to', int2str(round((LastFrame)./Frames_per_sec)),'secs',...
    '_ythresh_', int2str(Inputs_provided.Minimum_ythresh), '%to' , int2str(Inputs_provided.Maximum_ythresh), '%', ...
    '_leastfish_', int2str(Inputs_provided.Num_fish_close_to_subject)];

set(gcf, 'PaperPositionMode','auto','InvertHardCopy', 'off')
saveas(fs1, [Result_Folder_figures, filesep, name_file], 'tif');

% 2. Plot distance from centre of mass of groups and subject for selected
% frames
fs1 = figure(2);
set(fs1,'color','white')
subplot(2,1,1)
plot(Frames_grp1, Dist_subject_from_grp1_pos_mm, 'r.')
hold on
plot(Frames_grp2, Dist_subject_from_grp2_pos_mm, 'g.')
hold off

set(gca, 'TickDir','out', 'FontSize',12)
box off
xlabel(gca,'Frame Number', 'FontSize',12);
ylabel(gca,'distance (mm)', 'FontSize',12);
title(gca, 'Distance of subject from groups of fish', 'FontSize',12)

legend('Subject and Group1', 'Subject and Group2');

subplot(2,1,2)
[mean1, ci1] = normfit(Dist_subject_from_grp1_pos_mm);
errorbar(1,mean1,ci1,'r')
hold on
[mean2, ci2] = normfit(Dist_subject_from_grp2_pos_mm);
errorbar(2,mean2,ci2,'g')
hold off

set(gca, 'TickDir','out', 'FontSize',12)
box off
xlabel(gca,'Frame', 'FontSize',12);
ylabel(gca,'distance (mm)', 'FontSize',12);
title(gca, 'Distance of subject from groups of fish', 'FontSize',12)

legend('Subject and Group1', 'Subject and Group2');

%Save figure
name_file = ['Distance_sub_grp_', '_T=',int2str(round(FirstFrame./Frames_per_sec)), 'to', int2str(round((LastFrame)./Frames_per_sec)),'secs',...
    '_ythresh_', int2str(Inputs_provided.Minimum_ythresh), '%to' , int2str(Inputs_provided.Maximum_ythresh), '%', ...
    '_leastfish_', int2str(Inputs_provided.Num_fish_close_to_subject)];

set(gcf, 'PaperPositionMode','auto','InvertHardCopy', 'off')
saveas(fs1, [Result_Folder_figures, filesep, name_file], 'tif');

end
