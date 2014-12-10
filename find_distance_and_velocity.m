function find_distance_and_velocity(Tbin)

%% Find distance and velocity from trajectories. Calculate various parameters based on velocity. Save as excel. All trajectory files in a directory are taken
%% Input -
% Tbin : Time bins for which to calculate parameters. Output will be over
% each time bin and also over all timepoints taken together
%% Output -
% Excel file has outputs Above(s),Above(%),Bottom(s),Bottom(%),Pausing(s),Freezing(s),Freezing Eps,Darting Eps


%% Some default parameters
Frames_per_sec = 30; %Frames per second
Bin_size = 30; %Bin size
Pixels_to_mm = 0.94; %Pixels to mm
Freezing_in_bins = 5; %Sliding window size for freezing in bins
Pause_speed_threshold = 3.5; %If the fish goes below this speed, it will be considered pausing. 
    
%% Main Script

%Get directory from which to take all trajectories.mat file
PathName = uigetdir(pwd);
FolderNames = dir(PathName);

%Get sub folder in the directory
isub = [FolderNames(:).isdir];
FolderNames = {FolderNames(isub).name};
FolderNames(ismember(FolderNames,{'.','..'})) = [];

%Loop through subfolders
for ii = 1:length(FolderNames)
    TrajFile = dir([PathName,filesep,FolderNames{ii}]);
    
    %Only run analysis on those that contain modified_trajectories.mat
    for jj = 1:length(TrajFile)
        if ~isempty(strfind(TrajFile(jj).name, '_modified_trajectories.mat'))
            
            %Load trajectories.mat file and divide into groups and subjects
            traj = load([PathName,filesep,FolderNames{ii},filesep,TrajFile(jj).name]);
            %Concatenate all for easier processing
            trajectories_all = [traj.grp1_XY_mod, traj.grp2_XY_mod, traj.subject_XY_mod];
            
            %Get distance and velocity between each frame for each animal
            num_bins = mod(size(trajectories_all,1)-1,Bin_size);
            speed_pixel_by_bin = zeros(fix((size(trajectories_all,1)-1)/Bin_size)+(1*num_bins~=0),size(trajectories_all,2));
            
            for kk = 1:size(trajectories_all,2)
                dist1 = squeeze(trajectories_all(1:end-1,kk,:));
                dist2 = squeeze(trajectories_all(2:end,kk,:));
                euc_dist_perframe(:,kk) = sqrt((dist1(:,1)-dist2(:,1)).^2 + (dist1(:,2)-dist2(:,2)).^2); %euclidean distance
                
                %Find speed in bins - divide by number of bins to round up.
                speed_pixel_by_bin(1:end-(1*num_bins~=0),kk) = sum(reshape(euc_dist_perframe(1:end-num_bins,kk), Bin_size, fix(size(euc_dist_perframe,1)/Bin_size)));
                if num_bins~=0 %Sum remaining time points
                    speed_pixel_by_bin(end,kk) = sum(euc_dist_perframe(end-num_bins+1:end,kk));
                end
                speed_mm_by_bin(:,kk) = speed_pixel_by_bin(:,kk)/Pixels_to_mm; %Convert to mm
                speed_mm_by_sec(:,kk) = speed_mm_by_bin(:,kk)/(Bin_size/Frames_per_sec);                
            end
            A
            
            
        end
    end
end
