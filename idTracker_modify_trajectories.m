function idTracker_modify_trajectories
%% Load data, enter fish number for each group as seen during tracking and correct trajectories.
%% Select trajectories.mat file for the experiment, enter Fish Numbers (only for first time)
%  Group1 are fish on left hand side , Group 2 are on the right

warning off
close all

%User Input
threshold = 500;

%Ask user to load trajectories.mat file

[FileName, PathName] = uigetfile('*.mat', 'Select trajectories file');
Traj = load([PathName, FileName]);

numFish = size(Traj.trajectories,2);
numFrames = size(Traj.trajectories,1);
All_fish = 1:numFish;


%% Check if trajectories.mat ahs a field called fish number, else ask theuser and save it

if exist([PathName,FileName(1:end-4),'_FishNumber.mat'])
    load([PathName,FileName(1:end-4),'_FishNumber.mat'],'FishNumber');
    grp1_fish = FishNumber.grp1;
    grp2_fish = FishNumber.grp2;
    subject_fish = FishNumber.subject;
    disp(sprintf(['Grp1 Fish: ', num2str(grp1_fish), '\nGrp2 Fish: ', num2str(grp2_fish),'\nSubject fish: ', num2str(subject_fish)]));
    
else
    % Open dialogbox for user to enter grp 1, grp 2 and subject
    
    prompt = {'Enter fish#s for grp 1:','Enter fish# for subject:'};
    dlg_title = 'Input fish numbers as tracked';
    num_lines = 1;
    answer = inputdlg(prompt,dlg_title,num_lines);
    
    if size(answer,1)~=2
        disp('Cancelled Operation');
        return;
        
    else
        
        grp1_fish = str2num(answer{1});
        subject_fish = str2num(answer{2});
        grp2_fish = setdiff(setdiff(All_fish,grp1_fish),subject_fish);
        
        %Confirm with user
        answer{3} = sprintf('%.0f,' , grp2_fish);
        answer{3} = answer{3}(1:end-1);
        
        button = questdlg(sprintf(['Grp1 Fish: ', answer{1}, '\nGrp2 Fish: ', answer{3},'\nSubject fish: ', answer{2}]),...
            'You have chosen:', 'OK', 'Cancel','OK');
        
        disp(sprintf(['Grp1 Fish: ', answer{1}, '\nGrp2 Fish: ', answer{3},'\nSubject fish: ', answer{2}]));
        
        
        if strcmp(button,'Cancel')
            disp('Canceled file operation')
            return;
            
        elseif strcmp(button,'OK')
            FishNumber.grp1 = grp1_fish;
            FishNumber.grp2 = grp2_fish;
            FishNumber.subject = subject_fish;
            save([PathName,FileName(1:end-4),'_FishNumber.mat'],'FishNumber');
        end
    end
end

%% Analysis
%Start correcting errors and making heatmaps

grp1_XY = Traj.trajectories(:,grp1_fish,:);
grp2_XY = Traj.trajectories(:,grp2_fish,:);
subject_XY = Traj.trajectories(:,subject_fish,:);

%Approximate area that should be occupied by groups of fish.
%[xleft,xright,ytop,ybottom];
Group1_area =[100, 545,  300, 144];
Group2_area =[310, 580, 280, 50];
Subject_area =[100, 545,  300, 144];


% Remove any misidentification in trajectory by finding big jumps
% between consecutive points and interpolate -
% seperately for each group of fish

if size(grp1_XY,2)==0
    legend_string1 = {'grp2', 'subject'};
    legend_string2 = {'grp1','mean grp1','subject','mean subject'};
elseif size(grp2_XY,2)==0
    legend_string1 = {'grp1', 'subject'};
    legend_string2 = {'grp2','mean grp2', 'subject','mean subject'};
elseif size(subject_XY,2)==0
    legend_string1 = {'grp1', 'grp2'};
    legend_string2 = {'grp1','mean grp1', 'grp2','mean grp2'};
end

[grp1_XY_mod] = fix_trajectories(grp1_XY,threshold,1,Group1_area, legend_string1, legend_string2);
[grp2_XY_mod] = fix_trajectories(grp2_XY,threshold,2,Group2_area, legend_string1, legend_string2);
[subject_XY_mod] = fix_trajectories(subject_XY,threshold,3, Subject_area, legend_string1, legend_string2);

%Save Figures in This Folder
Result_Folder = [PathName, filesep, 'Figures', filesep, FileName(1:end-4), filesep];
mkdir(Result_Folder);

fs = figure(1);
print('-djpeg', [Result_Folder,  'Trajectories Before correction.jpeg']);
fs = figure(2);
print('-djpeg', [Result_Folder, 'Trajectories After correction.jpeg']);

save([PathName,FileName(1:end-4), '_modified_trajectories.mat'], 'grp1_XY_mod', 'grp2_XY_mod', 'subject_XY_mod');
