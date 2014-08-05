function idTracker_Heatmap(tbin, fps)
%% Get Directory from user - directory can contain one or more modified trajectories
%% Input :
%tbin - bin data over specified time bin and plot heatmaps for each
%fps - frames per second

%Radius
Rad = 8;

%Colorbar scaling
cmin = 0;
cmax = 2;

PathName = uigetdir(pwd, 'Select modified trajectories file');
FileName = dir([PathName, filesep,'*modified*.mat']);



for ii = 1:length(FileName)
    %Save Figures in This Folder
    SaveName = FileName(ii).name(1:strfind(FileName(1).name, 'modified')-2);
    Result_Folder = [PathName, filesep, 'Figures', filesep, SaveName];
    mkdir(Result_Folder);
    
    traj = load([PathName, filesep, FileName(ii).name]);
    NumFrames = size(traj.grp1_XY_mod, 1);
    frame_bin = tbin*fps;
    
    %Take number of frames divisible by the time bin and collapse the rest
    LastFrame = (fix(NumFrames/frame_bin)-1)*frame_bin;
    
    %Create a timespent variable spanning all x,y
    All_traj_X = [reshape(traj.grp1_XY_mod(:,:,1),size(traj.grp1_XY_mod,1)*size(traj.grp1_XY_mod,2),1); ...
        reshape(traj.grp2_XY_mod(:,:,1),size(traj.grp2_XY_mod,1)*size(traj.grp2_XY_mod,2),1);...
        traj.subject_XY_mod(:,1,1)];
    All_traj_Y = [reshape(traj.grp1_XY_mod(:,:,2),size(traj.grp1_XY_mod,1)*size(traj.grp1_XY_mod,2),1); ...
        reshape(traj.grp2_XY_mod(:,:,2),size(traj.grp2_XY_mod,1)*size(traj.grp2_XY_mod,2),1);...
        traj.subject_XY_mod(:,1,2)];
    
    min_traj_X = min(All_traj_X);
    max_traj_X = max(All_traj_X);
    min_traj_Y = min(All_traj_Y);
    max_traj_Y = max(All_traj_Y);
    
    clear All_traj_X All_traj_Y
    
    OverallTimeSpent = zeros(round(max_traj_X), round(max_traj_Y));
    
    
    %Go through each time bin and plot average heatmap
    for jj = 1:frame_bin:LastFrame
        
        TimeSpent = zeros(round(max_traj_X), round(max_traj_Y));
        
        %Group1
        
        [TimeSpent] = find_mean_time(traj.grp1_XY_mod(jj:jj+frame_bin-1,:,1), traj.grp1_XY_mod(jj:jj+frame_bin-1,:,2),TimeSpent,1);
        
        %Group2
        [TimeSpent] = find_mean_time(traj.grp2_XY_mod(jj:jj+frame_bin-1,:,1), traj.grp2_XY_mod(jj:jj+frame_bin-1,:,2),TimeSpent,2);
        
        %Subject
        [TimeSpent] = find_mean_time(traj.subject_XY_mod(jj:jj+frame_bin-1,:,1), traj.subject_XY_mod(jj:jj+frame_bin-1,:,2),TimeSpent,3);
        
        plot_heatmap(TimeSpent,Rad, cmin, cmax,min_traj_X, min_traj_Y, max_traj_X, max_traj_Y, fps, Result_Folder,jj, frame_bin)
        
        OverallTimeSpent = OverallTimeSpent + TimeSpent;
    end
    
    plot_heatmap(OverallTimeSpent,3, 0, 2, min_traj_X, min_traj_Y, max_traj_X, max_traj_Y, fps, Result_Folder)   
end

end


function plot_heatmap(TimeSpent,Rad,cmin,cmax, min_traj_X, min_traj_Y, max_traj_X, max_traj_Y, fps, Result_Folder, jj, frame_bin)

%Plot Heatmaps
S=+(bwdist(padarray(1,[1,1]*double(round(Rad*1.5))))<=Rad);
Filt_TimeSpent=double(convn(TimeSpent,S,'same'));
Filt_TimeSpent = smoothn(Filt_TimeSpent,5);

fs1 = figure(1);
set(fs1,'color','white')
pcolor((Filt_TimeSpent(min_traj_X:max_traj_X,min_traj_Y:max_traj_Y)./fps)')
caxis([cmin cmax]) % Colorbar setting
colorbar
colormap(jet(2000))
shading interp
set(gca, 'TickDir','out', 'FontSize',12)
box off
xlabel('x distance (mm)', 'FontSize',12);
ylabel('y distance (mm)', 'FontSize',12);

if nargin == 12
    name_file = ['HeatMap_T= ',int2str(round(jj./fps)), ' To ', int2str(round((jj+frame_bin-1)./fps)),'secs'];
else
    name_file = ['HeatMap_AllTimePoints'];
end

set(gcf, 'PaperPositionMode','auto','InvertHardCopy', 'off')
saveas(fs1, [Result_Folder, filesep, name_file], 'jpg');
end

