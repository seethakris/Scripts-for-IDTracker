function [modified_trajectories] = fix_trajectories(traj,threshold,flag, area, legend_string1, legend_string2)


numFish = size(traj,2);
numFrames = size(traj,1);

color = ['r','g','b'];

fs = figure(1);
set(fs, 'color','white');
hold on
temp = reshape(traj, size(traj,1)*size(traj,2),2);
plot(temp(:,1), temp(:,2),color(flag));
title('Before')
legend(legend_string1)

clear temp

if flag == 1
    disp('Processing Fish in Group1');
elseif flag==2
    disp('Processing Fish in Group2');
else
    disp('Processing Subject Fish');
end

%Go through each frame in each fish - First correct NaNs and then shootouts
for ii = 1:numFish
    
    
    %if first frame in NaN and then correct it to the first found value
    if isnan(traj(1,ii,1))
        ind1 = find(~(isnan(traj(:,ii,1))),1,'first');
        traj(1:ind1-1,ii,1) = squeeze(traj(ind1,ii,1));
        traj(1:ind1-1,ii,2) = squeeze(traj(ind1,ii,2));
    end
    
    %Check if it is NaN and correct with previous value
    for jj = 1:numFrames
        if isnan(traj(jj,ii,1))
            traj(jj,ii,:) = traj(jj-1,ii,:);
        end
    end
    
    %Check if first frame was misidentified and correct it to the first
    %correctly identified value
    if traj(1,ii,1) < area(1) || traj(1,ii,1) > area(2) || traj(1,ii,2) > area(3) || traj(1,ii,2) < area(4)
        ind1 = find(traj(:,ii,1) > area(1) & traj(:,ii,1) < area(2) & traj(:,ii,2) < area(3) & traj(:,ii,2) > area(4), 1, 'first');
        traj(1:ind1-1,ii,1) = squeeze(traj(ind1,ii,1));
        traj(1:ind1-1,ii,2) = squeeze(traj(ind1,ii,2));
    end
    
    %     %Find mean and std for each fish
    %     mean_x_fish(ii) = mean(traj(:,ii,1));
    %     std_x_fish(ii) = std(traj(:,ii,1));
    %     x_thresh1 = mean_x_fish(ii) + 2 * std_x_fish(ii);
    %     x_thresh2 = mean_x_fish(ii) - 2 * std_x_fish(ii);
    %
    %     %Find mean and std for each fish
    %     mean_y_fish(ii) = mean(traj(:,ii,2));
    %     std_y_fish(ii) = std(traj(:,ii,2));
    %     y_thresh1 = mean_y_fish(ii) +  2*std_y_fish(ii);
    %     y_thresh2 = mean_y_fish(ii) -  2*std_y_fish(ii);
    %
    
    %Check if there was huge jump between consecutive values and
    %correct with previous value
    for jj = 1:numFrames-1
        
        if (sqrt((traj(jj,ii,1)-traj(jj+1,ii,1)).^2+(traj(jj,ii,2)-traj(jj+1,ii,2)).^2))>threshold
            traj(jj+1,ii,:) = traj(jj,ii,:);
        end
        
        %One last check to remove anything above mean+/-2standarddeviation - for
        %the subject search in y only, for the groups search in x and y
        
        if traj(jj+1,ii,1) < area(1) || traj(jj+1,ii,1) > area(2) || traj(jj+1,ii,2) > area(3) || traj(jj+1,ii,2) < area(4)
            traj(jj+1,ii,:) = traj(jj,ii,:);
        end
        
        
        %         if flag == 1
        %             if traj(jj+1,ii,1) > 630 || traj(jj+1,ii,2) > y_thresh1 || traj(jj+1,ii,2) < y_thresh2
        %                 traj(jj+1,ii,:) = traj(jj,ii,:);
        %             end
        %         elseif flag == 3
        %             if traj(jj+1,ii,2) > 450
        %                 traj(jj+1,ii,:) = traj(jj,ii,:);
        %             end
        %         elseif flag == 2
        %             if traj(jj+1,ii,1) < x_thresh2 || traj(jj+1,ii,2) > y_thresh1 || traj(jj+1,ii,2) < y_thresh2
        %                 traj(jj+1,ii,:) = traj(jj,ii,:);
        %             end
        %         end
        
    end
end


modified_trajectories = traj;

fs = figure(2);
set(fs, 'color','white');
hold on
temp = reshape(modified_trajectories, size(modified_trajectories,1)*size(modified_trajectories,2),2);
plot(temp(:,1), temp(:,2),color(flag));
plot(mean(temp(:,1)), mean(temp(:,2)),'k*', 'MarkerSize',12)
title('After')
legend(legend_string2)
clear temp