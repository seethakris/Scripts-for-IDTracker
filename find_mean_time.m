function [TimeSpent] = find_mean_time(X,Y,TimeSpent,flag)

NumFish = size(X,2);
temp_TimeSpent = zeros(size(TimeSpent,1),size(TimeSpent,2), NumFish);

if flag == 1
    disp('Processing Fish in Group1');
elseif flag==2
    disp('Processing Fish in Group2');
else
    disp('Processing Subject Fish');
end


for ii = 1:NumFish
       
    %Extract time spent per fish and then combine
    XY = [squeeze(X(:,ii)),squeeze(Y(:,ii))];
    
    for jj = 1:length(XY)
        temp_TimeSpent(round(XY(jj,1)), round(XY(jj,2)), ii) = temp_TimeSpent(round(XY(jj,1)), round(XY(jj,2)), ii) + 1;
    end
    
end

TimeSpent = TimeSpent + mean(temp_TimeSpent,3);