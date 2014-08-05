function idTracker_Heatmap

[FileName, PathName] = uigetfile('*.mat', 'Select trajectories file');
Traj = load([PathName, FileName]);