% Author: Kevin Whitley
% Date created: 190710

% This function grabs trajectories of moving septa/bacteria/whatever from
% MicrobeJ data and extracts a small region around the centroids.

s_box = 30/2; % size of box for trimming around septa
n_frames_before = 8; % number of frames before start of trajectory to include (microbej won't get nascent septa)

path = uigetdir('D:\Data_Theia');
dirInf = dir([path '\*.mat']);

dat = load([path '\' dirInf.name]);
ncells = length(dat.Experiment.Lineage);

% for jj = 1:ncells
    for jj = 14
    if length(dat.Experiment.Lineage(jj).Trajectory)==1
        subim = zeros(s_box*2+1,s_box*2+1); newframes = subim; centx = []; centy = []; frames = [];
        box_x1 = []; box_x2 = []; xdim = []; box_y1 = []; box_y2 = []; ydim = [];
        
        for kk = 1:length(dat.Experiment.Lineage(jj).Trajectory.Bacteria)
            centx(kk,jj) = round(dat.Experiment.Lineage(jj).Trajectory.Bacteria(kk).LOCATION.x);
            centy(kk,jj) = round(dat.Experiment.Lineage(jj).Trajectory.Bacteria(kk).LOCATION.y);
            frames(kk,jj) = dat.Experiment.Lineage(jj).Trajectory.Bacteria(kk).TRAJECTORY.frame.abs;
            
            % create bounding box (need to cut if near border)
            box_x1(kk) = max(centy(kk,jj)-s_box, 1);
            box_x2(kk) = min(centy(kk,jj)+s_box, size(im,1));
            xdim(kk) = box_x2(kk) - box_x1(kk);
            box_y1(kk) = max(centx(kk,jj)-s_box, 1);
            box_y2(kk) = min(centx(kk,jj)+s_box, size(im,2));
            ydim(kk) = box_y2(kk) - box_y1(kk);
            
            subim(1:xdim(kk)+1,1:ydim(kk)+1,kk) = im(box_x1(kk):box_x2(kk), box_y1(kk):box_y2(kk), frames(kk,jj)+1);
        end
        
        % add another frame to the beginning to get more of nascent septa
        fbefore = min([frames(1,jj) n_frames_before]);
        if fbefore > 0
            newframes(1:xdim(1)+1,1:ydim(1)+1,1:fbefore) = im(box_x1(1):box_x2(1), box_y1(1):box_y2(1), (frames(1,jj)-fbefore+1):(frames(1,jj)));
            newmat = newframes;
            newmat(:,:,(fbefore+1):(size(subim,3)+fbefore)) = subim;
            subim = newmat;
        end
        
        % one more frame for luck (in case microbej didn't do a great job at
        % the end)
        if frames(end,jj)+2 <= size(im,3)
            subim(1:xdim(end)+1,1:ydim(end)+1,end+1) = im(box_x1(end):box_x2(end), box_y1(end):box_y2(end), frames(end,jj)+2);
        end
        
        subim1 = cast(round(subim),'uint16');
        imwritestack([dirInf.folder '\' dirInf.name(1:end-4) 'lin' num2str(jj) '.tif'], subim1)
    end
end