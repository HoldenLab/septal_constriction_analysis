% Author: Kevin Whitley
% Date created: 190710

% This function grabs trajectories of moving septa/bacteria/whatever from
% MicrobeJ data and extracts a small region around the centroids.

function [subim, centx, centy, frames, olaystack_cropped, xy0_subim] = separate_microbej_tracks(im, dat, num, s_box, n_frames_before, olaystack)

if nargin==0
    s_box = 30; % size of box for trimming around septa
    n_frames_before = 8; % number of frames before start of trajectory to include (microbej won't get nascent septa)
    path = uigetdir('D:\Data_Theia');
    dirInf = dir([path '\*.mat']);
    dat = load([path '\' dirInf.name]);
end

olaystack_cropped = [];
if length(dat.Experiment.Lineage(num).Trajectory)==1
    subim = zeros(s_box+1,s_box+1); newframes = subim; centx = []; centy = []; frames = [];
    box_x1 = []; box_x2 = []; xdim = []; box_y1 = []; box_y2 = []; ydim = [];
    xy0_subim=[];
    
    for kk = 1:length(dat.Experiment.Lineage(num).Trajectory.Bacteria)
        centx(kk) = round(dat.Experiment.Lineage(num).Trajectory.Bacteria(kk).LOCATION.x);
        centy(kk) = round(dat.Experiment.Lineage(num).Trajectory.Bacteria(kk).LOCATION.y);
        frames(kk) = dat.Experiment.Lineage(num).Trajectory.Bacteria(kk).TRAJECTORY.frame.abs;
        
        % create bounding box (need to cut if near border)
        box_x1(kk) = max(centy(kk)-s_box/2, 1);
        box_x2(kk) = min(centy(kk)+s_box/2, size(im,1));
        xdim(kk) = box_x2(kk) - box_x1(kk);
        box_y1(kk) = max(centx(kk)-s_box/2, 1);
        box_y2(kk) = min(centx(kk)+s_box/2, size(im,2));
        ydim(kk) = box_y2(kk) - box_y1(kk);
        
        subim(1:xdim(kk)+1,1:ydim(kk)+1,kk) = im(box_x1(kk):box_x2(kk), box_y1(kk):box_y2(kk), frames(kk)+1);
        if ~isempty(olaystack)
            olaystack_cropped(1:xdim(kk)+1,1:ydim(kk)+1,kk) = olaystack(box_x1(kk):box_x2(kk), box_y1(kk):box_y2(kk), frames(kk)+1);
        end
        
        %define coordinates of microbej peak within bounding box
        y0=centy(kk)-box_x1(kk)+1;
        x0=centx(kk)-box_y1(kk)+1;
        xy0_subim(kk,:)=[x0,y0];
    end
    
    % add a few frames to the beginning to get more of nascent septa
    fbefore = min([frames(1) n_frames_before]);
    if fbefore > 0
        newframes(1:xdim(1)+1,1:ydim(1)+1,1:fbefore) = im(box_x1(1):box_x2(1), box_y1(1):box_y2(1), (frames(1)-fbefore+1):(frames(1)));
        newmat = newframes;
        newmat(:,:,(fbefore+1):(size(subim,3)+fbefore)) = subim;
        subim = newmat;
        if ~isempty(olaystack)
            newframesOLAY(1:xdim(1)+1,1:ydim(1)+1,1:fbefore) = olaystack(box_x1(1):box_x2(1), box_y1(1):box_y2(1), (frames(1)-fbefore+1):(frames(1)));
            newmatOLAY = newframesOLAY;
            newmatOLAY(:,:,(fbefore+1):(size(subim,3)+fbefore)) = olaystack_cropped;
            olaystack_cropped= newmatOLAY;
        end
    end
    
    % one more frame for luck (in case microbej didn't do a great job at
    % the end). deprecated 190802
%     if frames(end,num)+2 <= size(im,3)
%         subim(1:xdim(end)+1,1:ydim(end)+1,end+1) = im(box_x1(end):box_x2(end), box_y1(end):box_y2(end), frames(end,num)+2);
%     end
    
    % deprecated 190802
    %subim1 = cast(round(subim),'uint16');
    %imwritestack([dirInf.folder '\' dirInf.name(1:end-4) 'lin' num2str(num) '.tif'], subim1)
end