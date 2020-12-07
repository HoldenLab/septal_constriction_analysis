% Author: Kevin Whitley
% Date created: 190710

% This function grabs trajectories of moving septa/bacteria/whatever from
% MicrobeJ data and extracts a small region around the centroids.

% v2 (201202):
%   - simplified input parameters
%   - initialized variables better
%   - cleaned some things up

function [subim, centx, centy, frames, olaystack_cropped, xy0_subim] = separate_microbej_tracks_v2(im, dat, num, param, olaystack)

if nargin==0
    param.s_box = 30; % size of box for trimming around septa
    param.n_frames_before = 8; % number of frames before start of trajectory to include (microbej won't get nascent septa)
    path = uigetdir('D:\Data_Theia');
    dirInf = dir([path '\*.mat']);
    dat = load([path '\' dirInf.name]);
end

n_frames = length(dat.Experiment.Lineage(num).Trajectory.Bacteria);

subim = zeros(param.s_box+1,param.s_box+1); newframes=subim; centx=zeros(n_frames,1); 
centy=centx; frames=centx; box_x1=centx; box_x2=centx; xdim=centx;
box_y1=centx; box_y2=centx; ydim=centx; xy0_subim=zeros(n_frames,2);
olaystack_cropped = [];

for kk = 1:n_frames
    centx(kk) = round(dat.Experiment.Lineage(num).Trajectory.Bacteria(kk).LOCATION.x);
    centy(kk) = round(dat.Experiment.Lineage(num).Trajectory.Bacteria(kk).LOCATION.y);
    frames(kk) = dat.Experiment.Lineage(num).Trajectory.Bacteria(kk).TRAJECTORY.frame.abs;
    
    % create bounding box (need to cut if near border)
    box_x1(kk) = max(centy(kk)-param.s_box/2, 1);
    box_x2(kk) = min(centy(kk)+param.s_box/2, size(im,1));
    xdim(kk) = box_x2(kk) - box_x1(kk);
    box_y1(kk) = max(centx(kk)-param.s_box/2, 1);
    box_y2(kk) = min(centx(kk)+param.s_box/2, size(im,2));
    ydim(kk) = box_y2(kk) - box_y1(kk);
    
    subim(1:xdim(kk)+1,1:ydim(kk)+1,kk) = im(box_x1(kk):box_x2(kk), box_y1(kk):box_y2(kk), frames(kk)+1);
    if ~isempty(olaystack)
        olaystack_cropped(1:xdim(kk)+1,1:ydim(kk)+1,kk) = olaystack(box_x1(kk):box_x2(kk), box_y1(kk):box_y2(kk), frames(kk)+1);
    end
    
    %define coordinates of microbej peak within bounding box
    y0 = centy(kk) - box_x1(kk) + 1;
    x0 = centx(kk) - box_y1(kk) + 1;
    xy0_subim(kk,:) = [x0, y0];
end

% add a few frames to the beginning to get more of nascent septa
fbefore = min([frames(1) param.n_frames_before]);
if fbefore > 0
    newframes(1:xdim(1)+1,1:ydim(1)+1,1:fbefore) = im(box_x1(1):box_x2(1), box_y1(1):box_y2(1), (frames(1)-fbefore+1):(frames(1)));
    newmat = newframes;
    newmat(:,:,(fbefore+1):(size(subim,3)+fbefore)) = subim;
    subim = newmat;
    if ~isempty(olaystack)
        newframesOLAY(1:xdim(1)+1,1:ydim(1)+1,1:fbefore) = olaystack(box_x1(1):box_x2(1), box_y1(1):box_y2(1), (frames(1)-fbefore+1):(frames(1)));
        newmatOLAY = newframesOLAY;
        newmatOLAY(:,:,(fbefore+1):(size(subim,3)+fbefore)) = olaystack_cropped;
        olaystack_cropped = newmatOLAY;
    end
end

