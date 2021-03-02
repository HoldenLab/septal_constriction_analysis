function imstack = imreadstack( path1)
%	imstack = imreadstack( path1)
% 
% Description:
% 	Automatically load all the images in a stack as a 
% 3D matrix
%
% Inputs:		
% 	path
%
% Outputs:
%	imstack
%
% History:
% 	080213 	- First alpha complete (SH)
%	080213 	- Function headers and help standardised (SH)
%	
% Author:
%	Seamus Holden
%

%number of frames in image and size of image
numFrames = numel(imfinfo(path1));
I = imread( path1, 1);

%preallocate array
imstack = zeros( [size(I) numFrames], class(I));
imstack(:,:,1) = I;

for p = 2:numFrames
	imstack(:,:,p) = imread(path1, p);
end 

