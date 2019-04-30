function [concat_img] = img_vertcat_whitepad(X1,X2)

% concatenates two images vertically that differ in their horizontal size
% by padding the smaller image with white space
% X1: image 1, larger image
% X2: imgae 2, smaller image

total_white = size(X1,2) - size(X2,2);
l_white = floor(0.5*total_white);
r_white = ceil(0.5*total_white);

white = permute([255 255 255],[1 3 2]); % make a white pixel
l_white = repmat(white,[size(X2,1) l_white 1]); % tile it
r_white = repmat(white,[size(X2,1) r_white 1]); % tile it
X2 = [l_white,X2,r_white];

concat_img = [X1;X2];