function [concat_img] = img_horzcat_whitepad(X1,X2)

% concatenates two images horizontally that differ in their vertical size
% by padding the smaller image with white space
% X1: image 1, larger image
% X2: imgae 2, smaller image

total_white = size(X1,1) - size(X2,1);
t_white = floor(0.5*total_white);
b_white = ceil(0.5*total_white);

white = permute([255 255 255],[1 3 2]); % make a white pixel
t_white = repmat(white,[t_white size(X2,2) 1]); % tile it
b_white = repmat(white,[b_white size(X2,2) 1]); % tile it
X2 = [t_white;X2;b_white];

concat_img = [X1,X2];