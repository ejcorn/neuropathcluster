function cropped_image = crop_brain_image_manual(cdata,crop_idx)

% crop_idx [xmin xmax ymin ymax]

if ~exist('marg','var')
    crop_idx = [45 375 35 157];
end
crop_left = crop_idx(1);
crop_right = crop_idx(2);
crop_top = crop_idx(3);
crop_bottom = crop_idx(4);

cropped_image = cdata(crop_top:crop_bottom,crop_left:crop_right,:);


