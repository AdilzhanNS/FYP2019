function bw = preprocessing(p) %bw stands for black-white, ie binary

im1 = p;
dim_check1 = size(im1); %check if the image needs to be converted to grayscale
[m_check1 n_check1] = size(dim_check1);
if (n_check1 > 2)
    grayim1 = rgb2gray(im1);
    bw = imbinarize(grayim1);
else
    bw = imbinarize(im1);
end
%check if the image needs to be complemented, ie reverse its 1s with 0s and
%other way around
check1_1 = size(find(bw(:,:)==1),1);
check0_1 = size(find(bw(:,:)==0),1);
if (check1_1 > check0_1)
    bw = imcomplement(bw);
end
%morphology to get rid of small holes inside or outside the fragment
bw = bwareaopen(bw, 50);
se = strel('disk',10);
bw = imclose(bw,se);

end



% im2 = q;
% dim_check2 = size(im2);
% [m_check2 n_check2] = size(dim_check2);
% if (n_check2 > 2)
%     grayim2 = rgb2gray(im2);
%     bw2 = imbinarize(grayim2);
% else
%     bw2 = imbinarize(im2);
% end
% 
% check1_2 = size(find(bw2(:,:)==1),1);
% check0_2 = size(find(bw2(:,:)==0),1);
% if (check1_2 > check0_2)
%     bw2 = imcomplement(bw2);
% end
% end




