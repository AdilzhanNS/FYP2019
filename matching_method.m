function [subset_contour_final, contour_final, min_check, wbox] = matching_method(contour1, contour2)
clc;
w=10;
workspace;  % Make sure the workspace panel is showing.
format long g;
format compact;
fontSize = 12;

tic; %start the stopwatch

%the following big "if" statement has three conditions: when contour of
%either is smaller than the other one and when they are exactly the same
%length in terms of the contour

if size(contour1(:,1),1)< size(contour2(:,1),1) 
    subset_contour = contour1; %create subset contour which will be windowing through the larger contour
    %xC1=sum(subset_contour(:,2))/size(subset_contour(:,2),1);
    %yC1=sum(subset_contour(:,1))/size(subset_contour(:,1),1);
    xC1 = mean(subset_contour(:,2));
    yC1 = mean(subset_contour(:,1));
    x1 = subset_contour(:, 2);
    y1 = subset_contour(:, 1);
    
    k1 = zeros(size(subset_contour(:,1),1),1); %preallocate memory
    distance1 = zeros(size(subset_contour(:,1),1),1); %preallocate memory
    for i = 1:size(subset_contour(:,1),1)
        distance1(i,1) = sqrt((x1(i)-xC1).^2 + (y1(i)-yC1).^2);
        k1(i,1) = i;
    end %all of the above is to find the distances of the pixels to the centroid of the smaller contour
    subplot(2,3,2);
    plot(transpose(k1), transpose(distance1));
    title('Diagram of pixel distances to centroid of Sample 1'); 
    differ_size = size(contour2(:,1),1)-size(subset_contour(:,1),1); %the size difference between two contours which is used at the next "for" loop as an iterative indicator
   
    %the following "for" loop essentially now checks each subset location
    %starting from the beginning of the larger contour and taking the size
    %equal to the smaller contour and performs operation of checking the
    %error scatter distribution between each iterative subset of the larger
    %contour and the smaller subset. The way I see it basically is that if
    %the shapes are the same then the summation of all the distances to the
    %centroid have to be ideally the same regardless of the position of
    %each contour/subset of the larger contour. This loop should break when
    %their sums of distances are the same (well that is for the ideal case,
    %in non-ideal there should be some small tolerance threshold which will
    %allow for the matching to still happen)
    count2 = 0;
    error_scatter_test = zeros(differ_size+1,1); %preallocate memory
    start_contour2 = zeros(differ_size+1,2); %preallocate memory
    start_subset_contour = zeros(differ_size+1,2); %preallocate memory
    size_check = zeros(differ_size+1,1); %preallocate memory
    for j = 0:differ_size
        contour2_test = contour2(j+1:calcadd(size(subset_contour(:,1),1) + j,size(contour2(:,1),1)), 1:2);
        xC2=sum(contour2_test(:,2))/size(contour2_test(:,2),1);
        yC2=sum(contour2_test(:,1))/size(contour2_test(:,1),1);
        x2 = contour2_test(:, 2);
        y2 = contour2_test(:, 1);
        k2 = zeros(size(contour2_test(:,1),1),1);
        distance2 = zeros(size(contour2_test(:,1),1),1);
        for i = 1:size(contour2_test(:,1),1)
            distance2(i,1) = sqrt((x2(i)-xC2).^2 + (y2(i)-yC2).^2);
            k2(i,1) = i;
        end
        %distance metric, precisely L_p distance, is the metric_err
        %metric_err = 0;
        %for h = 1:size(contour2_test(:,1),1)
        %    diff_dist = ((distance2(h,1) - distance1(h,1)).^t);
        %    metric_err = metric_err + diff_dist;
        %end
        %metric_err_fin = (metric_err).^(1/t);
        
        %m1 = mean(distance1);
        %disp(m1);
        %m2 = mean(distance2);
        %distance11 = distance1 - m1;
        %distance22 = distance2 - m2;
        %sum12 = distance11 + distance22;
        %sum123 = sum(sum12(:,1));
        %distance_metric = pdist2(transpose(distance2(:,1)),transpose(distance1(:,1)), 'euclidean');
        %distance_metric = pdist2(transpose(distance2(:,1)),transpose(distance1(:,1)), 'seucliden');
        %distance_metric = pdist2(transpose(distance2(:,1)),transpose(distance1(:,1)), 'minkowski',3);
        %distance_metric = pdist2(transpose(distance2(:,1)),transpose(distance1(:,1)), 'minkowski',4);
        %distance_metric = pdist2(transpose(distance2(:,1)),transpose(distance1(:,1)), 'chebyshev');
        %disp(distance_metric);
        %error_scatter_test(j+1,1) = metric_err_fin;
        
        [h,p,ks2stat] = kstest2(distance1(:,1),distance2(:,1)); %one-dimensional two sample Kolmogorov-Smirnov Test
        disp(ks2stat); %test statistic
        error_scatter_test(j+1,1) = ks2stat; %store the test statistic
        start_contour2(j+1:j+1,:) = contour2_test(1:1,:);
        start_subset_contour(j+1:j+1,:) = subset_contour(1:1,:);
        size_check(j+1, 1) = size(subset_contour(:,1),1);
    
        disp('executed for when contour2 is larger');
        count2 = count2 + 1;
        disp(count2);
    end
    cell_check = {error_scatter_test start_contour2 start_subset_contour size_check};
    disp(cell_check);    
    
    indx = find([cell_check{:}] == min(error_scatter_test));
    if (size(indx,1) > 1)
        ind = indx(1,1);
    else
        ind=indx;
    end
    %the following displays are for the debugging purposes
    disp(ind); 
    disp(error_scatter_test(ind));
    disp(start_contour2(ind,:));
    disp(start_subset_contour(ind,:));
    disp(size_check(ind));
    
    disp(size(error_scatter_test));
    min_check = min(error_scatter_test);
    disp(min(error_scatter_test));
    disp(max(error_scatter_test));
    disp(abs(1 - (max(error_scatter_test)-min(error_scatter_test))/max(error_scatter_test))*100);
    %the threshold to pass the matching 
    if (min(error_scatter_test)<0.0215)
        s=1;
        e=1;
        while s
            if(subset_contour(e,:)==start_subset_contour(ind,:))
                s=0;
                ind_row_sub_1 = e;
            end
            e=e+1;
        end
        disp(ind_row_sub_1);
        subset_contour_final = contour1(ind_row_sub_1:calcadd(size(subset_contour(:,1),1) + ind_row_sub_1 - 1,size(subset_contour(:,1),1)),1:2);
        a=1;
        ee = 1;
        while a
            if(contour2(ee,:) == start_contour2(ind,:))
                a=0;
                ind_row_2 = ee;
            end
            ee=ee+1;
        end
        disp(ind_row_2);
        disp(contour2(ind_row_2,:));
        contour_final = contour2(ind_row_2:size(subset_contour(:,1),1) + ind_row_2 - 1, 1:2);
        disp(size(contour_final));

        wbox = msgbox('contour2 is longer here and thus the contour2-related fragment takes contour_final values');
    else
        wbox = msgbox('These contours are most likely not matching');
    end
    
elseif size(contour1(:,1),1)> size(contour2(:,1),1) %this part of the bigger "if" statement should do the same operation as the aforementiond part but now the other contour is smaller than the other
   
    subset_contour = contour2;
    xC2 = mean(subset_contour(:,2));
    yC2 = mean(subset_contour(:,1));
    x2 = subset_contour(:, 2);
    y2 = subset_contour(:, 1);
    
    k2 = zeros(size(contour2(:,1),1),1); %preallocate memory
    distance2 = zeros(size(contour2(:,1),1),1); %preallocate memory
    for i = 1:size(contour2(:,1),1)
        distance2(i,1) = sqrt((x2(i)-xC2).^2 + (y2(i)-yC2).^2);
        k2(i,1) = i;
    end
    subplot(2,3,4);
    plot(transpose(k2), transpose(distance2));
    title('Diagram of pixel distances to centroid of Sample 2');
    differ_size = size(contour1(:,1),1)-size(subset_contour(:,1),1);
    count1 = 0;
    error_scatter_test = zeros(differ_size+1,1); %preallocate memory
    start_contour1 = zeros(differ_size+1,2); %preallocate memory
    start_subset_contour = zeros(differ_size+1,2); %preallocate memory
    size_check = zeros(differ_size+1,1); %preallocate memory
    for j = 0:differ_size 
        contour1_test = contour1(j+1:calcadd(size(subset_contour(:,1),1) + j,size(contour1(:,1),1)), 1:2);
        xC1=sum(contour1_test(:,2))/size(contour1_test(:,2),1);
        yC1=sum(contour1_test(:,1))/size(contour1_test(:,1),1);

        x1 = contour1_test(:, 2);
        y1 = contour1_test(:, 1);
        
        k1 = zeros(size(contour1_test(:,1),1),1); %preallocate memory
        distance1 = zeros(size(contour1_test(:,1),1),1); %preallocate memory
        for i = 1:size(contour1_test(:,1),1)
            distance1(i,1)= sqrt((x1(i)-xC1).^2 + (y1(i)-yC1).^2);
            k1(i,1) = i;
        end
        %distance metric, precisely L_p distance, is the metric_err   
        %metric_err = 0;
        %for h = 1:size(contour1_test(:,1),1)
        %    diff_dist = ((distance1(h,1) - distance2(h,1)).^t);
        %    metric_err = metric_err + diff_dist;
        %end
        
        %distance_metric = pdist2(transpose(distance2(:,1)),transpose(distance1(:,1)), 'Euclidean');
        %distance_metric = pdist2(transpose(distance2(:,1)),transpose(distance1(:,1)), 'seucliden');
        %distance_metric = pdist2(transpose(distance2(:,1)),transpose(distance1(:,1)), 'minkowski',3);
        %distance_metric = pdist2(transpose(distance2(:,1)),transpose(distance1(:,1)), 'minkowski',4);
        %distance_metric = pdist2(transpose(distance2(:,1)),transpose(distance1(:,1)), 'chebyshev');
        %disp(distance_metric);
        %[f1,x1] = ecdf(transpose(distance1(:,1)));
        %[f2,x2] = ecdf(transpose(distance2(:,1)));
        
        %m1 = mean(distance1(:,1));
        %disp(m1);
        %m2 = mean(distance2(:,1));
        %distance11 = distance1 - m1;
        %distance22 = distance2 - m2;
        %sum12 = distance11 + distance22;
        %sum123 = sum(sum12(:,1));
        %disp(sum123);
        %metric_err_fin = (metric_err).^(1/t);
        %error_scatter_test(j+1,1) = metric_err_fin;
        
        
        [h,p,ks2stat] = kstest2(distance1(:,1),distance2(:,1)); %one-dimensional two sample Kolmogorov-Smirnov Test
        disp(ks2stat); %test statistic
        error_scatter_test(j+1,1) = ks2stat; %store the test statistic

        start_contour1(j+1:j+1,:) = contour1_test(1:1,:);
        start_subset_contour(j+1:j+1,:) = subset_contour(1:1,:);
        size_check(j+1, 1) = size(subset_contour(:,1),1);
        disp('executed for when contour1 is larger');
        count1 = count1 + 1;
        disp(count1);
    end
    cell_check = {error_scatter_test start_contour1 start_subset_contour size_check};
    %the following displays are for debugging purposes
    disp(cell_check);
    indx = find([cell_check{:}] == min(error_scatter_test));
    if (size(indx,1) > 1)
        ind = indx(1,1);
    else
        ind=indx;
    end
    disp(ind);
    disp(error_scatter_test(ind));
    disp(start_contour1(ind,:));
    disp(start_subset_contour(ind,:));
    disp(size_check(ind));
    
    disp(size(error_scatter_test));
    min_check = min(error_scatter_test);
    disp(min(error_scatter_test));
    disp(max(error_scatter_test));
    disp(abs(1 - (max(error_scatter_test)-min(error_scatter_test))/max(error_scatter_test))*100); %the differene between the minimum error and large error

    %the threshold to pass the matching 
    if (min(error_scatter_test)<0.0215) 
        s=1;
        e=1;
        while s
            if(subset_contour(e,:)==start_subset_contour(ind,:))
                s=0;
                ind_row_sub_2 = e;
            end
            e=e+1;
        end
        disp(ind_row_sub_2);
        subset_contour_final = contour2(ind_row_sub_2:calcadd(size(subset_contour(:,1),1) + ind_row_sub_2 - 1,size(subset_contour(:,1),1)),1:2);
        a=1;
        ee = 1;
        while a
            if(contour1(ee,:) == start_contour1(ind,:))
                a=0;
                ind_row_1 = ee;
            end
            ee=ee+1;
        end
        disp(ind_row_1);
        disp(contour1(ind_row_1,:));
        contour_final = contour1(ind_row_1:size(subset_contour(:,1),1) + ind_row_1 - 1, 1:2);
        disp(size(contour_final));
        wbox = msgbox('contour1 is longer here and thus the contour1-related fragment takes contour_final values');
    else
        wbox = msgbox('These contours are most likely not matching');
    end
    
elseif  isequal(size(contour1(:,1),1),size(contour2(:,1),1)) %this part of "if" statement kicks in when the sizes are the same and straight up just performs all the operations on both contours 
    xC1=sum(contour1(:,2))/size(contour1(:,2),1);
    yC1=sum(contour1(:,1))/size(contour1(:,1),1);
    x1 = contour1(:, 2);
    y1 = contour1(:, 1);
%     labeledImage1 = bwlabel(bw1, 8);
%     measurements1 = regionprops(labeledImage1, 'Centroid');
%     xC1 = measurements1.Centroid(1);
%     yC1 = measurements1.Centroid(2);
%     x1 = contour1(:, 2);
%     y1 = contour1(:, 1);
    k1 = zeros(size(contour1(:,1),1),1);
    distance1 = zeros(size(contour1(:,1),1),1);
    for i = 1:size(contour1(:,1),1)
        distance1(i,1)= sqrt((x1(i)-xC1).^2 + (y1(i)-yC1).^2);
        k1(i,1) = i;

    end
    subplot(2,3,2);
    plot(transpose(k1), transpose(distance1));
    title('Diagram of pixel distances to centroid of Sample 1');
    
    xC2=sum(contour2(:,2))/size(contour2(:,2),1);
    yC2=sum(contour2(:,1))/size(contour2(:,1),1);
    x2 = contour2(:, 2);
    y2 = contour2(:, 1);
%     labeledImage2 = bwlabel(bw2, 8);
%     measurements2 = regionprops(labeledImage2, 'Centroid');
%     xC2 = measurements2.Centroid(1);
%     yC2 = measurements2.Centroid(2);
%     x2 = contour2(:, 2);
%     y2 = contour2(:, 1);
    k2 = zeros(size(contour2(:,1),1),1);
    distance2 = zeros(size(contour2(:,1),1),1);
    for i = 1:size(contour2(:,1),1)
        distance2(i,1)= sqrt((x2(i)-xC2).^2 + (y2(i)-yC2).^2);
        k2(i,1) = i;
    end
    subplot(2,3,4);
    plot(transpose(k2), transpose(distance2));
    title('Diagram of pixel distances to centroid of Sample 2');
    
    %metric_err = 0;
    %for h = 1:size(contour2(:,1),1)
    %    diff_dist = ((distance2(h,:) - distance1(h,:)).^t);
    %    metric_err = metric_err + diff_dist;
    %end
    %metric_err_fin = (metric_err).^(1/t);
    %error_scatter = metric_err_fin;
    
    [h,p,ks2stat] = kstest2(distance1(:,1),distance2(:,1)); %one-dimensional two sample Kolmogorov-Smirnov Test
    disp(ks2stat); %test statistic
    error_scatter = ks2stat; %store the test statistic

    if (error_scatter<0.0215) 
        subset_contout_final = contour1;
        contour_final = contour2;
    end

end
toc;
end
    

% function tempArray  = calcArray(array,startI,endI)
% size_in_row = size(array,1);
% size_in_col = size(array,2);
% CircStartI = calcadd(startI,size_in_row);
% CircEndI = calcadd(endI,size_in_row);
% if(CircEndI >= CircStartI)
%     
%     tempArray = array(CircStartI:CircEndI,1:size_in_col);
% elseif (CircEndI < CircStartI)
%     
%     tempArray(1:(size_in_row-CircStartI)+1,1:size_in_col) = array(CircStartI:size_in_row,1:size_in_col);
%     tempArray = [tempArray;    array(1:CircEndI,1:size_in_col)];
% end
% end

function newcircadd = calcadd(add,size)
if (mod(add,size)~=0)
    newcircadd = mod(add,size);
else
    newcircadd = size;
end
end
