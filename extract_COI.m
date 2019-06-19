function contour_extract = extract_COI(fragment, boundary)
%show the fragment with pixel info coordinates
imshow(fragment,[]);
hold on;
plot(boundary(:, 2),boundary(:,1),'g','LineWidth',1);
hold off;
impixelinfo;
%get the starting and ending pixel points 
prompt1 = 'Get the starting and ending coordinates of the contour of interest separately, please. \n Start with row-start and so on. \n Bear in mind that the contour tracing is clockwise.';
rstart = input(prompt1);
prompt2 = 'Get now the column start';
cstart = input(prompt2);
prompt3 = 'Get now the row end';
rend = input(prompt3);
prompt4 = 'Get now the column end';
cend = input(prompt4);
%create the contour of interest array based on the fragment contour
array_start = [rstart, cstart]; 
flag_idxs = 1;
for j=1:size(boundary(:,1),1)
    if (boundary(j:j,:)==array_start(1:1,:)) & (flag_idxs==1)
        idxstart = j;
        flag_idxs = 0;
    end
end
disp(idxstart);

array_end = [rend, cend];
flag_idxe = 1;
for j=1:size(boundary(:,1),1)
    if (boundary(j:j,:)==array_end(1:1,:)) & (flag_idxe==1)
        idxend = j;
        flag_idxe = 0;
    end
end
disp(idxend);

contour_extract = calcArray(boundary, idxstart, idxend);

call_prompt
    
end

function tempArray  = calcArray(array,startI,endI)
size_in_row = size(array,1);
size_in_col = size(array,2);
CircStartI = calcadd(startI,size_in_row);
CircEndI = calcadd(endI,size_in_row);
if(CircEndI >= CircStartI)
    
    tempArray = array(CircStartI:CircEndI,1:size_in_col);
elseif (CircEndI < CircStartI)
    
    tempArray(1:(size_in_row-CircStartI)+1,1:size_in_col) = array(CircStartI:size_in_row,1:size_in_col);
    tempArray = [tempArray;    array(1:CircEndI,1:size_in_col)];
end
end

function newcircadd = calcadd(add,size)
if (mod(add,size)~=0)
    newcircadd = mod(add,size);
else
    newcircadd = size;
end
end   

