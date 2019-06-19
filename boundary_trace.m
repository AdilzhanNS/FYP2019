function boundary = boundary_trace(fragment)
%find the first random whilte pixel (binary 1)
[r, c] = find(fragment, 1, 'first');
%create a 2D array by tracing the fragment
boundary = bwtraceboundary(fragment,[r c],'E',8,Inf,'clockwise'); %Moore-neighbourhood algorithm with Jacob's Stopping Criteria
end