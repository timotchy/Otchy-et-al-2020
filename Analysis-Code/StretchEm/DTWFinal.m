function [DTWdist p] = DTWFinal(F1,F2,Band)
% This function implements a Dynamic Time-Warping (DTW) algorithm. Inputs are two feature matrices
% F1 and F2 (In this application, typically spectrograms). F1 is and L-by-M1 matrix; F2 is and L-by-M2 matrix. 
% Each column represents a feature at a particular time step.
%
% F1 is the template; F2 will be warped to fit F1.
% Band sets the maximum distance (in columns of F1 and F2) that the alorgithm will seek for a best match.
%
%The outputs specify how much warping was done (DTWdist) and the warping path (p).
%
%
%Updated with comments 11/11/2015 TMO

%Check the Band value; if empty or missing, substitute a default (125 ms)
if nargin < 3 || isempty(Band)
    Band = 125;
end

%Confirm that the matrices have the same number of rows/features
[L1, M1] = size(F1);
[L2, M2] = size(F2);
if L1 ~= L2
    disp('Size mismatch; alignment aborted.');
    return;
end


% Accumulated Cost Matrix

%Initialize cost matrices with maximum cost/distance
D = inf*ones(M1,M2);
dM = D;

%Fill in the distance matrix (dM) for matching columns of F1 and F2
% Restrict calculation to the range allowed by Band; all other comparisons remain at inf
for i =1:M1
    for j = max(1,i-Band-1):min(M2,i+Band+1)
        %Calculate the euclidean distance between columns (i.e. timepoints) in F1 and F2 
        dM(i,j) = norm(F1(:,i)-F2(:,j));
    end
end

%Generate an accumulated distance matrix that sums along the entries of dM. There is some history dependent/edge effects to
%deal with, so the first two rows and columns of the accumlation matrix are coded separately.

%Update the first column (up to the length of the band)
for i = 1:min(M1,Band+1)
    %Initialize at 0
    D(i,1) = 0;
    
    %Sum all pairwise distances (in dM) up to the row you're currently counting
    for j = 1:i
        D(i,1) = D(i,1)+dM(j,1);
    end
end

%Update the first row (up to the length of the band) as just a straight copy from dM
for i = 1:min(M2,Band+1)   
    D(1,i) = dM(1,i); % Is the loop really faster than just a copy?
end

%Update the second column of the accumlation matrix
for i = 2:min(M1,Band+2)
        dij = dM(i,2); %Distance between current column of F1 and second column of F2
        D(i,2) = min(D(i-1,1)+2*dij, min(D(i-1,2)+dij, D(i,1)+dij)); %Pick shortest of the three distances nearest the current cell (1-cell memory b/c you're at the edge of the matrix)
end

%Update the second row of the accumlation matrix
for j = 2:min(M2,Band+2)
        dij = dM(2,j); %Distance between second column of F1 and current column of F2
        D(2,j) = min(D(1,j-1)+2*dij,min(D(1,j)+dij,D(2,j-1)+dij)); %Pick shortest of the three distances nearest the current cell (2-cell memory b/c you're at the edge of the matrix)
end
%Update the rest of the accumulation matrix (within the limits of the band) using the same basic rules
for j = 3:M2
    for i = max(3,j-Band):min(M1,j+Band)
        
        dij = dM(i,j);              %Distance between current column of F1 and current column of F2
        dim1j = dM(i-1,j);       %Distance between previous column of F1 and current column of F2
        dijm1 = dM(i,j-1);       %Distance between current column of F1 and previous column of F2
        
        %Pick shortest of the three distances nearest the current cell (n-cell memory up to the band limit)
        D(i,j) = min(2*dim1j+dij+D(i-2,j-1),min(2*dij+D(i-1,j-1), 2*dijm1+dij+D(i-1,j-2))); 
    end
end


% Implement a variable stepsize

%Find the optimal warping path by starting at the end of the matrix and working backward

%Locate the smallest value in the very last row of the accumulation matrix (up the the band limit)
[a b] = min(D(M1,(M2-Band):M2));
b = b+M2-Band-1; %convert to the actual index value
p(1,:) = [M1 b]; %The first point in the warping path matches the last point in M1 with the closest point nearby

%Initialize the pointer
count = 2;
while 1
    %The last set of matching points that were determined
    curindex = p(count-1,:);
    
    %If we're now at the first time step of F1, then we're done and we can end the while-loop
    if curindex(1) == 1
        break;
        
    %If we're on the second step of F1 and NOT the first of F2, OR if we're on the second step of F2    
    elseif ( curindex(1) == 2 && curindex(2)~= 1 ) || curindex(2) == 2
            
            %Accum distances along the three closest paths (with 1-cell memory b/c we're near the edge)
            dum1 = D(curindex(1)-1,curindex(2)-1);  %Accumulated distance on a move along the diagonal
            dum2 = D(curindex(1)-1,curindex(2));     %Accumulated distance on a move along F1
            dum3 = D(curindex(1),curindex(2)-1);     %Accumulated distance on a move along F2
            
            if (dum1<dum2 && dum1<dum3)
                %If the diagonal move makes the shortest path, update path coordinates appropriately
                p(count,:) = [curindex(1)-1 curindex(2)-1];
                
            elseif (dum2<=dum1 && dum2<dum3)
                %If move along F1 makes the shortest path, update path coordinates appropriately
                p(count,:) = [curindex(1)-1 curindex(2)];
                
            else
                %If move along F2 makes the shortest path, update path coordinates appropriately
                p(count,:) = [curindex(1) curindex(2)-1] ;
                
            end
            
        %if we're on the first step of F2
        elseif curindex(2) == 1
            p(count,:) = [curindex(1)-1 curindex(2)];
        
        %Otherwise...
        else
            %Accum distances along the three closest paths (with 2-cell memory b/c we have the room to run)
            dum1 = D(curindex(1)-1,curindex(2)-1);      %Accumulated distance on a move along the diagonal (no stretch)
            dum2 = D(curindex(1)-1,curindex(2)-2);      %Accumulated distance on a stretch along F1
            dum3 = D(curindex(1)-2,curindex(2)-1);      %Accumulated distance on a stretch along F2
            
            if (dum1<dum2 && dum1<dum3)
                %If the diagonal move makes the shortest path, update path coordinates appropriately
                p(count,:) = [curindex(1)-1 curindex(2)-1];
                
            elseif (dum2<=dum1 && dum2<dum3)
                %If move along F1 makes the shortest path, update path coordinates appropriately
                p(count,:) = [curindex(1) curindex(2)-1];
                count = count+1; %account for stretch
                p(count,:) = [curindex(1)-1 curindex(2)-2];
            else
                %If move along F2 makes the shortest path, update path coordinates appropriately
                p(count,:) = [curindex(1)-1 curindex(2)];
                count = count+1; %account for stretch
                p(count,:) = [curindex(1)-2 curindex(2)-1];
            end
    end
    
    %Update the point for the next go around
    count = count+1;
end

% Retrieve the total accumulated distance of the path chosen
DTWdist = D(M1,b);

%Since we worked from the back forward, we nnow need to flip the path vector to read from start to finish
p=flipud(p);

end

