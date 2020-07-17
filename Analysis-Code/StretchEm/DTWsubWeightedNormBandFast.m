%
% This function implements DTW algorithm. Inputs are two feature matrices
% F1 and F2. F1 is and L1-by-M1 matrix, each column represents the feature 
% at a particular time step. F1 is the template.
%
% This is a classical DTW. No extra conditions are assigned.
%


function [DTWdist p] = DTWsubWeightedNormBandFast(F1,F2, Band)

if nargin < 3 || isempty(Band)
    Band = 125;
end

[L1 M1] = size(F1);
[L2 M2] = size(F2);

Band = 200;

if L1 ~= L2
    disp('Size mismatch');
    return;
end

% Accumulated Cost Matrix


D = inf*ones(M1,M2);

for i =1:M1
    for j = 1:M2
   
        dM(i,j) = norm(F1(:,i)-F2(:,j))^2/norm(F1(:,i))/norm(F2(:,j));
    end
end

for i = 1:min(M1,Band+1)
    
    D(i,1) = 0;
    
    for j = 1:i
        D(i,1) = D(i,1)+dM(j,1);
    end
end


for i = 1:min(M2,Band+1)   
    D(1,i) = dM(1,i);
end

for i = 2:min(M1,Band+2)
        dij = dM(i,2);
        D(i,2) = min(D(i-1,1)+dij,min(D(i-1,2)+dij/2,D(i,1)+dij/2));
end

for j = 2:min(M2,Band+2)
        dij = dM(2,j);
        D(2,j) = min(D(1,j-1)+dij,min(D(1,j)+dij/2,D(2,j-1)+dij/2));
end


 
for j = 3:M2
    for i = max(3,j-Band):min(M1,j+Band)
        
        dij = dM(i,j);
        dim1j = dM(i-1,j);
        dim1jm1=dM(i-1,j-1);
        dijm1 = dM(i,j-1);
        
        D(i,j) = min(3/4*dij+3/8*dim1j+3/8*dim1jm1+D(i-2,j-1),min(dij+D(i-1,j-1),3/4*dij+3/8*dijm1+3/8*dim1jm1+D(i-1,j-2)));
    end
end


% Implement a variable stepsize



% OPTIMAL WARPING PATH

% 


[a b] = min(D(M1,(M2-Band):M2));
b = b+M2-Band-1;
p(1,:) = [M1 b];

count = 2;

while 1
    
    curindex = p(count-1,:);
    
    
    if curindex(1) == 1
            break;
    else if ( curindex(1) == 2 && curindex(2)~= 1 ) || curindex(2) == 2 
            
            dum1 = D(curindex(1)-1,curindex(2)-1);
            dum2 = D(curindex(1)-1,curindex(2));
            dum3 = D(curindex(1),curindex(2)-1);
                
            if (dum1<dum2 && dum1<dum3)
                    p(count,:) = [curindex(1)-1 curindex(2)-1];
            else if (dum2<=dum1 && dum2<dum3)
                    p(count,:) = [curindex(1)-1 curindex(2)];
                 else
                    p(count,:) = [curindex(1) curindex(2)-1] ;
                 end
            end
          
            
         else if curindex(2) == 1
                p(count,:) = [curindex(1)-1 curindex(2)];
              else
                
                dum1 = D(curindex(1)-1,curindex(2)-1);
                dum2 = D(curindex(1)-1,curindex(2)-2);
                dum3 = D(curindex(1)-2,curindex(2)-1);
                
                if (dum1<dum2 && dum1<dum3)
                    p(count,:) = [curindex(1)-1 curindex(2)-1];
                else if (dum2<=dum1 && dum2<dum3)
                     p(count,:) = [curindex(1)-1 curindex(2)];
                     count = count+1;
                     p(count,:) = [curindex(1)-1 curindex(2)-1];
                     count = count+1;
                     p(count,:) = [curindex(1)-1 curindex(2)-2];
                    else
                        p(count,:) = [curindex(1) curindex(2)-1];
                        count = count+1;
                        p(count,:) = [curindex(1)-1 curindex(2)-1];
                        count = count+1;
                        p(count,:) = [curindex(1)-2 curindex(2)-1];
                    end
                end
                
               % [dum1 dum2] = min([D(curindex(1)-1,curindex(2)-1) D(curindex(1)-1,curindex(2)) D(curindex(1),curindex(2)-1)]);
                % p(count,:) = (dum2==1)*[curindex(1)-1 curindex(2)-1] + (dum2==2)*[curindex(1)-1 curindex(2)] + (dum2==3)*[curindex(1) curindex(2)-1] ;
       
             end
        end
    end
        
    
    
    
    count = count+1;
    
    
    
end



DTWdist = D(M1,b);


p=flipud(p);



end

