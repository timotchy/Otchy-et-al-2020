%
% This function implements DTW algorithm. Inputs are two feature matrices
% F1 and F2. F1 is and L1-by-M1 matrix, each column represents the feature 
% at a particular time step. F1 is the template.
%
% This is a classical DTW. Unbanded...
%


function [DTWdist p] = DTWsub(F1,F2)


[L1 M1] = size(F1);
[L2 M2] = size(F2);

if L1 ~= L2
    disp('Size mismatch');
    return;
end

% Accumulated Cost Matrix


D = zeros(M1,M2);

for i = 1:M1
    
    D(i,1) = 0;
    
    for j = 1:i
        D(i,1) = D(i,1)+norm(F1(:,j)-F2(:,1));%cost(F1,F2,j,1);
    end
end


for i = 1:M2   
    D(1,i) = norm(F1(:,1)-F2(:,i));
end

 
for i = 2:M1
    for j = 2:M2
        D(i,j) = min(D(i-1,j-1),min(D(i-1,j),D(i,j-1)))+norm(F1(:,i)-F2(:,j));
    end
end


% Implement a variable stepsize



% OPTIMAL WARPING PATH


[a b] = min(D(M1,:));

p(1,:) = [M1 b];

count = 2;

while 1
    
    curindex = p(count-1,:);
    
    
    if curindex(1) == 1
            break;
        else if curindex(2) == 1
                p(count,:) = [curindex(1)-1 curindex(2)];
            else
                
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
                
               % [dum1 dum2] = min([D(curindex(1)-1,curindex(2)-1) D(curindex(1)-1,curindex(2)) D(curindex(1),curindex(2)-1)]);
                % p(count,:) = (dum2==1)*[curindex(1)-1 curindex(2)-1] + (dum2==2)*[curindex(1)-1 curindex(2)] + (dum2==3)*[curindex(1) curindex(2)-1] ;
       
            end
    end
        
    
    
    
    count = count+1;
    
    
    
end



DTWdist = D(M1,b);


p=flipud(p);



end

