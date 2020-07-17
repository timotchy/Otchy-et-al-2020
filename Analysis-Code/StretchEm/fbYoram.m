%%% Forward-Backward algorithm

% This function calculates marginal distributions for the anchor points


function marginals = fbYoram(rendition,template,anchorPoints,a,alpha,beta) % a global var, alpha local var, beta spectral noise


% Initialize the forward matriz

[dum nSequence] = size(rendition);
nAnchor = length(anchorPoints);

forward = zeros(nAnchor,nSequence);
backward = zeros(nAnchor,nSequence);


for i = 1:nSequence
    
   forward(1,i) = exp(-sum((rendition(:,i)-template(:,anchorPoints(1))).^2)/2/beta);
   backward(end,i) = 1;
   
end

normalization(1) = sum(forward(1,:));
forward(1,:) = forward(1,:)/sum(forward(1,:));


for i = 2:nAnchor

    deltaTau = anchorPoints(i)-anchorPoints(i-1);
    
    for j = 1:nSequence
    
        for k = max(1,round(j-deltaTau-a-12*sqrt(alpha))):1:min(max(1,j-deltaTau-a+12*sqrt(alpha)),j)
            
%             if i == 2
%                 -mean(sum((linearWarp(rendition(:,min(k,j):max(k,j)),deltaTau+1)-template(:,anchorPoints(i-1):anchorPoints(i))).^2))/2/beta
%             end
            
            forward(i,j) = forward(i,j) + forward(i-1,k)*exp(-(j-k-deltaTau-a)^2/2/alpha)*exp(-sum(sum((linearWarp(rendition(:,min(k,j):max(k,j)),deltaTau+1)-template(:,anchorPoints(i-1):anchorPoints(i))).^2))/2/beta);
        end
       % a = forward(i,j)
    end
    
    normalization(i) = sum(forward(i,:));
    forward(i,:) = forward(i,:)/sum(forward(i,:));
    
end

for i = 2:nAnchor

       bdeltaTau = anchorPoints(nAnchor-i+2)-anchorPoints(nAnchor-i+1);

       
    for j = 1:nSequence
    
        for k = max(j,round(j+deltaTau+a-12*sqrt(alpha))):1:min(nSequence,j+deltaTau+a+12*sqrt(alpha))
            
         backward(nAnchor-i+1,j) = backward(nAnchor-i+1,j) + backward(nAnchor-i+2,k)*exp(-(k-j-bdeltaTau-a)^2/2/alpha)*exp(-sum(sum((linearWarp(rendition(:,min(k,j):max(k,j)),bdeltaTau+1)-template(:,anchorPoints(nAnchor-i+1):anchorPoints(nAnchor-i+2))).^2))/2/beta);
        
        end
    end
    
    backward(nAnchor-i+1,:) = backward(nAnchor-i+1,:)/normalization(nAnchor-i+2);
    
end


marginals = forward.*backward;

end