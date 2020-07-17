function Vector = linearWarp(Vin,n)
%This function linearly warps columns of the the input matrix, Vin, to the length defined by n.
%
%Updated with comments 7/8/2015 TMO

%Determine the input matrix size
[k, l] = size(Vin);

%Predefine the output variable
Vector = zeros(k,n);

%If the warping is for something other than collapsing Vin to a single point, do this...
if n ~= 0 && n~= 1
    %Cycle through each bin of the output matrix
    for i = 1:n
        %Stretch/compress factor that determines how much of each interval in Vin will be in each interval of output
        %Compress: scFactor>1
        %Stretch: scFactor<1
        %Pass: scFactor==1
%         scFactor = (l-1)/(n-1);
        
        %Index for selecting bin of Vin to scale (pins the starts and ends together)
%         ind = (i-1)*scFactor+1;
        ind = (i-1)*(l-1)/(n-1)+1;
        
        %Extract fractional component of ind
        w = ind-floor(ind);
        
        %Warp linearly via weighted averaging
        Vector(:,i) = Vin(:,floor(ind))*(1-w) +Vin(:,ceil(ind))*w;
    end
    
else
    %...quickly collapse Vin to a single time point by averaging it all together across the whole Vin
    Vector = mean(Vin,2);
end


end