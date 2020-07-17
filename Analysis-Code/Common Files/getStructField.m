function out = getStructField(structure, field, pad)
%Set default
if nargin < 3
    pad = 0;
end

%Determine the size of the structure.
numSteps = length(structure);

%If padding on, run through structure to find max length
L = [];
if pad
    for i = 1:numSteps
        %Get the field
        s = getfield(structure(i), field);
        L(i) = length(s);
    end
    maxL = max(L);
end

%Strip out desired data in a data-type appropriate structure
out = [];
for i = 1:numSteps
    %Get the field
    s = getfield(structure(i), field);
    if isstr(s) || iscell(s)
 %if isstr(getfield(structure(1), field))
        %Output to a cell array
        out{i} = s;
    elseif isnumeric(getfield(structure(1), field))
        d = size(s);
        if ismember(1, d) && length(d)==1
            %Output as a numeric value
            out(i) = s;
        elseif (ismember(1, d) && length(d)==2) || (~ismember(1, d) && length(d)==1)
            %Added padding correction 12/14/15
            if pad
                s = padarray(s(:), maxL-L(i), NaN, 'post'); 
            end
            
            %Output as a numeric array
            out(i,:) = s;
        elseif ismember(1, d) && length(d)==3 || (~ismember(1, d) && length(d)==2)
            %Added padding correction 11/08/16
            if pad
                ss = [];
                for j = 1:size(s,2)
                    ss(:,j) = padarray(s(:,j), maxL-L(i), NaN, 'post');
                end
                s = ss;
            end
            
            %Output as a numeric matrix
            out(i,:, :) = s;
        elseif ismember(1, d) && length(d)==4 || (~ismember(1, d) && length(d)==3)
            %Output as a numeric matrix
            out(i,:, :, :) = s;
        end
        
    elseif islogical(getfield(structure(1), field))
        %Output to a logical array
        out(i) = s;
    elseif isstruct(getfield(structure(1), field))
        %Output to a structure array
        out = [out; s];
    end
end


