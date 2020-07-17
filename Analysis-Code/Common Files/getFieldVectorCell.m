function fieldVector=getFieldVectorCell(handles,field)

if isempty(getfield(handles, {1,1}, field))
    fieldVector = [];
else
    for k = 1:size(handles,2)
        fieldVector{k} = getfield(handles, {1,k}, field);
    end
 end