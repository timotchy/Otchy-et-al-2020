function fieldVector=getAnnotationVector(handles,field)

for k = 1:size(handles,2)
    AnnStruct=handles{1,k};
    if ~isstr(getfield(AnnStruct, field));
        fieldVector(k) = getfield(AnnStruct, field);
    else
        fieldVector{k} = getfield(AnnStruct, field);
    end
end