function [inc, index] =  annotationContainsFile(annotation,filename)
% Function checked in the passed annotation file for the filename
% specified.  If it is not found, inc is set to 0 and rec and ind will be
% set to []; if the file is included, the inc is set to 1, rec will
% contain the complete annotation record for the selected filename and
% index will return the index within the larger annotation file.

%Set to empty (the dafault value)
inc = false;
index =[];

% for i = 1:size(annotation,2)
%     %AnnStruct=annotation{k};
%     fieldVector(k) = getfield(annotation{k}, 'filename');
% end
index = strfind(annotation.keys, filename);

if ~isempty(index)
    inc=true;
end


