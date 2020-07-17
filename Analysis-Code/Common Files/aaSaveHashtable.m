function aaSaveHashtable(filename, hash)
if(~hash.isEmpty)
    keysOld = hash.keys;
    elementsOld = hash.elements;
    %sort the anotation vector
    filenums=getAnnotationVector(elementsOld,'filenum');
    [sortedFilenums,sortIndex] = sort(filenums);
    for i=1:length(filenums)
        elements{i}=elementsOld{sortIndex(i)};
        keys{i}=keysOld{sortIndex(i)};
    end
    save(filename, 'keys', 'elements', '-v7'); %Changed version for faster saving 12/30/14
else
    keys = {};
    elements = {};
    save(filename, 'keys', 'elements', '-v7'); %Changed version for faster saving 12/30/14
end