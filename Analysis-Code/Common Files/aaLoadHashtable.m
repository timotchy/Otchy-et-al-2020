function hash = aaLoadHashtable(filename)
load(filename); %keys and elements
hash = mhashtable;
for(nKey = 1:length(keys))
    hash.put(keys{nKey}, elements{nKey});
end