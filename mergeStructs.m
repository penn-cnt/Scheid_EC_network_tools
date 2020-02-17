function tostruct=mergeStructs(fromstruct, tostruct)
for f= fieldnames(fromstruct)'
    tostruct.(f{1})=fromstruct.(f{1});
end
end