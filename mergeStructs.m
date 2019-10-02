function s2=mergeStructs(s1, s2)
for f= fieldnames(s1)'
    s2.(f{1})=s1.(f{1});
end
end