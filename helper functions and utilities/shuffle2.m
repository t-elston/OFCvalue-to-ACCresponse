function shuffledV=shuffle2(v)
% pre-allocate
shuffledV = NaN(size(v));

% shuffle each column
numColumns = size(v,2);

for c = 1:numColumns
    
    thiscoldata = v(:,c);

     shuffledV(:,c)=thiscoldata(randperm(numel(thiscoldata)));
     
end % of cyling through columns

end