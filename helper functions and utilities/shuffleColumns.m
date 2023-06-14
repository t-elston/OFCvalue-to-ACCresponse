function shuffledV=shuffleColumns(v)
% pre-allocate
shuffledV = NaN(size(v));

% shuffle each column
numColumns = size(v,2);

for c = 1:numColumns   
     thiscoldata = v(:,c);     
     ShuffledColData=thiscoldata(randperm(numel(thiscoldata)));
     shuffledV(:,c)=ShuffledColData;
end % of cyling through columns

end % of function