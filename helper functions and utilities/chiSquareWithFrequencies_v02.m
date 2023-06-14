function [tbl,chi2stat,pval] = chiSquareWithFrequencies_v02(n1,N1,n2,N2)
 % calculates the chi square statistic with summarized frequencies
 % Thomas Elston
 % 12.Dec.2019
 
 
 % n1 = numerator for fraction 1
 % N1 = denomenator for fraction 1
 % n2 = numerator for fraction 2
 % N2 = denomenator for fraction 2

       
       x1 = [repmat('a',N1,1); repmat('b',N2,1)];
       x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
       [tbl,chi2stat,pval] = crosstab(x1,x2);
       
       return;