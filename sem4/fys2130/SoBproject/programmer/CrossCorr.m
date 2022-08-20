function CC = CrossCorr(f, g)
 % Give cross-correlation of two functions f and g. 
 CC = sum(abs(f.*g));
end