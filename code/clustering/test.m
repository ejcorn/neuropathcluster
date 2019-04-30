X = rand(100); X = 0.5*(X+X');
utr_mask = triu(X)>0; utr_mask(logical(eye(length(X)))) = 0;
mean_triu = mean(X(utr_mask));
X(logical(eye(length(X)))) = 0;
mean_wd = mean(sum(X,2)/99);