[JSD] = function JSD-metric(data, variable)
% calculate the Jensen-Shannon divergence metric across an emsemble of genes on a single variable (cell type, treatment condition, etc) where genes are organized in columns.
A1 = data(variable,:);
X1 = A1/max(A1);
% normalize all values in emsemble to maximum value (1.0). NOTE: this is not a formal probability.
% Alternative method: normalize by range (interval normalization) 
% 0.0 = minimum value, 1.0 = maximal value.
interval-whole = MAX(A)-MIN(A);
interval = MAX(A)-A(xi);
An = interval/interval-whole;
%%
X2 = log10(X1);
% for customized (base n) logarithm, use the following line: X2 = (log10(X1)/log10(n)).
c3 = -(sum(X1)*sum(X2));
c2 = 1/length(A1);
c1 = c3 * c2;
JSD = c1-c2.*c3;
% calculates JSD for a single ensemble of genes.




