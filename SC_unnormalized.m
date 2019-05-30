function clusterLabels = SC_unnormalized(adj,k)
%implementation of unnormalized SC as stated in Alg. 1 
%
%INPUT:
%adj ... (weighted) adjacency matrix of size n x n
%k ... number of clusters
%
%OUTPUT:
%clusterLabels ... vector of length n comprising the cluster label for each
%                  data point


n=size(adj,1);

degrees = sum(adj, 1);
D = diag(degrees);
L = D-adj;

try
    [H, eigValues] = eigs(L,k,'smallestabs','MaxIterations',500,'SubspaceDimension',min(n,max(2*k,25)));
catch   
    [H, eigValues] = eigs(L,k,'smallestreal','MaxIterations',1000,'SubspaceDimension',min(n,max(2*k,25)));
end

clusterLabels = kmeans(H,k,'Replicates',10);
end




