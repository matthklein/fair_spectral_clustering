function clusterLabels = SC_normalized(adj,k)
%implementation of normalized SC as described in Appendix A 
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
Dsqinv=diag(1./sqrt(degrees));

L = D-adj;

Msymm=Dsqinv*L*Dsqinv;
Msymm=(Msymm+Msymm')/2;


try
    [H, eigValues] = eigs(Msymm,k,'smallestabs','MaxIterations',500,'SubspaceDimension',min(n,max(2*k,25)));
catch
    [H, eigValues] = eigs(Msymm,k,'smallestreal','MaxIterations',1000,'SubspaceDimension',min(n,max(2*k,25)));
end

clusterLabels = kmeans(Dsqinv*H,k,'Replicates',10);
end




