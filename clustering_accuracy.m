function error=clustering_accuracy(labels,clustering)
%computes the error of a clustering with respect to the ground-truth
%clustering (see Section 4)
%
%INPUT:
%labels ... vector of length n comprising ground-truth cluster labels
%clustering ... vector of length n comprising cluster labels provided by a
%               clustering algorithm
%
%OUTPUT:
%error ... error of the clustering


n=length(labels);


if sum(size(labels)==[n,1])==2
    labels=reshape(labels,[1,n]);
end

if sum(size(clustering)==[n,1])==2
    clustering=reshape(clustering,[1,n]);
end


aa=unique(labels);
J=length(aa);

bb=unique(clustering);
K=length(bb);


if sum(aa==(1:J))<J
    labels_old=labels;
    temp=1;
    for ell=aa
        labels(labels_old==ell)=temp;
        temp=temp+1;
    end
end

if sum(bb==(1:K))<K
    clustering_old=clustering;
    temp=1;
    for ell=bb
        clustering(clustering_old==ell)=temp;
        temp=temp+1;
    end
end



permut=perms(1:max(K,J));
Kfac=size(permut,1);

error=Inf;
clustering_temp=clustering;

for ell=1:Kfac
    for mmm=1:K
        clustering_temp(clustering==mmm)=permut(ell,mmm);
    end    
    error_temp=sum(clustering_temp~=labels)/n;
    if error_temp<error
        error=error_temp;
    end
end
    
    