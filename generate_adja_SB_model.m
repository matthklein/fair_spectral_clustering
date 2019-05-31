function A=generate_adja_SB_model(n,a,b,c,d,k,h,block_sizes)
%generates
%
%INPUT:
%n ... number of elements
%a,b,c,d ... parameters / probabilities
%k ... number of clusters
%h ... number of groups
%block_sizes ... vector of length k*h with sum(block_sizes)=n; the 
%                first h entries correspond to the first cluster for 
%                the h groups, and so on ... 
%
%OUTPUT
%A ... adjacency matrix of size n x n


if (sum(block_sizes)~=n)||(length(block_sizes)~=(k*h))
    error('wrong input')
end

adja=random('Binomial',1,d,n,n);

for ell=1:k
    for mmm=1:k
        for ggg=1:h
            for fff=1:h
                if ell==mmm
                    if ggg==fff
                        adja((sum(block_sizes(1:((ell-1)*h+(ggg-1))))+1):(sum(block_sizes(1:((ell-1)*h+(ggg))))),...
                            (sum(block_sizes(1:((ell-1)*h+(ggg-1))))+1):(sum(block_sizes(1:((ell-1)*h+(ggg))))))=random('Binomial',1,a,block_sizes((ell-1)*h+(ggg)),block_sizes((ell-1)*h+(ggg)));
                    else
                        adja((sum(block_sizes(1:((ell-1)*h+(ggg-1))))+1):(sum(block_sizes(1:((ell-1)*h+(ggg))))),...
                            (sum(block_sizes(1:((ell-1)*h+(fff-1))))+1):(sum(block_sizes(1:((ell-1)*h+(fff))))))=random('Binomial',1,c,block_sizes((ell-1)*h+(ggg)),block_sizes((ell-1)*h+(fff)));
                    end
                    
                else
                    if ggg==fff
                        adja((sum(block_sizes(1:((ell-1)*h+(ggg-1))))+1):(sum(block_sizes(1:((ell-1)*h+(ggg))))),...
                                (sum(block_sizes(1:((mmm-1)*h+(ggg-1))))+1):(sum(block_sizes(1:((mmm-1)*h+(ggg))))))=random('Binomial',1,b,block_sizes((ell-1)*h+(ggg)),block_sizes((mmm-1)*h+(ggg)));
                    end
                end
            end
        end
    end
end
            



A=triu(adja,1);
A=A+A';