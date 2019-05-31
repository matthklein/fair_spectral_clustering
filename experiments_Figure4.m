anz_runs=10;
perturbation_range=0:0.05:1;
k=4;
h=2;



%%%% LEFT PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n=4000;

a=0.4;
b=0.3;
c=0.2;
d=0.1;


error_SC=zeros(anz_runs,length(perturbation_range));
error_SC_Normalized=zeros(anz_runs,length(perturbation_range));
error_Fair_SC=zeros(anz_runs,length(perturbation_range));
error_Fair_SC_Normalized=zeros(anz_runs,length(perturbation_range));

for mmm=1:length(perturbation_range)
    
    pp=perturbation_range(mmm);
    block_sizes=(n/(k*h))*ones(1,k*h);
    
        
    for ell=1:anz_runs
    
        sensitive=zeros(n,1);
        labels=zeros(n,1);
        for yyy=1:k
            for zzz=1:h
                sensitive(((n/k)*(yyy-1)+(n/(k*h))*(zzz-1)+1):((n/k)*(yyy-1)+(n/(k*h))*zzz))=zzz;
                labels(((n/k)*(yyy-1)+(n/(k*h))*(zzz-1)+1):((n/k)*(yyy-1)+(n/(k*h))*zzz))=yyy;
            end
        end
        
    
        temp1=random('binomial',1,pp,n,1);
        
        sensitive((temp1==1)&(sensitive==1))=2;
        
        adja=generate_adja_SB_model(n,a,b,c,d,k,h,block_sizes);
        
        clustering_SC=SC_unnormalized(adja,k);
        error_SC(ell,mmm)=clustering_accuracy(labels,clustering_SC);
        
        clustering_SC_NORMALIZED=SC_normalized(adja,k);
        error_SC_Normalized(ell,mmm)=clustering_accuracy(labels,clustering_SC_NORMALIZED);
        
        FAIR_clustering_SC=Fair_SC_unnormalized(adja,k,sensitive);
        error_Fair_SC(ell,mmm)=clustering_accuracy(labels, FAIR_clustering_SC);
        
        FAIR_clustering_SC_NORMALIZED=Fair_SC_normalized(adja,k,sensitive);
        error_Fair_SC_Normalized(ell,mmm)=clustering_accuracy(labels,FAIR_clustering_SC_NORMALIZED);
        
    end
end


%set default sizes for figures:
ulesfontsize = 24;
set(0, 'DefaultAxesFontSize', ulesfontsize);
set(0, 'DefaultTextFontSize', ulesfontsize);
set(0, 'DefaultUIControlFontSize', ulesfontsize);
set(0,'DefaultLineMarkerSize',ulesfontsize);
set(0,'DefaultLineLineWidth',1.5) 
set(gcf, 'PaperPositionMode','auto')
close all;

sfname=strcat('_SB_model_as_function_of_perturbation_with_n=',num2str(n),'_k=',num2str(k),'_h=',num2str(h),'_runs=',num2str(anz_runs),'_LeftPlot');

%save(strcat('DATA',sfname,'.mat'))

figure(1);clf;
plot(perturbation_range,mean(error_SC,1),'mo-.',perturbation_range,mean(error_SC_Normalized,1),'gs-.',perturbation_range,mean(error_Fair_SC,1),'rx-',perturbation_range,mean(error_Fair_SC_Normalized,1),'bx-')
legend('SC (Alg. 1)','Normalized SC','FAIR SC (Alg. 2)','FAIR Norm. SC (Alg. 3)')
xlabel('p')
ylabel('Error')
ylim([0,1])
title(strcat('n=',num2str(n),', k=',num2str(k),', h=',num2str(h),' --- a=',num2str(a),', b=',num2str(b),', c=',num2str(c),', d=',num2str(d)),'FontWeight','normal')
saveas(1,strcat('Error',sfname))
print(1,'-dpdf',strcat('Error',sfname))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%%% RIGHT PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n=2000;

a=0.5;
b=0.4;
c=0.3;
d=0.1;


error_SC=zeros(anz_runs,length(perturbation_range));
error_SC_Normalized=zeros(anz_runs,length(perturbation_range));
error_Fair_SC=zeros(anz_runs,length(perturbation_range));
error_Fair_SC_Normalized=zeros(anz_runs,length(perturbation_range));

for mmm=1:length(perturbation_range)
    
    pp=perturbation_range(mmm);
    block_sizes=(n/(k*h))*ones(1,k*h);
    
        
    for ell=1:anz_runs
    
        sensitive=zeros(n,1);
        labels=zeros(n,1);
        for yyy=1:k
            for zzz=1:h
                sensitive(((n/k)*(yyy-1)+(n/(k*h))*(zzz-1)+1):((n/k)*(yyy-1)+(n/(k*h))*zzz))=zzz;
                labels(((n/k)*(yyy-1)+(n/(k*h))*(zzz-1)+1):((n/k)*(yyy-1)+(n/(k*h))*zzz))=yyy;
            end
        end
        
    
        temp1=random('binomial',1,pp,n,1);
        
        sensitive((temp1==1)&(sensitive==1))=2;
        
        adja=generate_adja_SB_model(n,a,b,c,d,k,h,block_sizes);
        
        clustering_SC=SC_unnormalized(adja,k);
        error_SC(ell,mmm)=clustering_accuracy(labels,clustering_SC);
        
        clustering_SC_NORMALIZED=SC_normalized(adja,k);
        error_SC_Normalized(ell,mmm)=clustering_accuracy(labels,clustering_SC_NORMALIZED);
        
        FAIR_clustering_SC=Fair_SC_unnormalized(adja,k,sensitive);
        error_Fair_SC(ell,mmm)=clustering_accuracy(labels, FAIR_clustering_SC);
        
        FAIR_clustering_SC_NORMALIZED=Fair_SC_normalized(adja,k,sensitive);
        error_Fair_SC_Normalized(ell,mmm)=clustering_accuracy(labels,FAIR_clustering_SC_NORMALIZED);
        
    end
end


%set default sizes for figures:
ulesfontsize = 24;
set(0, 'DefaultAxesFontSize', ulesfontsize);
set(0, 'DefaultTextFontSize', ulesfontsize);
set(0, 'DefaultUIControlFontSize', ulesfontsize);
set(0,'DefaultLineMarkerSize',ulesfontsize);
set(0,'DefaultLineLineWidth',1.5) 
set(gcf, 'PaperPositionMode','auto')
close all;

sfname=strcat('_SB_model_as_function_of_perturbation_with_n=',num2str(n),'_k=',num2str(k),'_h=',num2str(h),'_runs=',num2str(anz_runs),'_RightPlot');

%save(strcat('DATA',sfname,'.mat'))

figure(1);clf;
plot(perturbation_range,mean(error_SC,1),'mo-.',perturbation_range,mean(error_SC_Normalized,1),'gs-.',perturbation_range,mean(error_Fair_SC,1),'rx-',perturbation_range,mean(error_Fair_SC_Normalized,1),'bx-')
legend('SC (Alg. 1)','Normalized SC','FAIR SC (Alg. 2)','FAIR Norm. SC (Alg. 3)')
xlabel('p')
ylabel('Error')
ylim([0,1])
title(strcat('n=',num2str(n),', k=',num2str(k),', h=',num2str(h),' --- a=',num2str(a),', b=',num2str(b),', c=',num2str(c),', d=',num2str(d)),'FontWeight','normal')
saveas(1,strcat('Error',sfname))
print(1,'-dpdf',strcat('Error',sfname))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

