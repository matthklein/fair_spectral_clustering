anz_runs=10;
h=4;
k_range=2:8;
factor_range=1:0.5:3;

af=0.25;
bf=0.2;
cf=0.15;
df=0.1;



%%%% LEFT PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_app=5000;

error_Fair_SC=zeros(anz_runs,length(k_range),length(factor_range));

for mmm=1:length(k_range)
    
    k=k_range(mmm);
    n=ceil(n_app/(k*h))*(k*h);
    block_sizes=(n/(k*h))*ones(1,k*h);
        
    sensitive=zeros(n,1);
    labels=zeros(n,1);
    for yyy=1:k
        for zzz=1:h
            sensitive(((n/k)*(yyy-1)+(n/(k*h))*(zzz-1)+1):((n/k)*(yyy-1)+(n/(k*h))*zzz))=zzz;
            labels(((n/k)*(yyy-1)+(n/(k*h))*(zzz-1)+1):((n/k)*(yyy-1)+(n/(k*h))*zzz))=yyy;
        end
    end
    
    
    
    for hg=1:length(factor_range)
        
        a=af*factor_range(hg);
        b=bf*factor_range(hg);
        c=cf*factor_range(hg);
        d=df*factor_range(hg);
    
        for ell=1:anz_runs

            adja=generate_adja_SB_model(n,a,b,c,d,k,h,block_sizes);

            FAIR_clustering_SC=Fair_SC_unnormalized(adja,k,sensitive);
            
            error_Fair_SC(ell,mmm,hg)=clustering_accuracy(labels, FAIR_clustering_SC);

        end    

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
    
sfname=strcat('_SB_model_as_function_of_k_with_n~',num2str(n_app),'_h=',num2str(h),'_runs=',num2str(anz_runs),'_LeftPlot');

%save(strcat('DATA',sfname,'.mat'))

figure(1);clf;
plot(k_range,mean(error_Fair_SC(:,:,1),1),'Marker','x','DisplayName',strcat('F=',num2str(factor_range(1))))
hold on
for hg=2:length(factor_range)
    plot(k_range,mean(error_Fair_SC(:,:,hg),1),'Marker','x','DisplayName',strcat('F=',num2str(factor_range(hg))))
end  
hold off
legend
xlabel('k')
ylabel('Error')
ylim([0,1])
title(strcat('FAIR SC (Alg. 2) --- n~',num2str(n_app),', h=',num2str(h)),'FontWeight','normal')
saveas(1,strcat('Error',sfname))
print(1,'-dpdf',strcat('Error',sfname))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







%%%% RIGHT PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
n_app=2000;

error_Fair_SC_Normalized=zeros(anz_runs,length(k_range));

for mmm=1:length(k_range)
    
    k=k_range(mmm);
    n=ceil(n_app/(k*h))*(k*h);
    block_sizes=(n/(k*h))*ones(1,k*h);
    
    sensitive=zeros(n,1);
    labels=zeros(n,1);
    for yyy=1:k
        for zzz=1:h
            sensitive(((n/k)*(yyy-1)+(n/(k*h))*(zzz-1)+1):((n/k)*(yyy-1)+(n/(k*h))*zzz))=zzz;
            labels(((n/k)*(yyy-1)+(n/(k*h))*(zzz-1)+1):((n/k)*(yyy-1)+(n/(k*h))*zzz))=yyy;
        end
    end
        
    
    
    for hg=1:length(factor_range)
        
        a=af*factor_range(hg);
        b=bf*factor_range(hg);
        c=cf*factor_range(hg);
        d=df*factor_range(hg);
    
        for ell=1:anz_runs

            adja=generate_adja_SB_model(n,a,b,c,d,k,h,block_sizes);

            FAIR_clustering_SC_NORMALIZED=Fair_SC_normalized(adja,k,sensitive);
            
            error_Fair_SC_Normalized(ell,mmm,hg)=clustering_accuracy(labels,FAIR_clustering_SC_NORMALIZED);

        end    

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

sfname=strcat('_SB_model_as_function_of_k_with_n~',num2str(n_app),'_h=',num2str(h),'_runs=',num2str(anz_runs),'_RightPlot');

%save(strcat('DATA',sfname,'.mat'))

figure(1);clf;
plot(k_range,mean(error_Fair_SC_Normalized(:,:,1),1),'Marker','x','DisplayName',strcat('F=',num2str(factor_range(1))))
hold on
for hg=2:length(factor_range)
    plot(k_range,mean(error_Fair_SC_Normalized(:,:,hg),1),'Marker','x','DisplayName',strcat('F=',num2str(factor_range(hg))))
end  
hold off
legend
xlabel('k')
ylabel('Error')
ylim([0,1])
title(strcat('FAIR Norm. SC (Alg. 3) --- n~',num2str(n_app),', h=',num2str(h)),'FontWeight','normal')
saveas(1,strcat('Error',sfname))
print(1,'-dpdf',strcat('Error',sfname))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    