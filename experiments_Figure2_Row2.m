anz_runs=10;
n_range=1000:1000:10000;
n_range_fair_norm=1000:1000:4000;



%%%% FIRST PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k=3;
h=5;

a=0.4;
b=0.3;
c=0.2;
d=0.1;



error_SC=zeros(anz_runs,length(n_range));
error_SC_Normalized=zeros(anz_runs,length(n_range));
error_Fair_SC=zeros(anz_runs,length(n_range));
error_Fair_SC_Normalized=zeros(anz_runs,length(n_range_fair_norm));

time_SC=zeros(anz_runs,length(n_range));
time_SC_Normalized=zeros(anz_runs,length(n_range));
time_Fair_SC=zeros(anz_runs,length(n_range));
time_Fair_SC_Normalized=zeros(anz_runs,length(n_range_fair_norm));

for mmm=1:length(n_range)
    
    n=n_range(mmm);
        
    disp('-------------------------------------------------')
    n
    
    block_sizes=[7*(n/10)*[0.2,0.2,0.2,0.2,0.2],15*(n/100)*[0.2,0.2,0.2,0.2,0.2],15*(n/100)*[0.2,0.2,0.2,0.2,0.2]];
    
    sensitive=zeros(n,1);
    labels=zeros(n,1);
    for yyy=1:k
        for zzz=1:h
            block=(yyy-1)*h+zzz;
            sensitive((sum(block_sizes(1:(block-1)))+1):sum(block_sizes(1:block)))=zzz;
            labels((sum(block_sizes(1:(block-1)))+1):sum(block_sizes(1:block)))=yyy;
        end
    end
        
    
    
    for ell=1:anz_runs
        
        adja=generate_adja_SB_model(n,a,b,c,d,k,h,block_sizes);
        
        tstart=tic;
        clustering_SC=SC_unnormalized(adja,k);
        time_SC(ell,mmm)=toc(tstart);
        error_SC(ell,mmm)=clustering_accuracy(labels,clustering_SC);
        
        
        tstart=tic;
        clustering_SC_NORMALIZED=SC_normalized(adja,k);
        time_SC_Normalized(ell,mmm)=toc(tstart);
        error_SC_Normalized(ell,mmm)=clustering_accuracy(labels,clustering_SC_NORMALIZED);
        
        
        tstart=tic;
        FAIR_clustering_SC=Fair_SC_unnormalized(adja,k,sensitive);
        time_Fair_SC(ell,mmm)=toc(tstart);
        error_Fair_SC(ell,mmm)=clustering_accuracy(labels, FAIR_clustering_SC);
        
        
        if mmm<=length(n_range_fair_norm)
            tstart=tic;
            FAIR_clustering_SC_NORMALIZED=Fair_SC_normalized(adja,k,sensitive);
            time_Fair_SC_Normalized(ell,mmm)=toc(tstart);
            error_Fair_SC_Normalized(ell,mmm)=clustering_accuracy(labels,FAIR_clustering_SC_NORMALIZED);
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
    

sfname=strcat('_SB_model_as_function_of_n_with_k=',num2str(k),'_h=',num2str(h),'_runs=',num2str(anz_runs),'_Row2Plot1');    

%save(strcat('DATA',sfname,'.mat'))

figure(1);clf;
plot(n_range,mean(error_SC,1),'mo-.',n_range,mean(error_SC_Normalized,1),'gs-.',n_range,mean(error_Fair_SC,1),'rx-',n_range_fair_norm,mean(error_Fair_SC_Normalized,1),'bx-')
legend('SC (Alg. 1)','Normalized SC','FAIR SC (Alg. 2)','FAIR Norm. SC (Alg. 3)')
xlabel('n')
ylabel('Error')
ylim([0,1])
title(strcat('k=',num2str(k),', h=',num2str(h),' --- a=',num2str(a),', b=',num2str(b),', c=',num2str(c),', d=',num2str(d)),'FontWeight','normal')
saveas(1,strcat('Error',sfname))
print(1,'-dpdf',strcat('Error',sfname))

figure(2);clf;
plot(n_range,mean(time_SC,1),'mo-.',n_range,mean(time_SC_Normalized,1),'gs-.',n_range,mean(time_Fair_SC,1),'rx-',n_range_fair_norm,mean(time_Fair_SC_Normalized,1),'bx-',...
    n_range_fair_norm,1.2*(time_Fair_SC_Normalized(1)/(n_range_fair_norm(1)^3))*n_range_fair_norm.^3,'c-')
legend('SC (Alg. 1)','Normalized SC','FAIR SC (Alg. 2)','FAIR Norm. SC (Alg. 3)','~ n^3')
xlabel('n')
ylabel('Running time [s]')
title(strcat('k=',num2str(k),', h=',num2str(h),' --- a=',num2str(a),', b=',num2str(b),', c=',num2str(c),', d=',num2str(d)),'FontWeight','normal')
saveas(2,strcat('Time',sfname))
print(2,'-dpdf',strcat('Time',sfname))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









%%%% SECOND PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k=2;
h=2;

af=10;
bf=7;
cf=4;
df=1;



error_SC=zeros(anz_runs,length(n_range));
error_SC_Normalized=zeros(anz_runs,length(n_range));
error_Fair_SC=zeros(anz_runs,length(n_range));
error_Fair_SC_Normalized=zeros(anz_runs,length(n_range_fair_norm));

time_SC=zeros(anz_runs,length(n_range));
time_SC_Normalized=zeros(anz_runs,length(n_range));
time_Fair_SC=zeros(anz_runs,length(n_range));
time_Fair_SC_Normalized=zeros(anz_runs,length(n_range_fair_norm));

for mmm=1:length(n_range)
    
    n=n_range(mmm);
    
    disp('-------------------------------------------------')
    n    

    a=af*(log(n)/n)^(2/3);
    b=bf*(log(n)/n)^(2/3);
    c=cf*(log(n)/n)^(2/3);
    d=df*(log(n)/n)^(2/3);
    
    expo1=2;
    expo2=3;
    
    
    block_sizes=[6*(n/10)*[0.6,0.4],4*(n/10)*[0.6,0.4]];
    
    sensitive=zeros(n,1);
    labels=zeros(n,1);
    for yyy=1:k
        for zzz=1:h
            block=(yyy-1)*h+zzz;
            sensitive((sum(block_sizes(1:(block-1)))+1):sum(block_sizes(1:block)))=zzz;
            labels((sum(block_sizes(1:(block-1)))+1):sum(block_sizes(1:block)))=yyy;
        end
    end
        
    
    
    for ell=1:anz_runs
    
        adja=generate_adja_SB_model(n,a,b,c,d,k,h,block_sizes);
        
        tstart=tic;
        clustering_SC=SC_unnormalized(adja,k);
        time_SC(ell,mmm)=toc(tstart);
        error_SC(ell,mmm)=clustering_accuracy(labels,clustering_SC);
        
        
        tstart=tic;
        clustering_SC_NORMALIZED=SC_normalized(adja,k);
        time_SC_Normalized(ell,mmm)=toc(tstart);
        error_SC_Normalized(ell,mmm)=clustering_accuracy(labels,clustering_SC_NORMALIZED);
        
        
        tstart=tic;
        FAIR_clustering_SC=Fair_SC_unnormalized(adja,k,sensitive);
        time_Fair_SC(ell,mmm)=toc(tstart);
        error_Fair_SC(ell,mmm)=clustering_accuracy(labels, FAIR_clustering_SC);
        
        
        if mmm<=length(n_range_fair_norm)
            tstart=tic;
            FAIR_clustering_SC_NORMALIZED=Fair_SC_normalized(adja,k,sensitive);
            time_Fair_SC_Normalized(ell,mmm)=toc(tstart);
            error_Fair_SC_Normalized(ell,mmm)=clustering_accuracy(labels,FAIR_clustering_SC_NORMALIZED);
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
    

sfname=strcat('_SB_model_as_function_of_n_with_k=',num2str(k),'_h=',num2str(h),'_runs=',num2str(anz_runs),'_Row2Plot2');    

%save(strcat('DATA',sfname,'.mat'))

figure(1);clf;
plot(n_range,mean(error_SC,1),'mo-.',n_range,mean(error_SC_Normalized,1),'gs-.',n_range,mean(error_Fair_SC,1),'rx-',n_range_fair_norm,mean(error_Fair_SC_Normalized,1),'bx-')
legend('SC (Alg. 1)','Normalized SC','FAIR SC (Alg. 2)','FAIR Norm. SC (Alg. 3)')
xlabel('n')
ylabel('Error')
ylim([0,1])
title(strcat('k=',num2str(k),', h=',num2str(h),' --- a,b,c,d ~ (log(n)/n)\^(',num2str(expo1),'/',num2str(expo2),')'),'FontWeight','normal')
saveas(1,strcat('Error',sfname))
print(1,'-dpdf',strcat('Error',sfname))

figure(2);clf;
plot(n_range,mean(time_SC,1),'mo-.',n_range,mean(time_SC_Normalized,1),'gs-.',n_range,mean(time_Fair_SC,1),'rx-',n_range_fair_norm,mean(time_Fair_SC_Normalized,1),'bx-',...
    n_range_fair_norm,1.2*(time_Fair_SC_Normalized(1)/(n_range_fair_norm(1)^3))*n_range_fair_norm.^3,'c-')
legend('SC (Alg. 1)','Normalized SC','FAIR SC (Alg. 2)','FAIR Norm. SC (Alg. 3)','~ n^3')
xlabel('n')
ylabel('Running time [s]')
title(strcat('k=',num2str(k),', h=',num2str(h),' --- a,b,c,d ~ (log(n)/n)\^(',num2str(expo1),'/',num2str(expo2),')'),'FontWeight','normal')
saveas(2,strcat('Time',sfname))
print(2,'-dpdf',strcat('Time',sfname))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%%% THIRD PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k=5;
h=3;

a=0.4;
b=0.3;
c=0.2;
d=0.1;



error_SC=zeros(anz_runs,length(n_range));
error_SC_Normalized=zeros(anz_runs,length(n_range));
error_Fair_SC=zeros(anz_runs,length(n_range));
error_Fair_SC_Normalized=zeros(anz_runs,length(n_range_fair_norm));

time_SC=zeros(anz_runs,length(n_range));
time_SC_Normalized=zeros(anz_runs,length(n_range));
time_Fair_SC=zeros(anz_runs,length(n_range));
time_Fair_SC_Normalized=zeros(anz_runs,length(n_range_fair_norm));

for mmm=1:length(n_range)
    
    n=n_range(mmm);
        
    disp('-------------------------------------------------')
    n
    
    block_sizes=[3*(n/10)*[0.5,0.3,0.2],3*(n/10)*[0.5,0.3,0.2],2*(n/10)*[0.5,0.3,0.2],1*(n/10)*[0.5,0.3,0.2],1*(n/10)*[0.5,0.3,0.2]];
    
    sensitive=zeros(n,1);
    labels=zeros(n,1);
    for yyy=1:k
        for zzz=1:h
            block=(yyy-1)*h+zzz;
            sensitive((sum(block_sizes(1:(block-1)))+1):sum(block_sizes(1:block)))=zzz;
            labels((sum(block_sizes(1:(block-1)))+1):sum(block_sizes(1:block)))=yyy;
        end
    end
        
    
    
    for ell=1:anz_runs
        
        adja=generate_adja_SB_model(n,a,b,c,d,k,h,block_sizes);
        
        tstart=tic;
        clustering_SC=SC_unnormalized(adja,k);
        time_SC(ell,mmm)=toc(tstart);
        error_SC(ell,mmm)=clustering_accuracy(labels,clustering_SC);
        
        
        tstart=tic;
        clustering_SC_NORMALIZED=SC_normalized(adja,k);
        time_SC_Normalized(ell,mmm)=toc(tstart);
        error_SC_Normalized(ell,mmm)=clustering_accuracy(labels,clustering_SC_NORMALIZED);
        
        
        tstart=tic;
        FAIR_clustering_SC=Fair_SC_unnormalized(adja,k,sensitive);
        time_Fair_SC(ell,mmm)=toc(tstart);
        error_Fair_SC(ell,mmm)=clustering_accuracy(labels, FAIR_clustering_SC);
        
        
        if mmm<=length(n_range_fair_norm)
            tstart=tic;
            FAIR_clustering_SC_NORMALIZED=Fair_SC_normalized(adja,k,sensitive);
            time_Fair_SC_Normalized(ell,mmm)=toc(tstart);
            error_Fair_SC_Normalized(ell,mmm)=clustering_accuracy(labels,FAIR_clustering_SC_NORMALIZED);
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
    

sfname=strcat('_SB_model_as_function_of_n_with_k=',num2str(k),'_h=',num2str(h),'_runs=',num2str(anz_runs),'_Row2Plot3');    

%save(strcat('DATA',sfname,'.mat'))

figure(1);clf;
plot(n_range,mean(error_SC,1),'mo-.',n_range,mean(error_SC_Normalized,1),'gs-.',n_range,mean(error_Fair_SC,1),'rx-',n_range_fair_norm,mean(error_Fair_SC_Normalized,1),'bx-')
legend('SC (Alg. 1)','Normalized SC','FAIR SC (Alg. 2)','FAIR Norm. SC (Alg. 3)')
xlabel('n')
ylabel('Error')
ylim([0,1])
title(strcat('k=',num2str(k),', h=',num2str(h),' --- a=',num2str(a),', b=',num2str(b),', c=',num2str(c),', d=',num2str(d)),'FontWeight','normal')
saveas(1,strcat('Error',sfname))
print(1,'-dpdf',strcat('Error',sfname))

figure(2);clf;
plot(n_range,mean(time_SC,1),'mo-.',n_range,mean(time_SC_Normalized,1),'gs-.',n_range,mean(time_Fair_SC,1),'rx-',n_range_fair_norm,mean(time_Fair_SC_Normalized,1),'bx-',...
    n_range_fair_norm,1.2*(time_Fair_SC_Normalized(1)/(n_range_fair_norm(1)^3))*n_range_fair_norm.^3,'c-')
legend('SC (Alg. 1)','Normalized SC','FAIR SC (Alg. 2)','FAIR Norm. SC (Alg. 3)','~ n^3')
xlabel('n')
ylabel('Running time [s]')
title(strcat('k=',num2str(k),', h=',num2str(h),' --- a=',num2str(a),', b=',num2str(b),', c=',num2str(c),', d=',num2str(d)),'FontWeight','normal')
saveas(2,strcat('Time',sfname))
print(2,'-dpdf',strcat('Time',sfname))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%%% FOURTH PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k=5;
h=3;

a=0.4;
b=0.3;
c=0.2;
d=0.1;



error_SC=zeros(anz_runs,length(n_range));
error_SC_Normalized=zeros(anz_runs,length(n_range));
error_Fair_SC=zeros(anz_runs,length(n_range));
error_Fair_SC_Normalized=zeros(anz_runs,length(n_range_fair_norm));

time_SC=zeros(anz_runs,length(n_range));
time_SC_Normalized=zeros(anz_runs,length(n_range));
time_Fair_SC=zeros(anz_runs,length(n_range));
time_Fair_SC_Normalized=zeros(anz_runs,length(n_range_fair_norm));

for mmm=1:length(n_range)
    
    n=n_range(mmm);
        
    disp('-------------------------------------------------')
    n
    
    block_sizes=[6*(n/10)*[0.5,0.3,0.2],1*(n/10)*[0.5,0.3,0.2],1*(n/10)*[0.5,0.3,0.2],1*(n/10)*[0.5,0.3,0.2],1*(n/10)*[0.5,0.3,0.2]];
    
    sensitive=zeros(n,1);
    labels=zeros(n,1);
    for yyy=1:k
        for zzz=1:h
            block=(yyy-1)*h+zzz;
            sensitive((sum(block_sizes(1:(block-1)))+1):sum(block_sizes(1:block)))=zzz;
            labels((sum(block_sizes(1:(block-1)))+1):sum(block_sizes(1:block)))=yyy;
        end
    end
        
    
    
    for ell=1:anz_runs
         
        adja=generate_adja_SB_model(n,a,b,c,d,k,h,block_sizes);
        
        tstart=tic;
        clustering_SC=SC_unnormalized(adja,k);
        time_SC(ell,mmm)=toc(tstart);
        error_SC(ell,mmm)=clustering_accuracy(labels,clustering_SC);
        
        
        tstart=tic;
        clustering_SC_NORMALIZED=SC_normalized(adja,k);
        time_SC_Normalized(ell,mmm)=toc(tstart);
        error_SC_Normalized(ell,mmm)=clustering_accuracy(labels,clustering_SC_NORMALIZED);
        
        
        tstart=tic;
        FAIR_clustering_SC=Fair_SC_unnormalized(adja,k,sensitive);
        time_Fair_SC(ell,mmm)=toc(tstart);
        error_Fair_SC(ell,mmm)=clustering_accuracy(labels, FAIR_clustering_SC);
        
        
        if mmm<=length(n_range_fair_norm)
            tstart=tic;
            FAIR_clustering_SC_NORMALIZED=Fair_SC_normalized(adja,k,sensitive);
            time_Fair_SC_Normalized(ell,mmm)=toc(tstart);
            error_Fair_SC_Normalized(ell,mmm)=clustering_accuracy(labels,FAIR_clustering_SC_NORMALIZED);
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
    

sfname=strcat('_SB_model_as_function_of_n_with_k=',num2str(k),'_h=',num2str(h),'_runs=',num2str(anz_runs),'_Row2Plot4');    

%save(strcat('DATA',sfname,'.mat'))

figure(1);clf;
plot(n_range,mean(error_SC,1),'mo-.',n_range,mean(error_SC_Normalized,1),'gs-.',n_range,mean(error_Fair_SC,1),'rx-',n_range_fair_norm,mean(error_Fair_SC_Normalized,1),'bx-')
legend('SC (Alg. 1)','Normalized SC','FAIR SC (Alg. 2)','FAIR Norm. SC (Alg. 3)')
xlabel('n')
ylabel('Error')
ylim([0,1])
title(strcat('k=',num2str(k),', h=',num2str(h),' --- a=',num2str(a),', b=',num2str(b),', c=',num2str(c),', d=',num2str(d)),'FontWeight','normal')
saveas(1,strcat('Error',sfname))
print(1,'-dpdf',strcat('Error',sfname))

figure(2);clf;
plot(n_range,mean(time_SC,1),'mo-.',n_range,mean(time_SC_Normalized,1),'gs-.',n_range,mean(time_Fair_SC,1),'rx-',n_range_fair_norm,mean(time_Fair_SC_Normalized,1),'bx-',...
    n_range_fair_norm,1.2*(time_Fair_SC_Normalized(1)/(n_range_fair_norm(1)^3))*n_range_fair_norm.^3,'c-')
legend('SC (Alg. 1)','Normalized SC','FAIR SC (Alg. 2)','FAIR Norm. SC (Alg. 3)','~ n^3')
xlabel('n')
ylabel('Running time [s]')
title(strcat('k=',num2str(k),', h=',num2str(h),' --- a=',num2str(a),', b=',num2str(b),', c=',num2str(c),', d=',num2str(d)),'FontWeight','normal')
saveas(2,strcat('Time',sfname))
print(2,'-dpdf',strcat('Time',sfname))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

