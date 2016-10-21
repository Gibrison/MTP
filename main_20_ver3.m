%% creating feature matrix : point based features
%========================================================================

clc;
clear all;
close all;
tic
load refine_20
M=refine_20;
[ro,co]=size(M);
M(:,7)=M(:,7)-30;

%% normalization      
%=====================================================================
% mean-variance method
% %mx=sum(M(:,1))/ro;
% %sig_x=sqrt(sum((M(:,1)-mx).*(M(:,1))/ro));

 M(:,1)=(M(:,1)-sum(M(:,1))/ro)/sqrt(sum((M(:,1)-sum(M(:,1))/ro).*(M(:,1))/ro));
% 
% %my=sum(M(:,2))/ro;
% %sig_y=sqrt(sum((M(:,2)-sum(M(:,2))/ro).*(M(:,2))/ro));

 M(:,2)=(M(:,2)-sum(M(:,2))/ro)/sqrt(sum((M(:,2)-sum(M(:,2))/ro).*(M(:,2))/ro));
 
%===================================================================
% % min-max method
% 
% M(:,1)=M(:,1)/(max(M(:,1))-min(M(:,1)));
% 
% M(:,2)=M(:,2)/(max(M(:,2))-min(M(:,2)));
%====================================================================
c=0;s1=0;
for i=1:20
    len1=find(M(:,7)==i);
   % len1(end)
    temp=M(len1(1):len1(end),6);
    len2=unique(temp(:,1));
   
    for j=1:len2(end)
        len3=find(temp(:,1)==j);
        %len3(end);
        s2=length(len3);
        s1=s1+s2;
        c=c+1;
        ends(c,1)=s1;
        
        t_gap(s1-s2+1:s1-1,1)=(M(s1-s2+2:s1,3)-M(s1-s2+1:s1-1,3));
        
        dist(s1-s2+1:s1-1,1)=sqrt((M(s1-s2+2:s1,1)-M(s1-s2+1:s1-1,1)).^2+(M(s1-s2+2:s1,2)-M(s1-s2+1:s1-1,2)).^2);
        cos_t(s1-s2+1:s1-1,1)=(M(s1-s2+2:s1,1)-M(s1-s2+1:s1-1,1));
        sin_t(s1-s2+1:s1-1,1)=(M(s1-s2+2:s1,2)-M(s1-s2+1:s1-1,2));
        %cos_p(s1-s2+1:s1-2,1)=cos_t(s1-s2+1:s1-2,1).*cos_t(s1-s2+2:s1-1,1)+sin_t(s1-s2+1:s1-2,1).*sin_t(s1-s2+2:s1-1,1);
        %sin_p(s1-s2+1:s1-2,1)=cos_t(s1-s2+1:s1-2,1).*sin_t(s1-s2+2:s1-1,1)-sin_t(s1-s2+1:s1-2,1).*cos_t(s1-s2+2:s1-1,1);
        % clear len3
        %label(s1-s2+1:s1-1,1)=M(s1-s2+1:s1-1,7);    % confirm that...assignment is corect 
        %
      
        
    end
    
    % len1
    % temp
     
end

 
 
 t_gap(ro,1)=0;
 dist(ro,1)=0;
 cos_t(ro,1)=0;
 sin_t(ro,1)=0;
 
 %zero_t=find(t_gap==0);
 
 % solving the issue of zero distance 
 
% zero_dist=find(dist==0);
 cos_t=cos_t./dist;
 sin_t=sin_t./dist;
 cos_t(find(dist==0))=1;
 sin_t(find(dist==0))=0;
 
 cos_p(:,1)=cos_t(1:ro-2,1).*cos_t(2:ro-1,1)+sin_t(1:ro-2,1).*sin_t(2:ro-1,1);
 sin_p(:,1)=cos_t(1:ro-2,1).*sin_t(2:ro-1,1)-sin_t(1:ro-2,1).*cos_t(2:ro-1,1);
 cos_p(ro-1:ro,1)=0;
 sin_p(ro-1:ro,1)=0;
 
%  cos_p(ends)=[];
%  sin_p(ends)=[];
  t_gap(ends)=10^7;
%  dist(ends)=[];
%  cos_t(ends)=[];
%  sin_t(ends)=[];
% %  
pt_ftr20(:,1) = dist./t_gap;      %% velocity feature
pt_ftr20(:,2) = cos_t;            %% cos theta angle
pt_ftr20(:,3) = sin_t;            %% sin theta angle
pt_ftr20(:,4) = cos_p;            %% cos phi angle
pt_ftr20(:,5) = sin_p;            %% sine phi angle
pt_ftr20(:,6) = M(1:ro,6);
pt_ftr20(:,7) = M(1:ro,7);
% pt_feature20(:,8) = M(1:ro-c,6);
% pt_feature20(:,9) = M(1:ro-c,7);
%pt_feature20(ro,:)=0;
%pt_feature20(ends,:)=[];
% %zz=find(M(11523:14104,6)==1);
pt_ftr20(ends,:)=[];
%=========================================================================
%% trainig and testing data generation from the faeture matrix (pt_ftr20)
%=========================================================================
%Tr=[];Ts=[];
TrainD_20=[];TestD_20=[];
%pt_ftr20=pt_ftr20;
%s3=0;s5=0;
for i=1:20
%     len4=find(D(:,7)==i);
%     
%     
%      
    l1=find(pt_ftr20(:,6)==1 & pt_ftr20(:,7)==i); l4=find(pt_ftr20(:,6)==4 & pt_ftr20(:,7)==i);
    l5=find(pt_ftr20(:,6)==5 & pt_ftr20(:,7)==i); l8=find(pt_ftr20(:,6)==8 & pt_ftr20(:,7)==i);
%     l1(1);
%     l4(end);
     TrainD_20=[TrainD_20;pt_ftr20(l1(1):l4(end),[1:5,7])];
     TestD_20=[TestD_20;pt_ftr20(l5(1):l8(end),:)];
end
% 
% TrainD_20(:,6)=[];
% TestD_20(:,6)=[];
% 
% 
%=========================================================================
%% Training from the writer data
%=========================================================================
tic
Tr=TrainD_20(:,1:5);

k=50;   % no of gaussian components


 opt = statset('MaxIter',1000);
% GMModel = gmdistribution.fit(Tr,k,'CovType','diagonal','Regularize',0.0002,'Options',opt); %'CovType','diagonal',
% 
% cov_m=GMModel.Sigma;
% u=GMModel.mu;
% w=GMModel.ComponentProportion;
%  %[idx,nlogl,P] = cluster(GMModel,TrainD_20(1:200,1:5));
% 
% 
% r=16;  % relevance factor
% 
%  pp=zeros(20000,k);
% disp('adaptation started....................');
% for i=1:20
%     
%      
%      l=find(TrainD_20(:,6)==i);
% %      [idx,nlogl,P] = cluster(GMModel,Tr(l(1):l(end),1:5));
%      c1=0;
%    for i1=l(1):l(end)
%      for m1=1:k
%          c1=c1+1;
%         pp(c1,m1)=mvnpdf(Tr(i1,:),u(m1,:),diag(cov_m(1,:,m1)))*w(m1);
%      end
%    end
%    gmm_writer= gmdistribution.fit(Tr(l(1):l(end),:),k,'CovType','diagonal','Regularize',0.0002,'Options',opt);
%    cov_c=gmm_writer.Sigma;
%    u_c=gmm_writer.mu;
%    w_c=gmm_writer.ComponentProportion;
%    for m=1:k
%        alpha=sum(pp(1:length(l),m))/(sum(pp(1:length(l),m))+r);
%       mu_c(m,:)=alpha*u_c(m,:)+(1-alpha)*u(m,:); 
%       cov_c2(1,:,m)=alpha*cov_c(1,:,m)+(1-alpha)*cov_m(1,:,m);
%       wt_c(m)=alpha*w_c(m)+(1-alpha)*w(m);
%    end
% 
%  
%     mean{i}=mu_c;      % total mean matrix for all writers 2D
%     covar{i}=cov_c2;   % total cov matrix for all writers 3D
%     wt{i}=wt_c;
%    writer= i
% end
%  toc
 for i=1:20
    l=find(TrainD_20(:,6)==i);
    obj= gmdistribution.fit(Tr(l,:),k,'CovType','diagonal','Regularize',0.0002,'Options',opt);
    covv{i}=obj.Sigma;
    muu{i}=obj.mu;
    wttt{i}=obj.ComponentProportion;
 end
 toc
%=========================================================================
 %% Testing from the remaining data
 tic
 Likelihood=[];
 pp2=zeros(20000,k);
for i=1:20
    %l=find(TestD_20(:,6)==i);
   
   for j=5:8
       
      l=find(TestD_20(:,6)==j & TestD_20(:,7)==i);
      
      for wr=1:20
          c1=0;
           mu1=muu{wr};  
            var=covv{wr};
            wtt=wttt{wr};
         for t=l(1):l(end) %length(find(TestD_20(:,6)==i))
           %norm2=0;
           
           for m=1:k
           
            c1=c1+1;
           % norm2=norm2+wt(wr,m)*mvnpdf(TestD_20(t,1:5),mean((wr-1)*5+1:wr*5,m)',diag(var(1,:,m)));
            pp2(c1,m)=wtt(m)*mvnpdf(TestD_20(t,1:5),mu1(m,:),diag(var(1,:,m)));
           end
           
          %norm1=norm1+log(norm2);
          
         end
         
        LL(j,2+wr)=sum(log(sum(pp2(1:length(l),:)')));
        
      end
      LL(j,1)=i;
      LL(j,2)=j;
   end
   Likelihood=[Likelihood;LL];
end

%Likelihood=Likelihood;
Likelihood( ~any(Likelihood,2), : ) = [] ;%rows
LLR=Likelihood(:,1:2);
ac=0;
for i=1:80
    [maxx,loc]=max(Likelihood(i,3:22));
    LLR(i,3)=loc;
    if(LLR(i,1)==LLR(i,3))
        ac=ac+1;
    end
end
%  [model] = trainmsvm(TrainD_20(:,1:5), TrainD_20(:,6)) ;
%  [labels, outputs] = predmsvm(model, TestD_20(:,1:5), TestD_20(:,7), ncpus);
%clear c co cos_t cos_p sin_t sin_p ends i j l1 l4 l5 l8 len1 len2 len3 m s1 s2 ro sig_i temp

toc









