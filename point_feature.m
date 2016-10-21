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

%% normalization      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%=========================================================================
% mean-variance method
% mx=sum(M(:,1))/ro;
% sig_x=sqrt(sum((M(:,1)-mx).*(M(:,1))/ro));
% M(:,1)=(M(:,1)-mx)/sig_x;
% 
% my=sum(M(:,2))/ro;
% sig_y=sqrt(sum((M(:,2)-my).*(M(:,2))/ro));
% M(:,2)=(M(:,2)-my)/sig_y;
%===================================================================
% min-max method

M(:,1)=M(:,1)/(max(M(:,1))-min(M(:,1)));

M(:,2)=M(:,2)/(max(M(:,2))-min(M(:,2)));
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
 
 cos_p(ends)=[];
 sin_p(ends)=[];
 t_gap(ends)=[];
 dist(ends)=[];
 cos_t(ends)=[];
 sin_t(ends)=[];
 
pt_ftr20(:,1) = dist./t_gap;      %% velocity feature
pt_ftr20(:,2) = cos_t;            %% cos theta angle
pt_ftr20(:,3) = sin_t;            %% sin theta angle
pt_ftr20(:,4) = cos_p;            %% cos phi angle
pt_ftr20(:,5) = sin_p;            %% sine phi angle
pt_ftr20(:,6) = M(1:ro-c,6);
pt_ftr20(:,7) = M(1:ro-c,7);
% pt_feature20(:,8) = M(1:ro-c,6);
% pt_feature20(:,9) = M(1:ro-c,7);
%pt_feature20(ro,:)=0;
%pt_feature20(ends,:)=[];
% %zz=find(M(11523:14104,6)==1);

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
     TestD_20=[TestD_20;pt_ftr20(l5(1):l8(end),[1:5,7])];
end
% 
% TrainD_20(:,6)=[];
% TestD_20(:,6)=[];
% 

%=========================================================================
%% Training from the writer data
%=========================================================================

% Tr=TrainD_20(:,1:5);
% %Ts=TestD_20;
% k=50;   % no of gaussian components
% 
% %GMModel1 = fitgmdist(Tr,k,'CovarianceType','diagonal','RegularizationValue',0.0002);
% GMModel2 = gmdistribution.fit(Tr,k,'CovType','diagonal','Regularize',0.0002);
% 
% cov_m=GMModel2.Sigma;
% u=GMModel2.mu;
% w=GMModel2.ComponentProportion;
% 
% r=16;  % relevance factor
% mean=[];
% covar=[];
% wt=[];
% for i=1:20
%      mu_i=zeros(5,k);
%      sig_i=zeros(1,5,k);
%      norm2=0;
%      l=find(TrainD_20(:,6)==i);
%     for m=1:k
%     
%      norm1=0;
%      for t=l(1):l(end)      %length(find(TrainD_20(:,6)==i))
%         mu_i(:,m)=mu_i(:,m)+ (mvnpdf(Tr(t,:),u(m,:),diag(cov_m(:,:,m)))*Tr(t,:)*w(m))';      % column format
%         sig_i(1,:,m)=sig_i(1,:,m)+ (mvnpdf(Tr(t,:),u(m,:),diag(cov_m(:,:,m)))*Tr(t,:)*w(m));   % row format
%         norm1=norm1+mvnpdf(Tr(t,:),u(m,:),diag(cov_m(:,:,m)))*w(m);          % constant number
%         norm2=norm2+mvnpdf(Tr(t,:),u(m,:),diag(cov_m(:,:,m)))*w(m);
%         
%      end
%     % mu_i(:,m)=mu_i(:,m)/norm1;
%     % sig_i(:,m)=sig_i(:,m)/norm1;
%      alpha=norm1/(norm1+r);
%      mu_cap(:,m)=(1-alpha)*u(m,:)'+(alpha)*mu_i(:,m)/norm1;   % gaussian component wise updated mean
%      sig_cap(1,:,m)=(1-alpha)*cov_m(:,:,m)+(alpha)*sig_i(1,:,m)/norm1-(mu_cap(:,m).^2)';      % '_ _ _'    updated var
%      w1(m)=(1-alpha)*w(m);
%      w2(m)=(alpha)*norm1;
%     end
%  
%     mean=[mean;mu_cap];      % total mean matrix for all writers 2D
%     covar{i}=sig_cap;   % total cov matrix for all writers 3D
%     wt=[wt;w1+w2/norm2];
% end
%  
% %=========================================================================
%  %% Testing from the remqaing data
%  
% %  Likelihood=zeros(1,20);
% %  
% % for i=1:20
% %     l=find(TestD_20(:,6)==i);
% %    norm1=0;
% %    for t=l(1):l(end) %length(find(TestD_20(:,6)==i))
% %        norm2=0;
% %        for m=1:k
% %            var=covar{i};
% %            norm2=norm2+wt(i,m)*mvnpdf(TestD_20(t,1:5),mean((i-1)*5+1:i*5,m)',diag(var(1,:,m)));
% %        end
% %        norm1=norm1+log(norm2);
% %    end
% %    Likelihood(i)=norm1;
% % end
% 
% [idx,ll,P]=cluster()
toc




