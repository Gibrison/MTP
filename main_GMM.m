clc;
close all;
clear all;
 
tic
load TrainD_20.mat;
load TestD_20.mat;

Tr=TrainD_20(:,1:5);
Ts=TestD_20;

[ro,co]=size(TestD_20);

%P=C';
%k=32;
k=50;
%r=0.002;

%[u,L]=myKmeans(P,k);

%[u1,L1] = vector_quantizer(P,N,r)  ;

%[u,cov_mat,w] = GMM(P,N,r);
%obj = gmdistribution.fit(Tr,k);
GMModel = fitgmdist(Tr,k,'CovarianceType','diagonal','RegularizationValue',0.002);
obj = gmdistribution.fit(Tr,k,'CovType','diagonal','Regularize',0.0002);

%% adaptation
    % mu of client=alpha*muofUBM+(1-alpha)*sumover T(p(i/xt)*xt )/sumover T(p(i/xt));
    % w of client=alpha* w of UBM + (1-alpha)*
%     alpha=0.5;    % alpha=ni/(ni+r) where r is relevance factor r=16;
%     Mu_wr=[];
%     Wt=[];
%     Cov_M=[];
%   for i=1:20
%       T1=find(C(:,6)==i);
%       X1=C(T1(1):T1(end),1:5);
%       Mv=zeros(5,N);
%       norm=0;
%       for m=1:N     
%          % X1=repmat(X1(j,:),N,1)';
%          for j=1:size(X1,1)
%           Mv(:,m)= Mv(:,m)+(mvnpdf(X1(j,:),u(:,m)',cov_matrix(:,:,m))*X1(j,:)*w(m))';
%           norm=norm+(mvnpdf(X1(j,:),u(:,m)',cov_matrix(:,:,m)));
%          end
%       Mu_up(:,m)=alpha* u(:,m) + (1-alpha)*Mv(:,m);
%       Wt_up(m)=alpha* w + (1-alpha)*sum(mvnpdf(X1,u,cov_matrix))/sum(sum(mvnpdf(X1,u,cov_matrix)));
%       %Cov_M=[Cov_M;alpha* u + (1-alpha)*sum(mvnpdf(X1,u,cov_matrix)*X1.*X1)/sum(mvnpdf(X1,u,cov_matrix))]
%       end
%   end
%  %check=(mvnpdf(C(1,:),u(:,1)',cov_matrix(:,:,1))*C(1,:))'
%  obj = gmdistribution.fit(C,N);
 %%
%  mu1=zeros(5,N);
%  norm1=0;
%  for m=1:N
%      mu1=zeros(5,N);
%      norm1=0;
%      for j=1:size(X1,1)
%         mu1(:,m)=mu1(:,m)+ (mvnpdf(C(j,:),u(:,m)',cov_matrix(:,:,m))*C(j,:)*w(m))';
%         norm1=norm1+(mvnpdf(C(j,:),u(:,m)',cov_matrix(:,:,m)))*w(m);
%      end
%      mu1(:,m)=mu1(:,m)/norm1;
%  end
 
 
% 
% for i=1:20
%     len=find(TrainD_20(:,6)==i);
%     C=TrainD_20(len(1):len(end),1:5);
%     P=C';
%     [mu,cov_mat,wt] = GMM(P,N,r);
%     for j=1:ro
%         for k=1:32
%  f(k)=(1/(sqrt(det(cov_mat(:,:,k)))*(2*pi)^5)*exp(-0.5*(T(j,:)-mu(:,k)')*inv(cov_mat(:,:,k))*(T(j,:)-mu(:,k)')'))*wt(k);
%         end
%         Val(j,i)=sum(f);
%     end
% 
% 
% 
% 
% 
% end
toc


