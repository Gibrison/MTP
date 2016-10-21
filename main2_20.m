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
%=========================================================================
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
%c=0;s1=0;
ftr_mat=[];
tgap=[];
for i=1:20
    len1=find(M(:,7)==i);
    %xx1=M(len1,6);
   % len1(end)
    %temp=M(len1,6);
    len2=unique(M(len1,6));
   
    for j=1:len2(end)
        
%          for i=len2        % to check the above concept of for loop
%            i
%          end
         %c=c+1;
        l3=find(M(:,6)==j & M(:,7)==i);
      %  check(c,1:2)=[l3(1),l3(end)];
        %len3(end);
        %s2=length(len3);
        %s1=s1+s2;
        %c=c+1;
        %ends(c,1)=s1;
        mat(1:length(l3)-1,1:7)=0;
        
        t_gap(1:length(l3)-1,1)=(M(l3(2):l3(end),3)-M(l3(1):l3(end)-1,3));
        % distance
        mat(1:length(l3)-1,1)=sqrt((M(l3(2):l3(end),1)-M(l3(1):l3(end)-1,1)).^2+(M(l3(2):l3(end),2)-M(l3(1):l3(end)-1,2)).^2); 
        mat(1:length(l3)-1,2)=(M(l3(2):l3(end),3)-M(l3(1):l3(end)-1,3));  % sine theta
        mat(1:length(l3)-1,3)=(M(l3(2):l3(end),2)-M(l3(1):l3(end)-1,2)); % cos Theta
       % mat(1:length(l3)-1,4:5)=0;
        mat(1:length(l3)-2,4)=mat(1:length(l3)-2,3).*mat(2:length(l3)-1,3)+mat(1:length(l3)-2,2).*mat(2:length(l3)-1,2);
        mat(1:length(l3)-2,5)=mat(1:length(l3)-2,3).*mat(2:length(l3)-1,2)+mat(1:length(l3)-2,2).*mat(2:length(l3)-1,3);
       
        mat(1:length(l3)-1,6:7)=M(l3(1):l3(end)-1,6:7);
        
        %sin_p(s1-s2+1:s1-2,1)=cos_t(s1-s2+1:s1-2,1).*sin_t(s1-s2+2:s1-1,1)-sin_t(s1-s2+1:s1-2,1).*cos_t(s1-s2+2:s1-1,1);
        % clear len3
        %label(s1-s2+1:s1-1,1)=M(s1-s2+1:s1-1,7);    % confirm that...assignment is corect 
        %
      ftr_mat=[ftr_mat;mat];
      tgap=[tgap;t_gap];
      t_gap=[];
      mat=[];
        
    end
    
    % len1
    % temp
     
end
clear t_gap mat
%zeros=find(tgap==0);
%ftr_mat(:,2)=ftr_mat(:,2)./ftr_mat(:,1);  % sine theta

%ftr_mat(:,3)=ftr_mat(:,3)./ftr_mat(:,1);   % cos theta
ftr_mat(find(ftr_mat(:,1)==0),2)=0;
ftr_mat(find(ftr_mat(:,1)==0),3)=1;

ftr_mat(:,4)=ftr_mat(:,4)./ftr_mat(:,1).^2;  % sine phi

ftr_mat(:,5)=ftr_mat(:,5)./ftr_mat(:,1).^2;  % cos phi

%ftr_mat(:,1)=ftr_mat(:,1)./tgap;  % velocity

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
    l1=find(ftr_mat(:,6)==1 & ftr_mat(:,7)==i); l4=find(ftr_mat(:,6)==4 & ftr_mat(:,7)==i);
    l5=find(ftr_mat(:,6)==5 & ftr_mat(:,7)==i); l8=find(ftr_mat(:,6)==8 & ftr_mat(:,7)==i);
%     l1(1);
%     l4(end);
     TrainD_20=[TrainD_20;ftr_mat(l1(1):l4(end),[1:5,7])];
     TestD_20=[TestD_20;ftr_mat(l5(1):l8(end),:)];
end


%% Training from the writer data
%=========================================================================
% 
% %Tr=TrainD_20(:,1:5);
% %Ts=TestD_20;
% k=50;   % no of gaussian components
% 
% %GMModel1 = fitgmdist(Tr,k,'CovarianceType','diagonal','RegularizationValue',0.0002);
% opt = statset('MaxIter',1000);
% 
% GMModel = gmdistribution.fit(TrainD_20(:,1:5),k,'CovType','diagonal','Regularize',0.0002,'Options',opt);
% 
% covm=GMModel.Sigma;
% mu=GMModel.mu;
% w=GMModel.ComponentProportion;
% % 
% %  
% %  t_gap(ro,1)=0;
% %  ftr_mat(ro,1:3)=0;
% %  ftr_mat(ro,1)=0;
% %  ftr_mat(ro,1)=0;
% %  
%  %zero_t=find(t_gap==0);
%  
%  % solving the issue of zero distance 
%  
% % zero_dist=find(dist==0);
%  ftr_mat(:,1)=ftr_mat(:,1)./ftr_mat;
%  ftr_mat=ftr_mat./ftr_mat;
%  ftr_mat(find(ftr_mat==0))=1;
%  ftr_mat(find(ftr_mat==0))=0;
%  
%  cos_p(:,1)=ftr_mat(1:ro-2,1).*ftr_mat(2:ro-1,1)+ftr_mat(1:ro-2,1).*ftr_mat(2:ro-1,1);
%  sin_p(:,1)=ftr_mat(1:ro-2,1).*ftr_mat(2:ro-1,1)-ftr_mat(1:ro-2,1).*ftr_mat(2:ro-1,1);
%  cos_p(ro-1:ro,1)=0;
%  sin_p(ro-1:ro,1)=0;
%  
% %  cos_p(ends)=[];
% %  sin_p(ends)=[];
% %  t_gap(ends)=[];
% %  dist(ends)=[];
% %  cos_t(ends)=[];
% %  sin_t(ends)=[];
% % %  
% pt_ftr20(:,1) = ftr_mat./t_gap;      %% velocity feature
% pt_ftr20(:,2) = ftr_mat;            %% cos theta angle
% pt_ftr20(:,3) = ftr_mat;            %% sin theta angle
% pt_ftr20(:,4) = cos_p;            %% cos phi angle
% pt_ftr20(:,5) = sin_p;            %% sine phi angle
% pt_ftr20(:,6) = M(1:ro,6);
% pt_ftr20(:,7) = M(1:ro,7);
% % pt_feature20(:,8) = M(1:ro-c,6);
% % pt_feature20(:,9) = M(1:ro-c,7);
% %pt_feature20(ro,:)=0;
% %pt_feature20(ends,:)=[];
% % %zz=find(M(11523:14104,6)==1);
% pt_ftr20(ends,:)=[];