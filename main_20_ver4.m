% %% creating feature matrix : point based features
% %========================================================================
% % in this we got 71% accuracy
% clc;
% clear all;
% close all;
% tic
% load refine_20
% M=refine_20;
% [ro,co]=size(M);
% M(:,7)=M(:,7)-30;
% 
% %% normalization      
% %=========================================================================
% % mean-variance method
% % %mx=sum(M(:,1))/ro;
% % %sig_x=sqrt(sum((M(:,1)-mx).*(M(:,1))/ro));
% 
%  %M(:,1)=(M(:,1)-sum(M(:,1))/ro)/sqrt(sum((M(:,1)-sum(M(:,1))/ro).*(M(:,1))/ro));
% % 
% % %my=sum(M(:,2))/ro;
% % %sig_y=sqrt(sum((M(:,2)-sum(M(:,2))/ro).*(M(:,2))/ro));
% 
%  %M(:,2)=(M(:,2)-sum(M(:,2))/ro)/sqrt(sum((M(:,2)-sum(M(:,2))/ro).*(M(:,2))/ro));
%  
% %===================================================================
% % % min-max method
% % 
%  M(:,1)=M(:,1)/(max(M(:,1))-min(M(:,1)));
% % 
%  M(:,2)=M(:,2)/(max(M(:,2))-min(M(:,2)));
% %====================================================================
% c=0;s1=0;
% for i=1:20
%     len1=find(M(:,7)==i);
%    % len1(end)
%     temp=M(len1(1):len1(end),6);
%     len2=unique(temp(:,1));
%    
%     for j=1:len2(end)
%         len3=find(temp(:,1)==j);
%         %len3(end);
%         s2=length(len3);
%         s1=s1+s2;
%         c=c+1;
%         ends(c,1)=s1;
%         
%         t_gap(s1-s2+1:s1-1,1)=(M(s1-s2+2:s1,3)-M(s1-s2+1:s1-1,3));
%         
%         dist(s1-s2+1:s1-1,1)=sqrt((M(s1-s2+2:s1,1)-M(s1-s2+1:s1-1,1)).^2+(M(s1-s2+2:s1,2)-M(s1-s2+1:s1-1,2)).^2);
%         cos_t(s1-s2+1:s1-1,1)=(M(s1-s2+2:s1,1)-M(s1-s2+1:s1-1,1));
%         sin_t(s1-s2+1:s1-1,1)=(M(s1-s2+2:s1,2)-M(s1-s2+1:s1-1,2));
%         %cos_p(s1-s2+1:s1-2,1)=cos_t(s1-s2+1:s1-2,1).*cos_t(s1-s2+2:s1-1,1)+sin_t(s1-s2+1:s1-2,1).*sin_t(s1-s2+2:s1-1,1);
%         %sin_p(s1-s2+1:s1-2,1)=cos_t(s1-s2+1:s1-2,1).*sin_t(s1-s2+2:s1-1,1)-sin_t(s1-s2+1:s1-2,1).*cos_t(s1-s2+2:s1-1,1);
%         % clear len3
%         %label(s1-s2+1:s1-1,1)=M(s1-s2+1:s1-1,7);    % confirm that...assignment is corect 
%         %
%       
%         
%     end
%     
%     % len1
%     % temp
%      
% end
% 
%  
%  
%  t_gap(ro,1)=0;
%  dist(ro,1)=0;
%  cos_t(ro,1)=0;
%  sin_t(ro,1)=0;
%  
%  %zero_t=find(t_gap==0);
%  
%  % solving the issue of zero distance 
%  
% % zero_dist=find(dist==0);
%  cos_t=cos_t./dist;
%  sin_t=sin_t./dist;
%  cos_t(find(dist==0))=1;
%  sin_t(find(dist==0))=0;
%  
%  cos_p(:,1)=cos_t(1:ro-2,1).*cos_t(2:ro-1,1)+sin_t(1:ro-2,1).*sin_t(2:ro-1,1);
%  sin_p(:,1)=cos_t(1:ro-2,1).*sin_t(2:ro-1,1)-sin_t(1:ro-2,1).*cos_t(2:ro-1,1);
%  cos_p(ro-1:ro,1)=0;
%  sin_p(ro-1:ro,1)=0;
%  
% %  cos_p(ends)=[];
% %  sin_p(ends)=[];
%   t_gap(ends)=10^7;
% %  dist(ends)=[];
% %  cos_t(ends)=[];
% %  sin_t(ends)=[];
% % %  
% pt_ftr20(:,1) = dist./t_gap;      %% velocity feature
% pt_ftr20(:,2) = cos_t;            %% cos theta angle
% pt_ftr20(:,3) = sin_t;            %% sin theta angle
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
% %=========================================================================
% %% trainig and testing data generation from the faeture matrix (pt_ftr20)
% %=========================================================================
% %Tr=[];Ts=[];
% TrainD_20=[];TestD_20=[];
% %pt_ftr20=pt_ftr20;
% %s3=0;s5=0;
% for i=1:20
% %     len4=find(D(:,7)==i);
% %     
% %     
% %      
%     l1=find(pt_ftr20(:,6)==1 & pt_ftr20(:,7)==i); l4=find(pt_ftr20(:,6)==4 & pt_ftr20(:,7)==i);
%     l5=find(pt_ftr20(:,6)==5 & pt_ftr20(:,7)==i); l8=find(pt_ftr20(:,6)==8 & pt_ftr20(:,7)==i);
% %     l1(1);
% %     l4(end);
%      TrainD_20=[TrainD_20;pt_ftr20(l1(1):l4(end),[1:5,7])];
%      TestD_20=[TestD_20;pt_ftr20(l5(1):l8(end),:)];
% end
% % 
% % TrainD_20(:,6)=[];
% % TestD_20(:,6)=[];
% 
% 


%% main program for writer identification

clc;
clear all;
close all;

disp('program : main_20_ver4 is running.............................','r');
disp(' ');

load nw_data2

l1=find(data_mat(:,7)==31);
l2=find(data_mat(:,7)==40);
M=data_mat(l1(1):l2(end),:);
M(:,7)=M(:,7)-30;
eps=0.01;
wn=10;
%% normalization
%==============================

tic

%================================
%mat(:,1)=sqrt((M(2:end,1)-M(1:end-1,1)).^2+(M(2:end,2)-M(1:end-1,2)).^2);
%zeros=find(mat(:,1)==0);

% tgap=M(2:end,3)-M(1:end-1,3); % no element is zero
% zeros=find(tgap==0);

% mat(:,2)=(M(2:end,2)-M(1:end-1,2))./mat(:,1);
% mat(:,3)=(M(2:end,1)-M(1:end-1,1))./mat(:,1);
% mat(zeros,2)=0;
% mat(zeros,3)=1;
c=0;
distt=[];
sin_t=[];
cos_t=[];
sin_p=[];
cos_p=[];
tgap=[];

for i=1:wn
   for j=1:8
       a1=find(M(:,7)==i& M(:,6)==j);
       r=unique(M(a1,5));
       for i1=1:length(r)
           s1=find(M(:,7)==i & M(:,6)==j & M(:,5)==i1);
           if length(s1)==1
%                tgp=1;
%                dist=0;
%                del_y=0;
%                del_x=1;
%                sinp=0;
%                cosp=1;
            c=c+1;
            ends(c)=s1;
           %elseif length(s1)==0
              continue;     %continue;
           else
               tgp=M(s1(2):s1(end),3)-M(s1(1):s1(end)-1,3);
 %M(s1(1):s1(end),1)=( M(s1(1):s1(end),1),1))/length(s1)).^2))/length(s1));
% M(s1(1):s1(end),2)=sum( M(s1(1):s1(end),2))/length(s1)).^2))/length(s1));
           dist=sqrt((M(s1(2):s1(end),1)-M(s1(1):s1(end)-1,1)).^2+(M(s1(2):s1(end),2)-M(s1(1):s1(end)-1,2)).^2);
           
           del_y=(M(s1(2):s1(end),2)-M(s1(1):s1(end)-1,2))./dist;
           del_x=(M(s1(2):s1(end),1)-M(s1(1):s1(end)-1,1))./dist;
           
           del_y(find(dist==0))=0;
           del_x(find(dist==0))=1;
           
           sinp=del_x(1:length(s1)-2).*del_y(2:length(s1)-1)-del_y(1:length(s1)-2).*del_x(2:length(s1)-1);
           cosp=del_y(1:length(s1)-2).*del_y(2:length(s1)-1)+del_x(1:length(s1)-2).*del_x(2:length(s1)-1);
           
           sinp(length(s1)-1,1)=0;
           cosp(length(s1)-1,1)=1;
           
           c=c+1;
           ends(c)=s1(end);
           
           dist=dist/150;     %(dist-sum(dist)/length(dist))/sum((dist-sum(dist)/length(dist)).^2)*length(dist);
           
           end
           tgap=[tgap;tgp];
           tgp(find(tgp==0))=0.05;
           distt=[distt;dist./(tgp)];
           sin_t=[sin_t;del_y];
           cos_t=[cos_t;del_x];
           sin_p=[sin_p;sinp];
           cos_p=[cos_p;cosp];
       end
   end
end
%velo=distt./tgap;

matt=M(:,7);
matt(ends)=[];
mat(:,7)=matt;
clear matt
matt=M(:,6);
matt(ends)=[];
mat(:,6)=matt;

mat(:,1)=distt;
mat(:,2)=sin_t;
mat(:,3)=cos_t;
mat(:,4)=sin_p;
mat(:,5)=cos_p;

clear del_x del_y l1 l2 tgp dist r cosp sinp i1 s1 matt

%mat(find(abs(mat(:,1))>20),:)=[];


toc
%A=find(M(:,7)==1 & M(:,6)==1 & M(:,5)==1);

%aa=find(distt>100);

%% trainig and testing data generation from the faeture matrix (pt_ftr20)
%=========================================================================
pt_ftr20=mat;
%Tr=[];Ts=[];
Tr20=[];Ts20=[];
%pt_ftr20=pt_ftr20;
%s3=0;s5=0;
for i=1:wn
%     len4=find(D(:,7)==i);
%     
%     
%      
    l1=find(pt_ftr20(:,6)==1 & pt_ftr20(:,7)==i); l4=find(pt_ftr20(:,6)==6 & pt_ftr20(:,7)==i);
    l5=find(pt_ftr20(:,6)==7 & pt_ftr20(:,7)==i); l8=find(pt_ftr20(:,6)==8 & pt_ftr20(:,7)==i);
%     l1(1);
%     l4(end);
     Tr20=[Tr20;pt_ftr20(l1(1):l4(end),[1:5,7])];
     Ts20=[Ts20;pt_ftr20(l5(1):l8(end),:)];
end
% 
% TrainD_20(:,6)=[];
% TestD_20(:,6)=[];
% 
%=========================================================================
%% Training from the writer data
%=========================================================================
disp('Training starts..................');
disp(' ');
Tr=Tr20(:,1:5);
%Ts=TestD_20;
k=50;   % no of gaussian components

%GMModel1 = fitgmdist(Tr,k,'CovarianceType','diagonal','RegularizationValue',0.0002);
 opt = statset('MaxIter',1000);
GMModel2 = gmdistribution.fit(Tr,k,'CovType','diagonal','Regularize',0.0002,'Options',opt); %'CovType','diagonal',

cov_m=GMModel2.Sigma;
u=GMModel2.mu;
w=GMModel2.ComponentProportion;

disp('adaptation started.........................................');
disp(' ');
r=16;  % relevance factor
mean=[];
covar=[];
wt=[];
for i=1:wn
     mu_i=zeros(5,k);
     sig_i=zeros(1,5,k);
     norm2=0;
     l=find(Tr20(:,6)==i);
    for m=1:k
    
     norm1=0;
     for t=l(1):l(end)      %length(find(TrainD_20(:,6)==i))
        mu_i(:,m)=mu_i(:,m)+ (mvnpdf(Tr(t,:),u(m,:),diag(cov_m(1,:,m)))*Tr(t,:)*w(m))';      % column format
        sig_i(1,:,m)=sig_i(1,:,m)+ (mvnpdf(Tr(t,:),u(m,:),diag(cov_m(1,:,m)))*Tr(t,:).^2*w(m));   % row format
        norm1=norm1+mvnpdf(Tr(t,:),u(m,:),(cov_m(:,:,m)))*w(m);          % constant number
        norm2=norm2+mvnpdf(Tr(t,:),u(m,:),(cov_m(:,:,m)))*w(m);
        
     end
    % mu_i(:,m)=mu_i(:,m)/norm1;
    % sig_i(:,m)=sig_i(:,m)/norm1;
     alpha=norm1/(norm1+r);
     mu_cap(:,m)=(1-alpha)*u(m,:)'+(alpha)*mu_i(:,m)/norm1;   % gaussian component wise updated mean
     sig_cap(1,:,m)=(1-alpha)*(cov_m(1,:,m)+u(m,:).^2)+(alpha)*(sig_i(1,:,m)/norm1+u(m,:).^2)-(mu_cap(:,m).^2)';      %  updated var
     w1(m)=(1-alpha)*w(m);
     w2(m)=(alpha)*norm1;
    end
 
    mean=[mean;mu_cap];      % total mean matrix for all writers 2D
    covar{i}=sig_cap;   % total cov matrix for all writers 3D
    wt=[wt;w1+w2/norm2];
    writer=i
end
 
%=========================================================================
 %% Testing from the remaining data
 disp('Testing starts......................');
 Likelihood=[];
 
for i=1:wn
    %l=find(TestD_20(:,6)==i);
   
   for j=7:8
       
      l=find(Ts20(:,6)==j & Ts20(:,7)==i);
      
      for wr=1:wn
          norm1=0;
         for t=l(1):l(end) %length(find(TestD_20(:,6)==i))
           norm2=0;
           
           for m=1:k
               
            var=covar{wr};
            norm2=norm2+wt(wr,m)*mvnpdf(Ts20(t,1:5),mean((wr-1)*5+1:wr*5,m)',diag(var(1,:,m)));
            
           end
           
          norm1=norm1+log(norm2);
          
         end
         
        LL(j,2+wr)=norm1;
        
      end
      LL(j,1)=i;
      LL(j,2)=j;
   end
   Likelihood=[Likelihood;LL];
   done_with_wrietr=i
end

%Likelihood=Likelihood;
Likelihood( ~any(Likelihood,2), : ) = [] ;%rows
LLR=Likelihood(:,1:2);
ac=0;
for i=1:2*wn
    [maxx,loc]=max(Likelihood(i,3:12));
    LLR(i,3)=loc;
    if(LLR(i,1)==LLR(i,3))
        ac=ac+1;
    end
end
%  [model] = trainmsvm(TrainD_20(:,1:5), TrainD_20(:,6)) ;
%  [labels, outputs] = predmsvm(model, TestD_20(:,1:5), TestD_20(:,7), ncpus);
clear c co cos_t cos_p sin_t sin_p ends i j l1 l4 l5 l8 len1 len2 len3 m s1 s2 ro sig_i temp

toc

