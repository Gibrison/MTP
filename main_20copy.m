%% creating feature matrix : point based features
%========================================================================


clc;
clear all;
close all;

disp('program : main_20copy is running...................');
disp(' ');

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

 %M(:,1)=(M(:,1)-sum(M(:,1))/ro)/sqrt(sum((M(:,1)-sum(M(:,1))/ro).*(M(:,1))/ro));
% 
% %my=sum(M(:,2))/ro;
% %sig_y=sqrt(sum((M(:,2)-sum(M(:,2))/ro).*(M(:,2))/ro));

 %M(:,2)=(M(:,2)-sum(M(:,2))/ro)/sqrt(sum((M(:,2)-sum(M(:,2))/ro).*(M(:,2))/ro));
 
%===================================================================
% % min-max method
% 
 %M(:,1)=M(:,1)/(max(M(:,1))-min(M(:,1)));
% 
% M(:,2)=M(:,2)/(max(M(:,2))-min(M(:,2)));
%====================================================================
one=[];two=[];three=[];four=[];five=[];six=[];seven=[];
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
        
        t_gap=(M(s1-s2+2:s1,3)-M(s1-s2+1:s1-1,3));
        
        dist=sqrt((M(s1-s2+2:s1,1)-M(s1-s2+1:s1-1,1)).^2+(M(s1-s2+2:s1,2)-M(s1-s2+1:s1-1,2)).^2);
        cos_t=(M(s1-s2+2:s1,1)-M(s1-s2+1:s1-1,1));
        sin_t=(M(s1-s2+2:s1,2)-M(s1-s2+1:s1-1,2));
        place=find(dist==0);
        cos_t(place)=1;
        sin_t(place)=0;
        cos_p=cos_t(1:s2-2,1).*cos_t(2:s2-1,1)+sin_t(1:s2-2,1).*sin_t(2:s2-1,1);
        sin_p=cos_t(1:s2-2,1).*sin_t(2:s2-1,1)-sin_t(1:s2-2,1).*cos_t(2:s2-1,1);
        % clear len3
        %label(s1-s2+1:s1-1,1)=M(s1-s2+1:s1-1,7);    % confirm that...assignment is corect 
        vel=dist./t_gap;
        
        one=[one;vel(1:s2-2)];
        two=[two;sin_t(1:s2-2)];
        three=[three;cos_t(1:s2-2)];
        four=[four;sin_p];
        five=[five;cos_p];
        six=[six;M(s1-s2+1:s1-2,6)];
        seven=[seven;M(s1-s2+1:s1-2,7)];
        
      
        
    end
    
    % len1
    % temp
     
end
% xx=[1 1 1 1]';
% y=[2 2 2 2]';
% yy=[xx;y]
 matx=[one,two,three,four,five,six,seven];
 clear one two three four five six seven sin_t cos_t sin_p cos_p vel s1 s2 temp len1 len2 len3
 zeroos=find(matx(:,1)>55);
 matx( ~any(matx,2), : ) = [] ;%rows
%  place=find(matx(:,1)==0);
%  matx(place,2)=0;
%  matx(place,3)=1;
%  matx(place,2)=0; 
 
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
%% trainig and testing data generation from the faeture matrix (pt_ftr20)
%=========================================================================
pt_ftr20=matx;
%Tr=[];Ts=[];
TrainD_20=[];TestD_20=[];
%pt_ftr20=pt_ftr20;
%s3=0;s5=0;
for i=1:20
%     len4=find(D(:,7)==i);
%     
%     
%      
    l1=find(pt_ftr20(:,6)==1 & pt_ftr20(:,7)==i); l4=find(pt_ftr20(:,6)==6 & pt_ftr20(:,7)==i);
    l5=find(pt_ftr20(:,6)==7 & pt_ftr20(:,7)==i); l8=find(pt_ftr20(:,6)==8 & pt_ftr20(:,7)==i);
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

disp('Training started........................................');
disp(' ')

Tr=TrainD_20(:,1:5);
%Ts=TestD_20;
k=50;   % no of gaussian components

%GMModel1 = fitgmdist(Tr,k,'CovarianceType','diagonal','RegularizationValue',0.0002);
 opt = statset('MaxIter',1000);
GMModel2 = gmdistribution.fit(Tr,k,'CovType','diagonal','Regularize',0.0002,'Options',opt);

cov_m=GMModel2.Sigma;
u=GMModel2.mu;
w=GMModel2.ComponentProportion;

r=10;  % relevance factor
mean=[];
covar=[];
wt=[];
w2=zeros(1,20);
w1=w2;

disp('adaptation started.................');
disp(' ');

for i=1:20
     mu_i=zeros(5,k);
     sig_i=zeros(1,5,k);
     norm2=0;
     l=find(TrainD_20(:,6)==i);
    for m=1:k
    
     norm1=0;
     for t=l(1):l(end)      %length(find(TrainD_20(:,6)==i))
        mu_i(:,m)=mu_i(:,m)+ (mvnpdf(Tr(t,:),u(m,:),diag(cov_m(:,:,m)))*Tr(t,:)*w(m))';      % column format
        sig_i(1,:,m)=sig_i(1,:,m)+ (mvnpdf(Tr(t,:),u(m,:),diag(cov_m(:,:,m)))*Tr(t,:).^2*w(m));   % row format
        norm1=norm1+mvnpdf(Tr(t,:),u(m,:),diag(cov_m(:,:,m)))*w(m);          % constant number
        norm2=norm2+mvnpdf(Tr(t,:),u(m,:),diag(cov_m(:,:,m)))*w(m);
        
     end
    % mu_i(:,m)=mu_i(:,m)/norm1;
    % sig_i(:,m)=sig_i(:,m)/norm1;
     alpha=norm1/(norm1+r);
     mu_cap(:,m)=(1-alpha)*u(m,:)'+(alpha)*mu_i(:,m)/norm1;   % gaussian component wise updated mean
     sig_cap(1,:,m)=(1-alpha)*(cov_m(:,:,m)+u(m,:).^2)+(alpha)*(sig_i(1,:,m)/norm1+u(m,:).^2)-(mu_cap(:,m).^2)';      % '_ _ _'    updated var
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
 %==========================================================================
 disp('Testing Started...................................................');
 
 Likelihood=[];
 
for i=1:20
    %l=find(TestD_20(:,6)==i);
   
   for j=7:8
       
      l=find(TestD_20(:,6)==j & TestD_20(:,7)==i);
      
      for wr=1:20
          norm1=0;
         for t=l(1):l(end) %length(find(TestD_20(:,6)==i))
           norm2=0;
           
           for m=1:k
               
            var=covar{wr};
            norm2=norm2+wt(wr,m)*mvnpdf(TestD_20(t,1:5),mean((wr-1)*5+1:wr*5,m)',diag(var(1,:,m)));
            
           end
           
          norm1=norm1+log(norm2);
          
         end
         
        LL(j,2+wr)=norm1;
        
      end
      LL(j,1)=i;
      LL(j,2)=j;
   end
   Likelihood=[Likelihood;LL];
   done=i
end

%Likelihood=Likelihood;
Likelihood( ~any(Likelihood,2), : ) = [] ;%rows
LLR=Likelihood(:,1:2);
ac=0;
for i=1:40
    [maxx,loc]=max(Likelihood(i,3:22));
    LLR(i,3)=loc;
    if(LLR(i,1)==LLR(i,3))
        ac=ac+1;
    end
end
%  [model] = trainmsvm(TrainD_20(:,1:5), TrainD_20(:,6)) ;
%  [labels, outputs] = predmsvm(model, TestD_20(:,1:5), TestD_20(:,7), ncpus);
clear c co cos_t cos_p sin_t sin_p ends i j l1 l4 l5 l8 len1 len2 len3 m s1 s2 ro sig_i temp

toc

 %fprintf(2,'my text...\n')     % print text in colors