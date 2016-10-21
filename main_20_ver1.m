%% main program for writer identification

clc;
clear all;
close all;

disp('program : main_20_ver1 is running............');
disp(' ');

load nw_data2

l1=find(data_mat(:,7)==31);
l2=find(data_mat(:,7)==50);
M=data_mat(l1(1):l2(end),:);
M(:,7)=M(:,7)-30;

eps=0.05;  % vulue assigned to zero time gaps

wn=20;  % no of writers

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
       
       a1=find(M(:,7)==i& M(:,6)==j);  % no of stroke in a file
       r=unique(M(a1,5));
       
       for i1=1:length(r)
           s1=find(M(:,7)==i & M(:,6)==j & M(:,5)==i1);  % no of points in a stroke
           
           if length(s1)==1   % if an stroke cosists single point
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
           
           dist=dist/10000;   % normalization
           %(dist-sum(dist)/length(dist))/sum((dist-sum(dist)/length(dist)).^2)*length(dist);
           
           end
           
           tgap=[tgap;tgp];     % time gap
           
           tgp(find(tgp==0))=eps;   % time gap zero instaces removal
           
           distt=[distt;dist./tgp];  % velocity
           
           sin_t=[sin_t;del_y];   % sine theta
           
           cos_t=[cos_t;del_x];    % cos theta
           
           sin_p=[sin_p;sinp];    % sine phi
           
           cos_p=[cos_p;cosp];   % cos phi
           
       end
       
   end
   
end
%velo=distt./tgap;
%distt=distt;
% storing every faeture in a single matrix name mat
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

toc
%A=find(M(:,7)==1 & M(:,6)==1 & M(:,5)==1);

%aa=find(distt>100);
mat(find(abs(mat(:,1))>20),:)=[];
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
     Tr20=[Tr20;pt_ftr20(l1(1):l4(end),[1:5,7])];  % Training data set 
     
     Ts20=[Ts20;pt_ftr20(l5(1):l8(end),:)];     % Testing data set
end
% 
% TrainD_20(:,6)=[];
% TestD_20(:,6)=[];
% 
% 
%%========================================================================
%% Training from the writer data
%=========================================================================
tic

disp('UBM modelling started........');
disp(' ')

Tr=Tr20(:,1:5);
%Ts=TestD_20;
k=50;   % no of gaussian components

% UBM paramaeters 
cov_m=zeros(1,5,k);
u=zeros(k,5);
w=zeros(1,k);

%GMModel1 = fitgmdist(Tr,k,'CovarianceType','diagonal','RegularizationValue',0.0002);

 opt = statset('MaxIter',1000);
GMModel = gmdistribution.fit(Tr,k,'CovType','diagonal','Regularize',0.0002,'Options',opt); %'CovType','diagonal',

cov_m=GMModel.Sigma;
u=GMModel.mu;
w=GMModel.ComponentProportion;

r=5;         % relevance factor

 %[idx,nlogl,P] = cluster(GMModel,TrainD_20(1:200,1:5));
 
%  [idx,nlogl,P] = cluster(GMModel,Ts20(1:100,1:5));
%  for i=1:100
%     for j=1:k
%        pp(i,j)=mvnpdf(Ts20(i,1:5),u(j,:),diag(cov_m(1,:,j)))*w(j); 
%     end
%  end
%  check=abs(P)-abs(pp);
%  x=sum(log(sum(pp')))
%  xx=log(sum(P'));
%PP=0;

disp('adaptation started..(parameter learning.....)...............');
covv=cell(1,wn);
mean=cell(1,wn);
wt=cell(1,wn);
alpha=zeros(1,k);
me=zeros(k,5);
covva=zeros(1,5,k);

 for i=1:wn
    l=find(Tr20(:,6)==i); 
   
    c=0;
    %clear PP
     PP=zeros(length(l),k);
    for i1=l(1):l(end)
        c=c+1;
       for j=1:k
           
          PP(c,j)=mvnpdf(Tr(i1,:),u(j,:),diag(cov_m(1,:,j)))*w(j); 
       end
    end
    alpha = sum(PP)./(sum(PP)+r);
    wt{i} = (1-alpha).*w + alpha.*(sum(PP)/length(l));
   % mean{i} = (1-alpha)*u + alpha.*((PP'*Tr(l(1):l(end),:))./sum(PP));
    for j=1:k
        
       me(j,:) = (1-alpha(j))*u(j,:) + alpha(j)*((PP(:,j)'*Tr(l(1):l(end),:))/sum(PP(:,j))); 
    covva(1,:,j)= alpha(j)* ((PP(:,j)'*(Tr(l(1):l(end),:).^2))/sum(PP(:,j)))+(1-alpha(j))*(cov_m(1,:,j)+u(j,:).^2)-(me(j,:)).^2;
    
    end
    covv{i} = covva;
    mean{i}=me;
    %clear PP
    writer=i  
    
 end

  toc
  
  clear l l1 l4 l5 l8 covva me 
  
 %%=========================================================================
 %% Testing from the remaining data 
 %==========================================================================
 tic
 
 disp('Testing started.........................');
 
 idx=0;
 LL=zeros(2*wn,3);
 LogLL=zeros(2*wn,wn);
 u1=zeros(k,5);
 var=zeros(1,5,k);
 wtt=zeros(1,k);
 
 for i=1:wn
     
    for j=7:8
        
       l=find(Ts20(:,6)==j & Ts20(:,7)==i);
       ar=zeros(1,wn);
       for id=1:wn
           PP2=zeros(length(l),k);
           c=0;
           
           for t=l(1):l(end)
               c=c+1;
              for j1=1:k
                  
                  u1=mean{id};
                  var=covv{id};
                  wtt=wt{id};
                 PP2(c,j1)=mvnpdf(Ts20(t,1:5),u1(j1,:),diag(var(1,:,j1)))*wtt(j1);
                 
              end
              
           end
           ar(1,id)=sum(log(sum(PP2')));
           
       end
       [mx,loc]=max(ar);
       idx=idx+1;
       LL(idx,1)=j;
       LL(idx,2)=i;
       LL(idx,3)=loc;
       LogLL(idx,:)=ar;
    end
    
    done_writer=i
    
 end
 
 
 result=find(LL(:,2)==LL(:,3));
 
 clear wtt u1 var idx loc t sin_t cos_t sin_p cos_p a1 c ends Tr w Tr cov_m alpha
 fprintf(2,'Task completed..!!')
 disp(' ');
 toc
%%=========================================================================
%% Results Discussion
%==========================================================================

% k=50 , r=16 , accuracy=70
% k=40 , r=16 , accuracy=75
% k=40 , r=5 , accuracy =80
% k=50 , r=5 , acc= 92.50 % for 20 writers and 100 % for 10 writers

%display('result')
%fprintf(2,'hello world..\n')



