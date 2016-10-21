clc;
clear all;
close all;


load nw_data2

l1=find(data_mat(:,7)==31);
l2=find(data_mat(:,7)==50);
M=data_mat(l1(1):l2(end),:);
M(:,7)=M(:,7)-30;
eps=0.01;

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

for i=1:20
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
           
           dist=(dist-sum(dist)/length(dist))/sum((dist-sum(dist)/length(dist)).^2)*length(dist);
           
           end
           tgap=[tgap;tgp];
           distt=[distt;dist./(tgp+eps)];
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

toc
%A=find(M(:,7)==1 & M(:,6)==1 & M(:,5)==1);

%aa=find(distt>100);
