clc;
clear all;
close all;

disp('Program : main_20_v1 Started.........');
disp(' ');

tic

load nw_data2

l1=find(data_mat(:,7)==31);
l2=find(data_mat(:,7)==40);
M=data_mat(l1(1):l2(end),:);
M(:,7)=M(:,7)-30;

eps=0.05;  % vulue assigned to zero time gaps

wn=10;  % no of writers

c=0;
distt=[];
 sin_t=[];
% cos_t=[];
% sin_p=[];
% cos_p=[];
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
%            del_x=(M(s1(2):s1(end),1)-M(s1(1):s1(end)-1,1))./dist;
%            
            del_y(find(dist==0))=0;
%            del_x(find(dist==0))=1;
%            
%            sinp=del_x(1:length(s1)-2).*del_y(2:length(s1)-1)-del_y(1:length(s1)-2).*del_x(2:length(s1)-1);
%            cosp=del_y(1:length(s1)-2).*del_y(2:length(s1)-1)+del_x(1:length(s1)-2).*del_x(2:length(s1)-1);
%            
%            sinp(length(s1)-1,1)=0;
%            cosp(length(s1)-1,1)=1;
           
           c=c+1;
           ends(c)=s1(end);
           
           dist=dist/max(dist);   % normalization
           %(dist-sum(dist)/length(dist))/sum((dist-sum(dist)/length(dist)).^2)*length(dist);
           
           end
           
           tgap=[tgap;tgp];     % time gap
           
           tgp(find(tgp==0))=eps;   % time gap zero instaces removal
           
           distt=[distt;dist./tgp];  % velocity
           
            sin_t=[sin_t;del_y];   % sine theta
%            
%            cos_t=[cos_t;del_x];    % cos theta
%            
%            sin_p=[sin_p;sinp];    % sine phi
%            
%            cos_p=[cos_p;cosp];   % cos phi
           
       end
       
   end
   
end
%velo=distt./tgap;
%distt=distt;
% storing every faeture in a single matrix name mat
matt=M(:,7);
matt(ends)=[];
mat(:,5)=matt;

clear matt
matt=M(:,6);
matt(ends)=[];
mat(:,4)=matt;

matt=M(:,5);
matt(ends)=[];
mat(:,3)=matt;


mat(:,1)=distt;
mat(:,2)=sin_t;
% mat(:,2)=sin_t;
% mat(:,3)=cos_t;
% mat(:,4)=sin_p;
% mat(:,5)=cos_p;

toc
subplot(1,2,1);
stem(mat(1:16,1))
subplot(122), plot(mat(1:16,1))
% figure
% plot(mat(1:4745,1),mat(1:4745,2),'.')
% hold
% plot(mat(38497:42971,1),mat(38497:42971,2),'.')
% hold
% plot(mat(66004:72352,1),mat(66004:72352,2),'.')

fprintf(2,'Task completed...!!')