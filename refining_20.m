%% doing for the 20 writers ..31 to 50..........

% 1. extract data of 20 writers from refine.
% 2. assign the point column to refine_20
% 3. calcilate time gap array and remove coinciding points
% 4. again perform point column assignment
% 5. now with that refine_20 get point based features
% 6. create featyre matrix for the 20 writers
% now your data is ready to fed to GMM(UBM)

clc;
clear all;
close all;
tic
load refine
tag_strt=find(refine(:,7)==31);
tag_end=find(refine(:,7)==50);
refine_20=refine(tag_strt(1):tag_end(end),:);
%refine_20(:,4)=column;
[ro ,co]=size(refine_20);
t_gap=refine_20(2:ro,3)-refine_20(1:ro-1,3);
zeros=find(t_gap==0);
refine_20(zeros+1,:)=0;

%% removing the zero rows and new column assignment

MM=refine_20;

 c=0;
 idx=0;
 writer=0;
for i=31:1:50
    files_per_Writer = size(getAllFiles(strcat('/home/suresh/Desktop/SUDHIR SHUKLA/MTP/IAM database by writer/original/',num2str(i))));
    %f=F(1:fileList(1));                                  
%     folderList = dir();
%     folderList(7).name
 fileList = getAllFiles(strcat('/home/suresh/Desktop/SUDHIR SHUKLA/MTP/IAM database by writer/original/',num2str(i)));
    for j=1:1:files_per_Writer(1)             % no of xml files per writer
       
          struct_data=xml2struct(fileList{j});
          data=struct_data.WhiteboardCaptureSession.StrokeSet.Stroke;
          n1=size(data); % no of strokes in j th file
         for k=1:n1(2)
              
             n2=size(data{k}.Point);  % no of points in k th stroke
             if n2(2)~=1
                 %for m=1:n2(2)
                    
                     
                     non_0=nnz(MM(1+c:c+n2(2)));
                     %l=(non_0);
                     %dif=n2(2)-l;
                     
                     column(1+idx:idx+non_0,1)=1:non_0;
                     idx=idx+non_0;
                     c=c+n2(2);
                     

                 %end
             else
                     non_0=nnz(MM(1+c:c+n2(2)));
                    % l=(non_0);
                     %dif=n2(2)-l;
                     
                     column(1+idx:idx+non_0,1)=1;
                     idx=idx+non_0;
                     c=c+n2(2);
                     

             end
            
             
             
         end
    
    end
    writer=writer+1
end

MM( ~any(MM,2), : ) = [] ;%rows
MM(:,4)=column(1:end-1,1);

refine_20=MM;
toc

%plot(MM2(28280:28320,1),-MM2(28280:28320,2),'.')

