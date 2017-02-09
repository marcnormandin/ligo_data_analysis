%Reading information of all the datafiles
clear; clc;
Files = dir('data_snr*.mat');

%Open a textfile to save filenames as strings. 
fid = fopen('filenames.txt', 'wt');

for i = 1:length(Files)
    %FileNames(i,:)=Files(i).name;
    %tmp=FileNames(i,:);
    tmp=Files(i).name;
    
    % .mat extension will be discared here
    %tmp_name(i,:)=tmp(1:end-4);
     tmp_name = tmp(1:end-4);
     
    fprintf(fid, '%s\n', tmp_name);
end
fclose(fid);