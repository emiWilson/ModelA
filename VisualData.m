clc; %clear the command window
close all; %close all figures
clear %erase all existing variables

%read data from file output.txt
%will make a array full of float values
fileID = fopen('output.txt','r');
formatSpec = '%f';
data = fscanf(fileID,formatSpec);
fclose(fileID);

%get sixe of file
s = size(data,1);

%extract important constants from file (appended to end)
N = data(s - 2);
T = data(s - 1);
dy = data(s);




m = zeros(N,N);
index = 1;

d = cell(N, 1);
for (k = 1:T)
    for (i = 1:N)
        for (j = 1:N)
            m(i,j) = data(index);
            index = index + 1;
        end
    end
    
    d{k} = m;
end

for(k = 1:T)
    figure 
    imshow(d{k,1},'InitialMagnification', 5000, 'DisplayRange',[-0.005 0.005]);
end
%for i = 1:N;
 %   imshow(d{N})
  %  pause(.001);
%end

