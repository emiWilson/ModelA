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
N = data(s - 4);
T = data(s - 3);
skipPrint = data(s - 2);
dt = data(s - 1);
dx = data(s);


timesteps = T/skipPrint;

m = zeros(N,N);
index = 1;

d = cell(N, 1);
for (k = 1:timesteps)
    for (i = 1:N)
        for (j = 1:N)
            m(i,j) = data(index);
            index = index + 1;
        end
    end
    
    d{k} = m;
end

for(k = 1:timesteps)
    figure 
    imshow(d{k,1},'InitialMagnification', 200, 'DisplayRange',[-0.005 0.005]);
end

%for(k = 1:T)
   % figure
   % fname = sprintf('myfile%d.jpg', k);
   % s = imshow(d{k,1},'InitialMagnification', 100, 'DisplayRange',[-0.005 0.005]);
   % t = sprintf('dx=%0.2d, dt=%d, N=%d, T=%d', dx, 5, N, T);
   % title(t)
   % saveas(s,fname);
    
  
    


%end

%figure
%imshow(d{1,1},'InitialMagnification', 100, 'DisplayRange',[-0.005 0.005]);