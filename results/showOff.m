%fid = fopen('10x16.txt');
%miVals = fscanf(fid,'%lf',[10,16]);
%fclose(fid);
fileName = 'resOf256_new.txt';
miVals=importdata(fileName);
maxOfEachByte = eye(1,16);
disp(miVals(1));
for i=1:16
    col = miVals(:,i);
    [maax,index] = max(col);
    maxOfEachByte(i)=index-1;
    
end
disp(maxOfEachByte);
hold on;
xAxiss = 1:16;

grid on;
for i=1:16
    yAxiss = transpose(miVals(i,:));
    plot(xAxiss,yAxiss);
end
