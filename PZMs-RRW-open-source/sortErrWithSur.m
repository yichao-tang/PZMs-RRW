function [img_pe,img_err,err_sorted,seq_map] = sortErrWithSur(img,xo)
[r,c] = size(img);
seq_map = ones(r,c)*r*c;
img = double(img);
img_pe = zeros(r,c);
img_err = img_pe;
err = zeros(1,(r-2) * (c-2)/2);
sur = err;
err_sorted = err;

if xo == 1
    offset = [0,0;1,1];
else
    offset = [1,0;0,1];
end

k = 1;
for i = 1:round(r/2) - 1
    for j = 1:round(c/2) - 1
        row = 2*i + offset(1,1);
        col = 2*j + offset(1,2);
        
        surround = [img(row,col-1),img(row-1,col),img(row,col+1),img(row+1,col)];
        dist_near = abs(diff([surround,img(row,col-1)]));
        img_pe(row,col) = round(mean(surround));%round(sum(surround)/4);
        img_err(row,col) = img(row,col) - img_pe(row,col);
        err(1,k) = img_err(row,col);%将预测误差存储成一个向量，方便直方图统计
%         sur(1,k) = mean((surround - img_pe(row,col)).^2);
        sur(1,k) = var(dist_near);
%         sur(1,k) = var(surround,1);
%         sur(1,k) = err(1,k);
        k = k + 1;
        
        row = 2*i + offset(2,1);
        col = 2*j + offset(2,2);
        
        surround = [img(row,col-1),img(row-1,col),img(row,col+1),img(row+1,col)];
        dist_near = abs(diff([surround,img(row,col-1)]));
        img_pe(row,col) = round(mean(surround));%round(sum(surround)/4);
        img_err(row,col) = img(row,col) - img_pe(row,col);
        err(1,k) = img_err(row,col);%将预测误差存储成一个向量，方便直方图统计
%         sur(1,k) = mean((surround - img_pe(row,col)).^2);
        sur(1,k) = var(dist_near);
%         sur(1,k) = var(surround,1);
%         sur(1,k) = err(1,k);
        k = k + 1;
    end
end
[sur,pos] = sort(sur);
for i = 1:(r-2)*(c-2)/2
    err_sorted(1,i) = err(1,pos(i));
    l = ceil(pos(i)/(c-2));
    j = ceil(mod(pos(i),c-2)/2);
    if j == 0
        j = (c-2)/2;
    end
    k = mod(pos(i)+1,2) + 1;
    row = 2*l + offset(k,1);
    col = 2*j + offset(k,2);
    seq_map(row,col) = i;
end