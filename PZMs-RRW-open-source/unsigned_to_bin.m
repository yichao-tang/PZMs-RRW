function [ d_2 ] = unsigned_to_bin( d,n )
%无符号数d转二进制
%n为位数
% d=d(:);
l=length(d);
k=1;
for i=1:l
    a=dec2bin(d(i),n);
    b=str2num(regexprep(a,'(?<=\d)(\d)',' $1')); %将二进制数字拆分到数组中
    d_2(k:k+n-1)=b;
    k=k+n;
end

