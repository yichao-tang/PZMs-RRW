function [ d ] =bin_to_unsigned( d_2,n )
%d_2转十进制有符号数
%n为位数
l=length(d_2);
k=1;
for i=1:n:l
    a=d_2(i:i+n-1);
    b=num2str(a);
    c=strrep(b, ' ', '');%去掉空格
    d(k)=bin2dec(c);
    k=k+1;
end
