function [ d_2 ] = signed_to_bin( d,n )
%�з�����dת������
%nΪλ��
d=d(:);
l=length(d);
k=1;
for i=1:l
    if d(i)>=0
        a=dec2bin(d(i),n);
    else
        a=dec2bin(2^n+d(i));
    end
    b=str2num(regexprep(a,'(?<=\d)(\d)',' $1'));
    d_2(k:k+n-1)=b;
    k=k+n;
end

