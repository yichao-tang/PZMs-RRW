function [ d_2 ] = unsigned_to_bin( d,n )
%�޷�����dת������
%nΪλ��
% d=d(:);
l=length(d);
k=1;
for i=1:l
    a=dec2bin(d(i),n);
    b=str2num(regexprep(a,'(?<=\d)(\d)',' $1')); %�����������ֲ�ֵ�������
    d_2(k:k+n-1)=b;
    k=k+n;
end

