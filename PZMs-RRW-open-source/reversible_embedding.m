function [xf,PSNR,Ok,per1,per2,P1,Z1,P2,Z2] = reversible_embedding(f, reversible_watermarking, Rows, Cols, mask) 
%ʹ��vasiliy�ķ�������Ƕ��
%���룺ͼ��f������ˮӡ��Ϣpl
%�������ˮӡͼ��xf��PSNR���Ƿ�ȫ��Ƕ��Ok��per1��per2�ֱ��ʾ��X����O����Ƕ���ʣ�P1��P2�Լ�Z1��Z2��ʱ���岻��
% xf = double(f);
% pl = rand(1,rate) > 0.5;
rate = numel(reversible_watermarking);
[xf,okC,per1,P1,Z1] = embedCrossOrDot...
    (f,reversible_watermarking(1,1:round(rate/2)),1, Rows, Cols, mask);                    %X��Ƕ��
[xf,okD,per2,P2,Z2] = embedCrossOrDot...
    (xf,reversible_watermarking(1,round(rate/2)+1:numel(reversible_watermarking)),2, Rows, Cols, mask);          %O��Ƕ��
if okD && okC
    Ok = 1;
else
    Ok = 0;
end
[PSNR,mse] = psnr(f,xf);
%%
function [xf,fullEmbed,per,P,Z] = embedCrossOrDot(f,pl,CorD, r, c, mask)
%��X������O������Ƕ��
%���룺ͼ��f����Ƕ����Ϣpl��XO����־CorD
%�������ˮӡͼ��xf���ܷ�ȫ��Ƕ��fullEmbed
xf = f;
per = 0;
% [r,c] = size(xf);
Pnum = numel(pl);
Pnum = Pnum + 34;                                                          %ʵ��Ƕ��pl��34��lsb
fullEmbed = 1;

[img_pre,img_err,err_sorted,seq_map] = sortErrWithSur(xf,CorD);            %���X��O��Ԥ��ֵ��Ԥ���ֵ��������Ԥ���ֵ��Ԥ���ֵ�����
% [img_pre,img_err,err_sorted,seq_map,Tu] = impRhom(xf,CorD);
seqToRC = zeros( (r-2)*(c-2)/2 , 2);                                       %�����ص���������ŵõ�������
tmp = r*c;

if Pnum>28000
    for i = 2:r-1
        for j = 2:c-1
            if seq_map(i,j) ~= tmp %&& mask(i,j) == 0 %�˴�����ж�������ֻȡ����Բ������ص� %����������Բ�����ɵ�ֵ֮����޷����������ˣ���ʱע��
                seqToRC( seq_map(i,j) , : ) = [i,j];
            end
        end
    end
else
    for i = 2:r-1
        for j = 2:c-1
            if seq_map(i,j) ~= tmp && mask(i,j) == 0 %�˴�����ж�������ֻȡ����Բ������ص� %����������Բ�����ɵ�ֵ֮����޷����������ˣ���ʱע��
                seqToRC( seq_map(i,j) , : ) = [i,j];
            end
        end
    end
end
id = find(seqToRC(:,1)==0 & seqToRC(:,2)==0);% ɾ������Բ�е���
seqToRC(id,:)=[];
err_sorted(:,id)=[];

lsb34 = zeros(1,34);                                                       %��ȡ�����ǰ34�����ص�lsb
for i = 1:34
    lsb34(i) = bitget(xf(seqToRC(i,1),seqToRC(i,2)),1);
end

Tn = int16(0);                                                             %��ʼ��Tn��Tp
Tp = int16(0);
H = hist(err_sorted(1,35:numel(err_sorted)),-255:255);
while sum(H(1,Tn+256:Tp+256)) < Pnum
    if abs(Tn) <= Tp
        Tn = Tn -1;
    else
        Tp = Tp + 1;
    end
end

L = zeros(3,10000);                                                        %��¼��Ҫλͼ�ĵ����ţ����޸Ĵ������Լ���Ƕ��λ����ƽ��λ
E = zeros(2,numel(err_sorted) - 34);                                       %��¼���޸����εĵ����ţ��Լ���Ƕ��λ����ƽ��λ
numE = 0;
numES = 0;
numL = 0;
while numE < Pnum + numL                                                   %��������Ƕ��ĵ���С��Ƕ����Ϣ��С+λͼ��Сʱ����
    for i = 35:numel(err_sorted)                                           %�Դ����Ϊ35�ĵ㣬��ʼ����
        [mdNum,EorS] = modifiedNum(img_pre(seqToRC(i,1),seqToRC(i,2)),err_sorted(i),Tn,Tp);%���㵱ǰ���ܹ��޸ļ��Σ��Ǵ���Ƕ��λ����ƽ��λ
        if mdNum == 2                                                      %�����޸�����
            numES = numES + 1;                                             %�����¼��E��
            E(1,numES) = i;
            E(2,numES) = EorS;
            if EorS == 1                                                   %�����ǿ�Ƕ��λ����numE����1
                numE = numE + 1;
            end
        else                                                               %�����޸�һ�λ����޸�
            numL = numL + 1;                                               %�����¼��L��
            L(1,numL) = i;                                                 
            L(2,numL) = mdNum;
            L(3,numL) = EorS;
        end
        if Pnum + numL == numE
            per = double(i)/numel(err_sorted);
            break;
        end
    end
    if Pnum + numL ~= numE                                                 %���������㣬������Tp����abs(Tn)
        if abs(Tn) <= Tp
            Tn = Tn - 1;
        else
            Tp = Tp + 1;
        end
        numE = 0;
        numES = 0;
        numL = 0;
    end
end
H = hist(err_sorted(1,35:34+numES+numL),-255:255);
P = Tn:1:Tp;
Z = [];
tmp = 1;
while H(tmp) == 0
    tmp=tmp+1;
end
tmp = tmp - 1;
if Tn ~= 0
    Z = tmp+Tn+1:1:tmp;
    Z = Z - 256;
end
tmp = 511;
while H(tmp) == 0
    tmp=tmp-1;
end
tmp = tmp+1;
tz = tmp:1:tmp+Tp;
tz = tz-256;
Z = [Z,tz];

msg = [L(2,1:numL),pl,lsb34];                                              %��ҪǶ���ȫ����Ϣ = λͼ+pl+lsb
count = 0;
for i = 1:numES                                                            %�Կ��޸����εĵ㣬����Ƕ���ƽ��
    err = err_sorted(E(1,i));
    if E(2,i) == 1
        count = count + 1;
        err = err*2 + msg(count);
    else
        if err < 0
            err = err + Tn;
        else
            err = err + Tp + 1;
        end
    end
    row = seqToRC(E(1,i),1);
    col = seqToRC(E(1,i),2);
    xf(row,col) = img_pre(row,col) + err;                                  %ʵ���޸�ͼ��
end
if count ~= numel(msg)
    fullEmbed = 0;
end
for i = 1:numL                                                             %����Ҫ��¼λͼ�ĵ�
    err = err_sorted(L(1,i));
    if L(2,i) == 1                                                         %�����޸�һ�Σ���������һ�η��������޸�
        if L(3,i) == 1                                                     %��Ƕ��λ
            if err < 0
                err = err*2;
            else
                err = err*2 + 1;
            end
        else                                                               %��ƽ��λ
            if err < 0 
                err = err + Tn;
            else
                err = err + Tp + 1;
            end
        end
    end
    row = seqToRC(L(1,i),1);
    col = seqToRC(L(1,i),2);
    xf(row,col) = img_pre(row,col) + err;                                  %ʵ���޸�ͼ��
end
lsb34 = [toBitSeq(uint16(abs(Tn)),7),toBitSeq(uint16(Tp),7),toBitSeq(uint32(Pnum - 34),20)];%�����滻��lsb = abs��Tn�� + Tp + Pnum - 34
for i = 1:34                                                               %lsb�滻���޸�ͼ��
    xf(seqToRC(i,1),seqToRC(i,2)) = bitset(xf(seqToRC(i,1),seqToRC(i,2)),1,lsb34(i));
end
%%
function [BS]  = toBitSeq(x,l)
%�������������Ϊ����������
%���룺����x�����������г���l
%���������������BS
BS = zeros(1,l);
x = abs(x);
for i=1:l
    BS(1,i) = bitget(x,i);
end
%%
function [num,EorS]  = modifiedNum(xpre,err,Tn,Tp)
%����һ��������޸ĵĴ������Լ���һ���޸���Ƕ�뻹��ƽ��
%���룺Ԥ��ֵxpre��Ԥ���ֵerr����ֵTn��Tp
%��������޸Ĵ���num��Ƕ��ƽ�Ʊ�־EorS��1����Ƕ�룬2����ƽ�ƣ�
num = 0;
EorS = 1;%E
if err > -1 && err <= Tp
    err1 = err*2 + 1;
    if xpre + err1 <= 255
        num = 1;
        
        if err1 <= Tp
            err2 = err1*2 +1;
        else
            err2 = err1 + Tp +1;
        end
        
        if xpre + err2 <= 255
            num = 2;
        end
    end
elseif err < 0 && err >= Tn
    err1 = err*2;
    if xpre + err1 >= 0
        num = 1;
        
        if err1 >= Tn
            err2 = err1*2;
        else
            err2 = err1 + Tn;
        end
        
        if xpre + err2 >=0
            num = 2;
        end
    end
elseif err > -1
    EorS = 2;
    err1 = err + Tp + 1;
    if xpre + err1 <= 255
        num = 1;
        err2 = err1 + Tp +1;
        if xpre + err2 <= 255
            num = 2;
        end
    end
else
    EorS = 2;
    err1 = err + Tn;
    if xpre + err1 >= 0
        num = 1;
        err2 = err1 + Tn;
        if xpre + err2 >= 0
            num = 2;
        end
    end
end
