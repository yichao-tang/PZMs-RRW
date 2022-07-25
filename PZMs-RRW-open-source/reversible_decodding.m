function [xf,pl] = reversible_decodding(f, Rows, Cols, mask)
%ʹ��vasiliy�ķ������н���
%���룺��ˮӡͼ��f
%�����ԭͼ��xf����Ϣpl
xf = double(f);
[xf,pl2] = decodeCrossOrDot(xf,2, Rows, Cols, mask);                                         %��O���Ͻ��н���
[xf,pl1] = decodeCrossOrDot(xf,1, Rows, Cols, mask);                                         %��X���Ͻ��н���
pl =[pl1,pl2];
%%
function [xf,pl] = decodeCrossOrDot(f,CorD, r, c, mask)
%��X����O���Ͻ��н���
%���룺��ˮӡͼ��f��XO��ָʾCorD��1����X���ϣ�
%�����ȥ��X��O����ˮӡ��ͼ��xf,��Ϣpl
xf = f;
% [r,c] = size(xf);

[img_pre,img_err,err_sorted,seq_map] = sortErrWithSur(xf,CorD);            %���X��O��Ԥ��ֵ��Ԥ���ֵ��������Ԥ���ֵ��Ԥ���ֵ�����

seqToRC = zeros( (r-2)*(c-2)/2 , 2);                                       %�����ص���������ŵõ�������
tmp = r*c;
for i = 2:r-1
    for j = 2:c-1
        if seq_map(i,j) ~= tmp && mask(i,j) == 0 %�˴�����ж�������ֻȡ����Բ������ص�
            seqToRC( seq_map(i,j) , : ) = [i,j];            
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
Tn = deBitSeq(lsb34(1,1:7));                                               %��lsb�������ֵTn��Tp��Pnum
Tn = -int32(Tn);
Tp = int32(deBitSeq(lsb34(1,8:14)));
Pnum = deBitSeq(lsb34(1,15:34));

msg = zeros(1,(((r-2)*(c-2))-nnz(mask))/2);                                                      %������������Ϣ λͼ+pl+lsb
L = zeros(2,1000000);                                                      %��¼��Ҫλͼ�ĵ㣬�Լ����Ǳ�Ƕ�뻹�Ǳ�ƽ�Ƶ�
numE = 0;
numL = 0;
count = 34;                                                                %�ӵ�35���㿪ʼ����
while numE ~= (Pnum + 34 + numL) %&& count<= (((r-2)*(c-2)-nnz(mask))/2-1)     %�� Ƕ����Ϣ������=Pnum+34+��λͼ���� ʱ����
    count = count + 1;
    [mdNum,EorS] = modifiedNum(img_pre(seqToRC(count,1),seqToRC(count,2)),err_sorted(count),Tn,Tp);%���㵱ǰ����޸Ĵ������Լ����Ǳ�Ƕ�뻹�Ǳ�ƽ�Ƶ�
    if mdNum == 1                                                          %�����޸�һ�Σ�����ʵ�ʵ�Ƕ������ƽ�Ƶ㣩
        if EorS == 1                                                       %����Ƕ���
            numE = numE + 1;                                               %Ƕ����Ϣ��������1 
            msg(numE) = bitget(abs(err_sorted(count)),1);                  %��ȡǶ����Ϣ
            err_sorted(count) = floor(err_sorted(count)/2);                %�ָ�ԭԤ���ֵ
        elseif err_sorted(count) > -1                                      %����ƽ�Ƶ㣬��ָ�Ԥ���ֵ
            err_sorted(count) = err_sorted(count) - Tp - 1;
        else
            err_sorted(count) = err_sorted(count) - Tn;
        end
    else                                                                   %�������޸ģ����¼��L��
        numL = numL + 1;
        L(1,numL) = count;
        L(2,numL) = EorS;
    end
end

LM = msg(1,1:numL);                                                        %����ȡ������Ϣ���õ�λͼ��pl���Լ�lsb
pl = msg(1,numL+1:numE - 34);
lsb34 = msg(1,numE - 33:numE);

for i = 1:numL                                                             %����Ҫλͼ�ĵ���в���
    if LM(i)                                                               %�����޸Ĺ��ĵ㣬��ָ�ԭ��ֵ
        if L(2,i) == 1                                                     
            err_sorted(L(1,i)) = floor(err_sorted(L(1,i))/2);
        elseif err_sorted(L(1,i)) > -1
            err_sorted(L(1,i)) = err_sorted(L(1,i)) - Tp - 1;
        else
            err_sorted(L(1,i)) = err_sorted(L(1,i)) - Tn;
        end
    end
end

for i = 1:34                                                               %�ָ�ǰ34�����lsb
    xf(seqToRC(i,1),seqToRC(i,2)) = bitset(xf(seqToRC(i,1),seqToRC(i,2)),1,lsb34(i));
end
for i = 35:count                                                           %ʵ�ʻָ�ͼ��
    row = seqToRC(i,1);
    col = seqToRC(i,2);
    xf(row,col) = img_pre(row,col) + err_sorted(i);
end


%%
function [x] = deBitSeq(seq)
%�Զ��������н��н��룬�õ�������
%����:����������seq
%�����������x
len = numel(seq);
x = uint32(0);
for i=1:len
    x = bitset(x,i,seq(i));
end
%%
function [num,EorS] = modifiedNum(xpre,err,Tn,Tp)
%����һ����Ŀ��޸Ĵ������Լ�Ƕ��ʱ���޸���Ƕ�뻹��ƽ��
%���룺Ԥ��ֵxpre��Ԥ���ֵerr����ֵTn��Tp
%��������޸Ĵ���num����ʱ���޸���Ƕ�뻹��ƽ��ָʾEorS��1����Ƕ�룩
num = 0;
EorS = 1;
if err > Tp*2 + 1 || err < Tn*2
    EorS = 2;
end
if err > -1 && err <= Tp
    if xpre + err*2 + 1 <= 255
        num =1;
    end
elseif err < 0 && err >= Tn
    if xpre + err*2 >= 0
        num = 1;
    end
elseif err > Tp
    if xpre + err + Tp + 1 <= 255
        num = 1;
    end
else
    if xpre + err + Tn >= 0 
        num = 1;
    end
end