function [xf,PSNR,Ok,per1,per2,P1,Z1,P2,Z2] = reversible_embedding(f, reversible_watermarking, Rows, Cols, mask) 
%使用vasiliy的方法镜像嵌入
%输入：图像f，可逆水印信息pl
%输出：带水印图像xf，PSNR，是否全部嵌入Ok，per1和per2分别表示与X集和O集的嵌入率，P1、P2以及Z1、Z2暂时意义不明
% xf = double(f);
% pl = rand(1,rate) > 0.5;
rate = numel(reversible_watermarking);
[xf,okC,per1,P1,Z1] = embedCrossOrDot...
    (f,reversible_watermarking(1,1:round(rate/2)),1, Rows, Cols, mask);                    %X集嵌入
[xf,okD,per2,P2,Z2] = embedCrossOrDot...
    (xf,reversible_watermarking(1,round(rate/2)+1:numel(reversible_watermarking)),2, Rows, Cols, mask);          %O集嵌入
if okD && okC
    Ok = 1;
else
    Ok = 0;
end
[PSNR,mse] = psnr(f,xf);
%%
function [xf,fullEmbed,per,P,Z] = embedCrossOrDot(f,pl,CorD, r, c, mask)
%对X集或者O集进行嵌入
%输入：图像f，需嵌入信息pl，XO集标志CorD
%输出：带水印图像xf，受否全部嵌入fullEmbed
xf = f;
per = 0;
% [r,c] = size(xf);
Pnum = numel(pl);
Pnum = Pnum + 34;                                                          %实际嵌入pl和34个lsb
fullEmbed = 1;

[img_pre,img_err,err_sorted,seq_map] = sortErrWithSur(xf,CorD);            %获得X或O集预测值，预测差值，排序后的预测差值，预测差值的序号
% [img_pre,img_err,err_sorted,seq_map,Tu] = impRhom(xf,CorD);
seqToRC = zeros( (r-2)*(c-2)/2 , 2);                                       %由像素点排序后的序号得到其坐标
tmp = r*c;

if Pnum>28000
    for i = 2:r-1
        for j = 2:c-1
            if seq_map(i,j) ~= tmp %&& mask(i,j) == 0 %此处添加判定条件，只取内切圆外测像素点 %当超过外切圆能容纳的值之后就无法正常运行了，暂时注释
                seqToRC( seq_map(i,j) , : ) = [i,j];
            end
        end
    end
else
    for i = 2:r-1
        for j = 2:c-1
            if seq_map(i,j) ~= tmp && mask(i,j) == 0 %此处添加判定条件，只取内切圆外测像素点 %当超过外切圆能容纳的值之后就无法正常运行了，暂时注释
                seqToRC( seq_map(i,j) , : ) = [i,j];
            end
        end
    end
end
id = find(seqToRC(:,1)==0 & seqToRC(:,2)==0);% 删除内切圆中的行
seqToRC(id,:)=[];
err_sorted(:,id)=[];

lsb34 = zeros(1,34);                                                       %获取排序后前34个像素的lsb
for i = 1:34
    lsb34(i) = bitget(xf(seqToRC(i,1),seqToRC(i,2)),1);
end

Tn = int16(0);                                                             %初始化Tn和Tp
Tp = int16(0);
H = hist(err_sorted(1,35:numel(err_sorted)),-255:255);
while sum(H(1,Tn+256:Tp+256)) < Pnum
    if abs(Tn) <= Tp
        Tn = Tn -1;
    else
        Tp = Tp + 1;
    end
end

L = zeros(3,10000);                                                        %记录需要位图的点的序号，可修改次数，以及是嵌入位还是平移位
E = zeros(2,numel(err_sorted) - 34);                                       %记录可修改两次的点的序号，以及是嵌入位还是平移位
numE = 0;
numES = 0;
numL = 0;
while numE < Pnum + numL                                                   %当可用于嵌入的点数小于嵌入信息大小+位图大小时继续
    for i = 35:numel(err_sorted)                                           %对从序号为35的点，开始操作
        [mdNum,EorS] = modifiedNum(img_pre(seqToRC(i,1),seqToRC(i,2)),err_sorted(i),Tn,Tp);%计算当前点能够修改几次，是处于嵌入位还是平移位
        if mdNum == 2                                                      %若可修改两次
            numES = numES + 1;                                             %则将其记录到E中
            E(1,numES) = i;
            E(2,numES) = EorS;
            if EorS == 1                                                   %若且是可嵌入位，则numE自增1
                numE = numE + 1;
            end
        else                                                               %若可修改一次或不能修改
            numL = numL + 1;                                               %则将其记录到L中
            L(1,numL) = i;                                                 
            L(2,numL) = mdNum;
            L(3,numL) = EorS;
        end
        if Pnum + numL == numE
            per = double(i)/numel(err_sorted);
            break;
        end
    end
    if Pnum + numL ~= numE                                                 %若容量不足，则增大Tp或者abs(Tn)
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

msg = [L(2,1:numL),pl,lsb34];                                              %需要嵌入的全部信息 = 位图+pl+lsb
count = 0;
for i = 1:numES                                                            %对可修改两次的点，进行嵌入和平移
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
    xf(row,col) = img_pre(row,col) + err;                                  %实际修改图像
end
if count ~= numel(msg)
    fullEmbed = 0;
end
for i = 1:numL                                                             %对需要记录位图的点
    err = err_sorted(L(1,i));
    if L(2,i) == 1                                                         %若可修改一次，则对其进行一次幅度最大的修改
        if L(3,i) == 1                                                     %对嵌入位
            if err < 0
                err = err*2;
            else
                err = err*2 + 1;
            end
        else                                                               %对平移位
            if err < 0 
                err = err + Tn;
            else
                err = err + Tp + 1;
            end
        end
    end
    row = seqToRC(L(1,i),1);
    col = seqToRC(L(1,i),2);
    xf(row,col) = img_pre(row,col) + err;                                  %实际修改图像
end
lsb34 = [toBitSeq(uint16(abs(Tn)),7),toBitSeq(uint16(Tp),7),toBitSeq(uint32(Pnum - 34),20)];%用于替换的lsb = abs（Tn） + Tp + Pnum - 34
for i = 1:34                                                               %lsb替换，修改图像
    xf(seqToRC(i,1),seqToRC(i,2)) = bitset(xf(seqToRC(i,1),seqToRC(i,2)),1,lsb34(i));
end
%%
function [BS]  = toBitSeq(x,l)
%将输入的整数变为二进制序列
%输入：整数x，二进制序列长度l
%输出：二进制序列BS
BS = zeros(1,l);
x = abs(x);
for i=1:l
    BS(1,i) = bitget(x,i);
end
%%
function [num,EorS]  = modifiedNum(xpre,err,Tn,Tp)
%计算一个点可以修改的次数，以及第一次修改是嵌入还是平移
%输入：预测值xpre，预测差值err，阈值Tn，Tp
%输出：可修改次数num，嵌入平移标志EorS（1代表嵌入，2代表平移）
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
