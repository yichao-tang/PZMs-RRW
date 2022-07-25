function [xf,pl] = reversible_decodding(f, Rows, Cols, mask)
%使用vasiliy的方法进行解码
%输入：带水印图像f
%输出：原图像xf，信息pl
xf = double(f);
[xf,pl2] = decodeCrossOrDot(xf,2, Rows, Cols, mask);                                         %对O集合进行解码
[xf,pl1] = decodeCrossOrDot(xf,1, Rows, Cols, mask);                                         %对X集合进行解码
pl =[pl1,pl2];
%%
function [xf,pl] = decodeCrossOrDot(f,CorD, r, c, mask)
%对X或者O集合进行解码
%输入：带水印图像f，XO集指示CorD（1代表X集合）
%输出：去除X或O集合水印的图像xf,信息pl
xf = f;
% [r,c] = size(xf);

[img_pre,img_err,err_sorted,seq_map] = sortErrWithSur(xf,CorD);            %获得X或O集预测值，预测差值，排序后的预测差值，预测差值的序号

seqToRC = zeros( (r-2)*(c-2)/2 , 2);                                       %由像素点排序后的序号得到其坐标
tmp = r*c;
for i = 2:r-1
    for j = 2:c-1
        if seq_map(i,j) ~= tmp && mask(i,j) == 0 %此处添加判定条件，只取内切圆外测像素点
            seqToRC( seq_map(i,j) , : ) = [i,j];            
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
Tn = deBitSeq(lsb34(1,1:7));                                               %由lsb解码出阈值Tn，Tp，Pnum
Tn = -int32(Tn);
Tp = int32(deBitSeq(lsb34(1,8:14)));
Pnum = deBitSeq(lsb34(1,15:34));

msg = zeros(1,(((r-2)*(c-2))-nnz(mask))/2);                                                      %保存解码出的信息 位图+pl+lsb
L = zeros(2,1000000);                                                      %记录需要位图的点，以及其是被嵌入还是被平移的
numE = 0;
numL = 0;
count = 34;                                                                %从第35个点开始搜索
while numE ~= (Pnum + 34 + numL) %&& count<= (((r-2)*(c-2)-nnz(mask))/2-1)     %当 嵌入信息点数！=Pnum+34+需位图点数 时继续
    count = count + 1;
    [mdNum,EorS] = modifiedNum(img_pre(seqToRC(count,1),seqToRC(count,2)),err_sorted(count),Tn,Tp);%计算当前点可修改次数，以及其是被嵌入还是被平移的
    if mdNum == 1                                                          %若可修改一次（即是实际的嵌入点或者平移点）
        if EorS == 1                                                       %若是嵌入点
            numE = numE + 1;                                               %嵌入信息点数自增1 
            msg(numE) = bitget(abs(err_sorted(count)),1);                  %提取嵌入信息
            err_sorted(count) = floor(err_sorted(count)/2);                %恢复原预测差值
        elseif err_sorted(count) > -1                                      %若是平移点，则恢复预测差值
            err_sorted(count) = err_sorted(count) - Tp - 1;
        else
            err_sorted(count) = err_sorted(count) - Tn;
        end
    else                                                                   %若不能修改，则记录到L中
        numL = numL + 1;
        L(1,numL) = count;
        L(2,numL) = EorS;
    end
end

LM = msg(1,1:numL);                                                        %由提取出的信息，得到位图，pl，以及lsb
pl = msg(1,numL+1:numE - 34);
lsb34 = msg(1,numE - 33:numE);

for i = 1:numL                                                             %对需要位图的点进行操作
    if LM(i)                                                               %若是修改过的点，则恢复原差值
        if L(2,i) == 1                                                     
            err_sorted(L(1,i)) = floor(err_sorted(L(1,i))/2);
        elseif err_sorted(L(1,i)) > -1
            err_sorted(L(1,i)) = err_sorted(L(1,i)) - Tp - 1;
        else
            err_sorted(L(1,i)) = err_sorted(L(1,i)) - Tn;
        end
    end
end

for i = 1:34                                                               %恢复前34个点的lsb
    xf(seqToRC(i,1),seqToRC(i,2)) = bitset(xf(seqToRC(i,1),seqToRC(i,2)),1,lsb34(i));
end
for i = 35:count                                                           %实际恢复图像，
    row = seqToRC(i,1);
    col = seqToRC(i,2);
    xf(row,col) = img_pre(row,col) + err_sorted(i);
end


%%
function [x] = deBitSeq(seq)
%对二进制序列进行解码，得到正整数
%输入:二进制序列seq
%输出：正整数x
len = numel(seq);
x = uint32(0);
for i=1:len
    x = bitset(x,i,seq(i));
end
%%
function [num,EorS] = modifiedNum(xpre,err,Tn,Tp)
%计算一个点的可修改次数，以及嵌入时的修改是嵌入还是平移
%输入：预测值xpre，预测差值err，阈值Tn，Tp
%输出：可修改次数num，入时的修改是嵌入还是平移指示EorS（1代表嵌入）
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