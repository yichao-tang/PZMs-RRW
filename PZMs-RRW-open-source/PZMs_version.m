function [ psnr_Iw2,psnr_Iw3, BER_no_attack ]= PZMs_version(I, nMax, Delta, num, T_start,slide)
img=double(I);
rng(0,'twister');

%% Initialize
% nMax =26; %256bits
% nMax =18; %128bits

[Rows, Cols]=size(img);
xstep = 2/(Rows-1);
ystep = 2/(Cols-1);
[X,Y] = meshgrid(-1:xstep:1, 1:-ystep:-1);
[theta,rho] = cart2pol(X, Y);
inside = find(rho<=1);
mask = zeros(Rows,Cols);
mask(inside) = ones(size(inside));
moment_area = nnz(mask);

CenX = (Rows / 2);
CenY = (Cols / 2);

rmin = Rows * Cols;
radial = CenY;
rmin = min(radial, rmin);
radial = CenX;
rmin = min(radial, rmin);
radial = Cols - CenX;
rmin = min(radial, rmin);
radial = Rows - CenY;
rmin = min(radial, rmin);
RMax = rmin;

if Rows == Cols
    DMax = Rows;
else
    DMax = 2*RMax;
end

%% 
temp_i=1;
map(temp_i)=0; % n
map_2(temp_i)=0; % m
for n=1 : nMax
    for m = n : -1 : 0
        temp_i=temp_i+1;
        map(temp_i)=n;
        map_2(temp_i)=m;
    end
end
for k = 2 : length(map)
    K_3(k) = ((-2)*(2*map_2(k)+3)*(map_2(k)+1))/((map(k)+map_2(k)+2)*(map(k)-map_2(k)));
    K_2(k) = ((K_3(k)*(map(k)+map_2(k)+3)*(map(k)-map_2(k)-1))/(2*(map_2(k)+2)))+(2*map_2(k)+3);
    K_1(k) = ((2*map_2(k)+5)*(map_2(k)+2))-((2*map_2(k)+5)*K_2(k))+(K_3(k)*(map(k)+map_2(k)+4)*(map(k)-map_2(k)-2))/2;
end

%% PZMs of original image 
tic;
Znm = zeros(1,length(map));
for i = 1 : DMax
    for j = 1 : DMax
        if (mask(i, j) == 0)
            continue
        end
        if (mask(i, j) == 1)
            r(1) = 1;
            r(3) = rho(i,j)^2;
            r(2) = rho(i,j);
            for rh = 4 : nMax+1
                r(rh) = r(rh-2)*r(3);
            end
            cos_C(1) = 1;
            sin_C(1) = 0;
            cos_C(2) = cos(theta(i,j));
            sin_C(2) = sin(theta(i,j));
            tan_C = tan(theta(i,j));
            cos_C(3) = (1-tan_C^2)/(1+tan_C^2);
            sin_C(3) = (2*tan_C)/(1+tan_C^2);
            for th = 4 : nMax+1
                cos_C(th) = cos_C(th-2)*cos_C(3) - sin_C(th-2)*sin_C(3);
                sin_C(th) = sin_C(th-2)*cos_C(3) + sin_C(3)*cos_C(th-2);
            end
            for k = 1 : length(map)
                if map_2(k) == map(k)
                    Rad(k) = r(map(k)+1);
                elseif map_2(k) == map(k) - 1
                    Rad(k) = (2*map(k)+1)*r(map(k)+1) - 2*map(k)*r(map(k));
                else
                    Rad(k) = K_1(k)*Rad(k-2) + (K_2(k)+(K_3(k)/r(2)))*Rad(k-1);
                end
                Znm(k) = Znm(k) + ((map(k)+1)/(moment_area))...
                    * Rad(k) * complex(img(i,j)) * ...
                    (cos_C(map_2(k)+1)-1i*sin_C(map_2(k)+1));
            end
        end
    end
end
t1 = toc;

%% Robust Embedding Stage
tic;
Irw1 = zeros(Rows, Cols);
Irw2 = zeros(Rows, Cols);
Vnm = zeros(Rows, Cols);
Rnm = zeros(Rows, Cols);
Rnm_1 = zeros(Rows, Cols);
Rnm_2 = zeros(Rows, Cols);

moment = 0;
w =randi([0,1],1,num);
d_0 = 7/8*Delta*ones(1,num);
d_1 = d_0 + (Delta/2);
for k = 1 : length(map)
    if mod(map_2(k), 4) ~= 0
        moment = moment + 1;
        tem(moment)=k;
    end
end
insert_palace = tem(1:num);

% adaptive normalized weight
for k = 1 : length(map)
    T(k) = T_start-map(k)*slide;
end 

moment = 0;
for k = 1 : length(map)
        if ismember(k,insert_palace)
            moment = moment + 1;
            ZRnm(k)=(Znm(k)/Znm(1))*T(k);
            ZRnm_j(1,k)=round((abs(ZRnm(k))-d_0(moment))/Delta)*Delta+d_0(moment);
            ZRnm_j(2,k)=round((abs(ZRnm(k))-d_1(moment))/Delta)*Delta+d_1(moment);
            d_int(moment)  = round(abs(ZRnm_j(w(moment)+1,k))-abs(ZRnm(k)));
            dec_R(moment) = abs(ZRnm_j(w(moment)+1,k))-abs(ZRnm(k)) - d_int(moment);
            ZRw(k) = abs(ZRnm_j(w(moment)+1,k)) - dec_R(moment);
            Zw(k) = (abs(ZRw(k))/abs(ZRnm(k)))*Znm(k);
        end
    for i = 1 : Rows
        for j = 1 : Cols
            if (mask(i, j) == 0)
                Vnm(i, j) = NaN;
            end
            if (mask(i, j) == 1)
                r(1) = 1;
                r(3) = rho(i,j)^2;
                r(2) = rho(i,j);
                for rh = 4 : nMax+1
                    r(rh) = r(rh-2)*r(3);
                end
                if map_2(k) == map(k)
                    Rnm(i,j) = r(map(k)+1);
                    if k ==1
                        Rnm_1(i,j) = Rnm(i,j);
                        Rnm_2(i,j) = Rnm(i,j);
                    else
                        Rnm_2(i,j)=Rnm_1(i,j);
                        Rnm_1(i,j) = Rnm(i,j);
                    end
                elseif map_2(k) == map(k) - 1
                    Rnm(i,j) = (2*map(k)+1)*r(map(k)+1) - 2*map(k)*r(map(k));
                    Rnm_2(i,j)=Rnm_1(i,j);
                    Rnm_1(i,j) = Rnm(i,j);
                else
                    Rnm(i,j) = K_1(k)*Rnm_2(i,j) + (K_2(k)+(K_3(k)/r(2)))*Rnm_1(i,j);
                    Rnm_2(i,j)=Rnm_1(i,j);
                    Rnm_1(i,j) = Rnm(i,j);
                end
                Vnm(i,j) = Rnm(i,j).*exp(1i*map_2(k)*theta(i, j));
            end
        end
    end
        if ismember(k,insert_palace)
            Irw1 = Irw1 + (Zw(k)-Znm(k))*Vnm...
                + (conj(Zw(k))-conj(Znm(k)))*conj(Vnm);
        end
end

Irw1(isnan(Irw1))=0;
Iw1 = round(img + Irw1);
Iw1(Iw1>255)=255;
Iw1(Iw1<0)=0;
psnr_Iw1 = psnr(uint8(img),uint8(Iw1));

%% rounded error compensation
Znm_w = zeros(1,length(map));
for i = 1 : DMax
    for j = 1 : DMax
        if (mask(i, j) == 0)
            continue
        end
        if (mask(i, j) == 1)
            r(1) = 1;
            r(3) = rho(i,j)^2;
            r(2) = rho(i,j);
            for rh = 4 : nMax+1
                r(rh) = r(rh-2)*r(3);
            end
            cos_C(1) = 1;
            sin_C(1) = 0;
            cos_C(2) = cos(theta(i,j));
            sin_C(2) = sin(theta(i,j));
            tan_C = tan(theta(i,j));
            cos_C(3) = (1-tan_C^2)/(1+tan_C^2);
            sin_C(3) = (2*tan_C)/(1+tan_C^2);
            for th = 4 : nMax+1
                cos_C(th) = cos_C(th-2)*cos_C(3) - sin_C(th-2)*sin_C(3);
                sin_C(th) = sin_C(th-2)*cos_C(3) + sin_C(3)*cos_C(th-2);
            end
            for k = 1 : length(map)
                if map_2(k) == map(k)
                    Rad(k) = r(map(k)+1);
                elseif map_2(k) == map(k) - 1
                    Rad(k) = (2*map(k)+1)*r(map(k)+1) - 2*map(k)*r(map(k));
                else
                    Rad(k) = K_1(k)*Rad(k-2) + (K_2(k)+(K_3(k)/r(2)))*Rad(k-1);
                end
                Znm_w(k) = Znm_w(k) + ((map(k)+1)/(moment_area))...
                    * Rad(k) * complex(Iw1(i,j)) * ...
                    (cos_C(map_2(k)+1)-1i*sin_C(map_2(k)+1));
            end
        end
    end
end
for k = 1 : length(map)
    for i = 1 : Rows
        for j = 1 : Cols
            if (mask(i, j) == 0)
                Vnm(i, j) = NaN;
            end
            if (mask(i, j) == 1)
                r(1) = 1;
                r(3) = rho(i,j)^2;
                r(2) = rho(i,j);
                for rh = 4 : nMax+1
                    r(rh) = r(rh-2)*r(3);
                end
                if map_2(k) == map(k)
                    Rnm(i,j) = r(map(k)+1);
                    if k ==1
                        Rnm_1(i,j) = Rnm(i,j);
                        Rnm_2(i,j) = Rnm(i,j);
                    else
                        Rnm_2(i,j)=Rnm_1(i,j);
                        Rnm_1(i,j) = Rnm(i,j);
                    end
                elseif map_2(k) == map(k) - 1
                    Rnm(i,j) = (2*map(k)+1)*r(map(k)+1) - 2*map(k)*r(map(k));
                    Rnm_2(i,j)=Rnm_1(i,j);
                    Rnm_1(i,j) = Rnm(i,j);
                else
                    Rnm(i,j) = K_1(k)*Rnm_2(i,j) + (K_2(k)+(K_3(k)/r(2)))*Rnm_1(i,j);
                    Rnm_2(i,j)=Rnm_1(i,j);
                    Rnm_1(i,j) = Rnm(i,j);
                end
                Vnm(i,j) = Rnm(i,j).*exp(1i*map_2(k)*theta(i, j));
            end
        end
    end
    if ismember(k,insert_palace)
        Irw2 = Irw2 + (Zw(k)-Znm_w(k)).*Vnm...
            + (conj(Zw(k))-conj(Znm_w(k)))*conj(Vnm);
    end
end
Irw2(isnan(Irw2))=0;
Iw2 = round(Iw1 + Irw2);
Iw2(Iw2>255)=255;
Iw2(Iw2<0)=0;
psnr_Iw2 = psnr(uint8(img),uint8(Iw2))
% imwrite(uint8(Iw2),sprintf("Inserted_image.bmp"))
%% PZMs of Iw 
Znm_w2 = zeros(1,length(map));
for i = 1 : DMax
    for j = 1 : DMax
        if (mask(i, j) == 0)
            continue
        end
        if (mask(i, j) == 1)
            r(1) = 1;
            r(3) = rho(i,j)^2;
            r(2) = rho(i,j);
            for rh = 4 : nMax+1
                r(rh) = r(rh-2)*r(3);
            end
            cos_C(1) = 1;
            sin_C(1) = 0;
            cos_C(2) = cos(theta(i,j));
            sin_C(2) = sin(theta(i,j));
            tan_C = tan(theta(i,j));
            cos_C(3) = (1-tan_C^2)/(1+tan_C^2);
            sin_C(3) = (2*tan_C)/(1+tan_C^2);
            for th = 4 : nMax+1
                cos_C(th) = cos_C(th-2)*cos_C(3) - sin_C(th-2)*sin_C(3);
                sin_C(th) = sin_C(th-2)*cos_C(3) + sin_C(3)*cos_C(th-2);
            end
            for k = 1 : length(map)
                if map_2(k) == map(k)
                    Rad(k) = r(map(k)+1);
                elseif map_2(k) == map(k) - 1
                    Rad(k) = (2*map(k)+1)*r(map(k)+1) - 2*map(k)*r(map(k));
                else
                    Rad(k) = K_1(k)*Rad(k-2) + (K_2(k)+(K_3(k)/r(2)))*Rad(k-1);
                end
                Znm_w2(k) = Znm_w2(k) + ((map(k)+1)/(moment_area))...
                    * Rad(k) * complex(Iw2(i,j)) * ...
                    (cos_C(map_2(k)+1)-1i*sin_C(map_2(k)+1));
            end
        end
    end
end

%% Watermark-Removed Image Generation
moment = 0;
for k = 1 : length(map)
        if ismember(k,insert_palace)
            moment = moment + 1;
            ZRnm_w(k) = (Znm_w2(k)/Znm_w2(1))*T(k);
            ZRnm_w_j(1,k)=round((abs(ZRnm_w(k))-d_0(moment))/Delta)*Delta+d_0(moment);
            ZRnm_w_j(2,k)=round((abs(ZRnm_w(k))-d_1(moment))/Delta)*Delta+d_1(moment);
            if (abs(ZRnm_w(k))-abs(ZRnm_w_j(1,k)))^2 <= (abs(ZRnm_w(k))-abs(ZRnm_w_j(2,k)))^2
                w_ex(moment) = 0;
            else
                w_ex(moment) = 1;
            end
            w_be(moment)=w(moment);
            ZRnm_re(k) = abs(ZRnm_w(k))-d_int(moment);
            Znm_re(k) = abs(ZRnm_re(k)/abs(ZRnm_w(k)))*Znm_w2(k);
        end
end 

%
isequal(w_ex,w_be)
te=sum(abs(w_be-w_ex));
BER_no_attack = te/length(w_ex);

%
Irw_re=zeros(size(I));
for k = 1 : length(map)
    for i = 1 : Rows
        for j = 1 : Cols
            if (mask(i, j) == 0)
                Vnm(i, j) = NaN;
            end
            if (mask(i, j) == 1)
                r(1) = 1;
                r(3) = rho(i,j)^2;
                r(2) = rho(i,j);
                for rh = 4 : nMax+1
                    r(rh) = r(rh-2)*r(3);
                end
                if map_2(k) == map(k)
                    Rnm(i,j) = r(map(k)+1);%rho(i,j)^map(k);%
                    if k ==1
                        Rnm_1(i,j) = Rnm(i,j);
                        Rnm_2(i,j) = Rnm(i,j);
                    else
                        Rnm_2(i,j)=Rnm_1(i,j);
                        Rnm_1(i,j) = Rnm(i,j);
                    end
                elseif map_2(k) == map(k) - 1
                    Rnm(i,j) = (2*map(k)+1)*r(map(k)+1) - 2*map(k)*r(map(k));%(2*map(k)+1)*rho(i,j)^map(k) - 2*map(k)*rho(i,j)^(map(k)-1);%
                    Rnm_2(i,j)=Rnm_1(i,j);
                    Rnm_1(i,j) = Rnm(i,j);
                else
                    Rnm(i,j) = K_1(k)*Rnm_2(i,j) + (K_2(k)+(K_3(k)/r(2)))*Rnm_1(i,j);%K_1(k)*Rnm_2(i,j) + (K_2(k)+(K_3(k)/rho(i,j)^2))*Rnm_1(i,j);%
                    Rnm_2(i,j)=Rnm_1(i,j);
                    Rnm_1(i,j) = Rnm(i,j);
                end
                Vnm(i,j) = Rnm(i,j).*exp(1i*map_2(k)*theta(i, j));%(cos_C(map_2(k)+1)-1i*sin_C(map_2(k)+1));
            end
        end
    end
        if ismember(k,insert_palace)
            Irw_re = Irw_re + (Znm_re(k)-Znm_w2(k)).*Vnm...
                + (conj(Znm_re(k))-conj(Znm_w2(k)))*conj(Vnm);
        end
end
Irw_re(isnan(Irw_re))=0;
I_re = round(Irw_re + Iw2);
I_re(I_re>255)=255;
I_re(I_re<0)=0;
%% Computing dr
temp = 0;
for i = 1 : DMax
    for j = 1 : DMax
        if (mask(i, j) == 1)
            temp = temp + 1;
            err_incircle(temp) = I_re(i,j)-img(i,j);
            continue
        end
    end
end
%% Encode
tem_length1 = length(dec2bin(max(abs(err_incircle))));
[ dq_1 ] = signed_to_bin(err_incircle,tem_length1+1);
dq_C = cell(1,1);
dq_C{1} = dq_1;
dq_data = Arith07(dq_C);
[ dq_D ] = unsigned_to_bin( dq_data,8 );
dq_beta=dq_D;

%% Decode
% data_2=dq_beta;
% [ data ] = bin_to_unsigned( data_2,8);
% xR = Arith07(data');
% d1_2 = xR{1};
% [ d1 ] =  bin_to_signed( d1_2(:)',tem_length1+1);

%% Reversible Embedding Stage
err = dq_beta;
err=logical(err);
[Iw_reversible, ~, Ok, ~, ~, ~, ~, ~, ~] ...
        = reversible_embedding(Iw2, err, Rows, Cols, mask) ;
psnr_Iw3 = psnr(uint8(img),uint8(Iw_reversible))
% imwrite(uint8(Iw_reversible),'Iw_reversible.bmp');

%%
%% Reversible Extraction
[Iw_1,err_info] = reversible_decodding(Iw_reversible, Rows, Cols, mask);
check_img=isequal(Iw_1,Iw2); %chect image extraction
check_err= isequal(err_info,err); %%chect watermarking extraction


end
%% 原版伪zernike基计算
function [Vnm, Rnm] =  CmpVnm(n, m, rho, theta)
% Rnm: 径向多项式
% Vnm: 正交多项式（zernike基）
% rho: 半径ρ
% theta: 像素点在单位圆上的夹角
% n: 阶数
% m: 重复数
m1 = abs(m);
nCtrl = n-m1;
Rnm = zeros(size(rho));
for s = 0 : nCtrl
    dCoeff = factorial(s)*factorial(n+m1+1-s)*...
        factorial(n-m1-s);
    dCoeff = factorial(2*n+1-s)/dCoeff;
    Rnm = Rnm + (-1)^s*dCoeff*rho.^(n-s); % Rnm(ρ)
end
Vnm = Rnm.*exp(1i*m*theta);
end