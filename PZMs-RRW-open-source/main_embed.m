clc;clear;
%% Parameters setting
nMax=18;% nMax =18;（128bits  % nMax =26; （256bits
num=128;
Delta=32;
slide=10;

%% Images Reading
file_path =  'image\';
img_path_list = dir(strcat(file_path,'*.bmp'));
img_num = length(img_path_list);

%% Watermark Embedding
if img_num > 0
    temp=0;
    for j = 1:4
        image_name = img_path_list(j).name;
        image =  imread(strcat(file_path,image_name));
        mysize=size(image);
        if numel(mysize)>2
            if mysize(3) ==2
                image = image(:,:,1);
            else
                image=rgb2gray(image);
            end
        end
        [image_Rows, image_Cols]=size(image);
        if image_Rows~=512 || image_Cols~=512
            image =imresize(uint8(image),[512,512]);
        end
        for T_start =2000 %128bit %2400 %256bit
                temp=temp+1;
                [psnr1(temp,1),psnr2(temp,1), BER_no_attack(temp,1)]...
                    = PZMs_version(image, nMax, Delta, num, T_start, slide);
                Card(temp,:) = [j;seed_num];
        end
    end
end
