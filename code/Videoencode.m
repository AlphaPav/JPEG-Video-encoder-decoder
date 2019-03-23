function Videoencode()
    clear all;
    video_path='mayun.avi';    
    video_obj=VideoReader(video_path);   % read the video and save it as video_obj
    global number; %the global variable 'number' is used to name every txt document
    number = 1;
    frame_number=video_obj.NumberOfFrames; % get the frame number of the video 
    frame_rate=video_obj.FrameRate; %get the frame rate of the video
    
    new_folder1 = 'Figure'; 
    new_folder2 = 'NewFigure';
    mkdir(new_folder1); %generate a new folder named 'Figure'
    mkdir(new_folder2); %generate a new folder named 'NewFigure'
    
    for i=1:frame_number
        oriframe=read(video_obj,i); % read the frame one by one from the video 
        cd(strcat('Figure'));
        myfilename = strcat('pic_',num2str(i),'.txt');
        save(myfilename, 'oriframe'); %save the original frame as a ".txt"
        cd ..
        Imagecompress(oriframe); % compress and encode the current frame
    end
    long=size(oriframe,1);% get the size of the frame and save them as long and height
    height=size(oriframe,2);
    save('VideoAttribute.txt','long','height','frame_number','frame_rate');
    % save video's attributes in order to recover the video in the decoding preocess
end



function Imagecompress(img)
global number;

%define quantization matrix to Normalize the DCT Matrix. Norm_Mat1 is used
%for Y layer while Norm_Mat2 is used for U and V layer.
Norm_Mat1=[16 11 10 16 24 40 51 61      
          12 12 14 19 26 58 60 55
          14 13 16 24 40 57 69 56
          14 17 22 29 51 87 80 62
          18 22 37 56 68 109 103 77
          24 35 55 64 81 104 113 92
          49 64 78 87 103 121 120 101
          72 92 95 98 112 100 103 99];
Norm_Mat2=[17 18 24 47 99 99 99 99;  
         18 21 26 66 99 99 99 99;  
         24 26 56 99 99 99 99 99;  
         47 66 99 99 99 99 99 99;  
         99 99 99 99 99 99 99 99;  
         99 99 99 99 99 99 99 99;  
         99 99 99 99 99 99 99 99;  
         99 99 99 99 99 99 99 99;]; 
img=double(img);
RGB = img;
R = RGB(:,:,1);  
G = RGB(:,:,2);  
B = RGB(:,:,3);
% RGB to YUV Conversion
Y=0.299*double(R)+0.587*double(G)+0.114*double(B); % Y is the luma component which represents the intensity of the image and look likes a gray scale version.
U=-0.169*double(R)-0.3316*double(G)+0.5*double(B); % U and V are the chrominance components which represent the color information in the image
V=0.5*double(R)-0.4186*double(G)-0.0813*double(B);

% Chrominance Down-Samlping
[m,n]=size(U);
for i =1:m
    for j=1:n/2
        U(i,2*j)=U(i,2*j-1); % for the U image layer, make the adjacent pixels share the same chromatic value
        V(i,2*j)=V(i,2*j-1); % for the V image layer, make the adjacent pixels share the same chromatic value
    end
end

% Compress the Y,U,V matrix
[huff_DCY,dict_DCY,huff_ACY,dict_ACY]=Compression(Y,Norm_Mat1);
[huff_DCU,dict_DCU,huff_ACU,dict_ACU]=Compression(U,Norm_Mat2);
[huff_DCV,dict_DCV,huff_ACV,dict_ACV]=Compression(V,Norm_Mat2);

long=size(img,1);
height=size(img,2);
% Save the encoding results to a 'txt'
cd(strcat('NewFigure'));
myfilename = strcat('pic_',num2str(number),'.txt');
save(myfilename ,'huff_DCY','dict_DCY','huff_ACY','dict_ACY','huff_DCU','dict_DCU','huff_ACU','dict_ACU','huff_DCV','dict_DCV','huff_ACV','dict_ACV');
cd ..
number=number+1;


end

%----------------------------------------------------------
% Image Compression Overview
%----------------------------------------------------------
function [huff_DC,dict_DC,huff_AC,dict_AC]=Compression(img,Norm_Mat)
quanti_img=DCT(img,Norm_Mat);
DPCM=DPCMencode(quanti_img);
[huff_DC,dict_DC]= HUFFMANencodeDC(DPCM);
zig_seq=ZIGZAGscan(quanti_img);
ac_rle=RunlengthEncode(zig_seq);
[huff_AC,dict_AC]= HUFFMANencodeAC(ac_rle);
end

%----------------------------------------------------------
% DCT 
%----------------------------------------------------------
 function quanti_img=DCT(img,Norm_Mat)
 T = dctmtx(8); %  create a 8*8 two-dimension DCT conversion matrix
 dct_img = blkproc(img,[8 8],'P1*x*P2',T,T');  % use the conversion matrix to converse a 8*8 block of the image
 quanti_img = blkproc(dct_img,[8 8],'round(x./P1)',Norm_Mat);  % divide the block by the quantization matrix and reserve the aliquot part
 end
 
%----------------------------------------------------------
% DPCM on DC components
%----------------------------------------------------------
function DPCM=DPCMencode(quanti_img)
blocknumber=size(quanti_img,1)/8*size(quanti_img,2)/8;
a=1;
for i=1:8:size(quanti_img,1)
    for j=1:8:size(quanti_img,2)
         DctResult(1:8,a:a+7)=quanti_img(i:i+7,j:j+7);%transform the quanti_img matrix to a 8*8n matix, every 8*8 matix is a block.
         a=a+8;
    end 
end
count=2;%copy the first dc value of the 1st 8x8 block into dpcm1(1,1)
DPCM=zeros(1,blocknumber);
DPCM(1,1)=DctResult(1,1);
 for i=1:8:size(DctResult,2)-7
     TempDc=DctResult(1,i);%  copy first element of the current  8x8 block 
     if i==size(DctResult,2)-7
         break;
     end
     NextDc=DctResult(1,i+8); % copy the first element of  the next 8x8 block 
     DPCM(1,count)=NextDc-TempDc;% subtract the next element from the current element 
     count=count+1;

 end
end


%----------------------------------------------------------
% Huffman encoding for DPCM results
%----------------------------------------------------------
function [huff_DC,dict_DC]= HUFFMANencodeDC(DPCM)
p=tabulate(DPCM(:)); % use matlab's function 'tabulate()' to get the elements and their probability.
element=(p(:,1))'; %get the element in the DPCM sequence.
prob=(p(:,3)/100)'; % get the probability corresponding to its element. 
if sum(prob)~=1     % make sure the sum of the probability is 1.
    prob = prob /sum(prob);
end

[dict_DC, avglen1] = huffmandict(element,prob); % use the element and probability to generate the huffman dictionary
huff_DC = huffmanenco(DPCM, dict_DC);%use the original sequence and the dictionary to encode the sequence
end


%----------------------------------------------------------
% ZIGZAG scan for AC coefficients
%----------------------------------------------------------
function zig_seq=ZIGZAGscan(quanti_img)
blocknumber=size(quanti_img,1)/8*size(quanti_img,2)/8;
%define a sequence z which contains the Zigzag order to read the AC coefficients
z=[9 2 3 10 17 25 18 11 4 5 12 19 26 33 41 34 27 20 13 6 7 14 21 28 35 42 49 57 50 43 36 29 22 15 8 16 23 30 37 44 51 58 59 52 45 38 31 24 32 39 46 53 60 61 54 47 40 48 55 62 63 56 64]; 
v=1;
tempzig_seq=zeros(blocknumber,63);
for i=1:8:size(quanti_img,1)
    for j=1:8:size(quanti_img,2)
         A(i:i+7,j:j+7)=quanti_img(i:i+7,j:j+7);
         B(1:8,1:8)=A(i:i+7,j:j+7);% in the circulation,every time B is a 8*8 block of the quanti_img.
         tempzig_seq(v,1:63)=B(z); %read B's AC coefficients in order of the value of element in sequence z. 
         % in the tempzig_seq who is a n*63 matix, every row save the zigzag ordered 63 AC coefficients of a 8*8 block. 
         v=v+1;
    end
end

for i=1:size(tempzig_seq,1)
    zig_seq(1,(i-1)*63+1:i*63)=tempzig_seq(i,1:63); %transform the tempzig_seq to into the one-dimensional vector format. 
end
end

%----------------------------------------------------------
% Run-Length encoding for AC coefficients
%----------------------------------------------------------
function ac_rle=RunlengthEncode(zig_seq)
count=1; %represent the current position in run-length-encode sequence 
runlength=0; % represent the current number of zeros after the last non-zero component and before the next non-zero component
    for i=1:size(zig_seq,2)
            temp=zig_seq(1,i);
            if  temp==0
                if i==size(zig_seq,2)&&runlength>0 % if it is the last element of zig_seq
                    ac_rle(1,count)=runlength;
                    ac_rle(1,count+1)=0; 
                else
                    runlength=runlength+1;
                end                
            else 
                 ac_rle(1,count)=runlength; %save the number of zeros to the ac_rle sequence
                 ac_rle(1,count+1)=temp; %save the non-zero component to the ac_rle sequence
                 count = count+2;
                 runlength=0;
             end
     end
end

%----------------------------------------------------------
 %Huffman encoding for Run-length-encoding results
%----------------------------------------------------------
function [huff_AC,dict_AC]= HUFFMANencodeAC(ac_rle)
p=tabulate(ac_rle(:)); % use matlab's function 'tabulate()' to get the elements and their probability.
element=(p(:,1))'; %get the element in the ac_rle sequence.
prob=(p(:,3)/100)'; % get the probability corresponding to its element. 
if sum(prob)~=1     % make sure the sum of the probability is 1.
    prob = prob /sum(prob);
end

[dict_AC, avglen1] = huffmandict(element,prob); % use the element and probability to generate the huffman dictionary
huff_AC = huffmanenco(ac_rle, dict_AC);%use the original sequence and the dictionary to encode the sequence
end



