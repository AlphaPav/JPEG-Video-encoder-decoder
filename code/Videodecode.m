function Videodecode()

    readattribute=importdata('VideoAttribute.txt'); %read the video attribute
    frame_number=readattribute.frame_number; %define some attributes for recovering video afterwards
    frame_rate=readattribute.frame_rate;
    long=readattribute.long;
    height=readattribute.height;
    
    myObj = VideoWriter('recover.avi'); %create a new 'recover.avi' video
    myObj.FrameRate = frame_rate;
    
    cd(strcat('NewFigure')); %enter the folder 'NewFigure'
    open(myObj); %open the video 
    for i=1:frame_number
        myfilename = strcat('pic_',num2str(i),'.txt');
        ReadResult=importdata(myfilename); % load the encoding data from the document
        rec_frame = ImageDecompress(ReadResult,long,height); % recover the frame by recalling a function to process the data
        writeVideo(myObj,rec_frame); %write the recoverd frame into the new video
    end
    cd .. %return to the last folder
    close(myObj); %close the video
        
end



function rec_frame =ImageDecompress(ReadResult,long,height)

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
%load data from 'ReadResult'
    huff_DCY=ReadResult.huff_DCY;
    dict_DCY=ReadResult.dict_DCY;
    huff_ACY=ReadResult.huff_ACY;
    dict_ACY=ReadResult.dict_ACY;
    huff_DCU=ReadResult.huff_DCU;
    dict_DCU=ReadResult.dict_DCU;
    huff_ACU=ReadResult.huff_ACU;
    dict_ACU=ReadResult.dict_ACU;
    huff_DCV=ReadResult.huff_DCV;
    dict_DCV=ReadResult.dict_DCV;
    huff_ACV=ReadResult.huff_ACV;
    dict_ACV=ReadResult.dict_ACV;

% Decompression 
rec_Y=Decompression(huff_DCY,dict_DCY,huff_ACY,dict_ACY,long,height,Norm_Mat1);
rec_U=Decompression(huff_DCU,dict_DCU,huff_ACU,dict_ACU,long,height,Norm_Mat2);
rec_V=Decompression(huff_DCV,dict_DCV,huff_ACV,dict_ACV,long,height,Norm_Mat2);

% YUV to RGB
rec_img=zeros(long,height,3);
rec_img(:,:,1)=double(rec_Y)-0.001*double(rec_U)+1.402*double(rec_V);
rec_img(:,:,2)=double(rec_Y)-0.344*double(rec_U)-0.714*double(rec_V);
rec_img(:,:,3)=double(rec_Y)+1.772*double(rec_U)+0.001*double(rec_V);

% return the recovered image
rec_frame=uint8(rec_img);

end

%----------------------------------------------------------
% Image Decompression Overview
%----------------------------------------------------------
function rec_img=Decompression(huff_DC,dict_DC,huff_AC,dict_AC,long,height,Norm_Mat)
blocknumber=long/8*height/8; % calculate the total number of 8*8 size blocks in a long*height size image
re_RLE= huffmandeco(huff_AC, dict_AC);  %call the huffman decoding function 
re_zig_seq=RunlengthDecode(re_RLE,blocknumber); 
re_quanti_img=reverseZIGZAGscan(re_zig_seq,long,height);
re_DPCM= huffmandeco(huff_DC, dict_DC);
re_quanti_img=DPCMdecode(re_DPCM,re_quanti_img);
rec_img=IDCT(re_quanti_img,Norm_Mat);
end

%----------------------------------------------------------
% Run-Length decoding for AC coefficients
%----------------------------------------------------------
function re_zig_seq=RunlengthDecode(re_RLE,blocknumber)

count=1;
re_zig_seq=zeros(1,blocknumber*63);
for i=1:size(re_RLE,2)/2
    fir=re_RLE(1,2*i-1); %represents the number of zeros
    sec=re_RLE(1,2*i); %represents the next non-zero element
    if fir==0
        re_zig_seq(1,count)=sec;
        count=count+1;
    else
        for j=0:fir-1
            re_zig_seq(1,count+j)=0; % there are 'fir' zeros
        end
        re_zig_seq(1,count+fir)=sec;
        count=count+fir+1;
    end
end
end

%----------------------------------------------------------
%reverse Zigzag scaning for AC coefficients
%----------------------------------------------------------
function re_quanti_img=reverseZIGZAGscan(re_zig_seq,long,height)
blocknumber=long/8*height/8;
re_tempzig_seq=zeros(blocknumber,63);
for i=1:blocknumber
    for j=1:63
         re_tempzig_seq(i,j)=re_zig_seq(1,63*(i-1)+j); %transform the one-dimensional vector re_zig_seq  into the two-dimensional vector re_tempzig_seq.
    end
end
%define a sequence z which contains the Zigzag order to recover the AC coefficients
z=[9 2 3 10 17 25 18 11 4 5 12 19 26 33 41 34 27 20 13 6 7 14 21 28 35 42 49 57 50 43 36 29 22 15 8 16 23 30 37 44 51 58 59 52 45 38 31 24 32 39 46 53 60 61 54 47 40 48 55 62 63 56 64]; 
re_quanti_img=zeros(long,height);
count=1; % count will be from 1 to blocknumber
B=zeros(8,8);
for i=1:8:long
    for j=1:8:height
        for v=1:63
         B(z(v))=re_tempzig_seq(count,v);
        end
        re_quanti_img(i:i+7,j:j+7)=B(1:8,1:8);
        count=count+1;
    end
end
end


%----------------------------------------------------------
% DPCM decoding for DC coefficients
%----------------------------------------------------------
function re_quanti_img=DPCMdecode(re_DPCM,re_quanti_img)
blocknumber=size(re_quanti_img,1)/8*size(re_quanti_img,2)/8;
re_DC=zeros(1,blocknumber);
re_DC(1,1)=re_DPCM(1,1);
for i=2:blocknumber
    re_DC(1,i)=re_DPCM(1,i)+re_DC(1,i-1); % recover the DC coefficients
end

count=1;
for i=1:8:size(re_quanti_img,1)
    for j=1:8:size(re_quanti_img,2)
        re_quanti_img(i,j)=re_DC(1,count); %assign the DC coefficients to every corresponding  8*8 blocks in re_quanti_img
        count=count+1;
    end
end
end

%----------------------------------------------------------
% IDCT
%----------------------------------------------------------
 function rec_img=IDCT(re_quanti_img,Norm_Mat)
 T = dctmtx(8);   % create a 8*8 two-dimension IDCT conversion matrix T.
 redct_img = blkproc(re_quanti_img, [8 8], 'x .* P1', Norm_Mat);  % multiply the block by the quantization matrix.
 rec_img= blkproc(redct_img, [8 8], 'P1 * x * P2', T', T);% use the conversion matrix to converse the 8*8 block.
 end
 

