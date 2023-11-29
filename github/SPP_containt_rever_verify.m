% restoredefaultpath
clc;
clear all;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%信道%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for SNR = -20:5: -5;
for hunxiang = 30:20:50

hunxiangfile = ['h3_0' num2str(hunxiang) '_newnew.mat']
c = 340;                    % Sound velocity (m/s)
fs = 16e3;                 % Sample frequency (samples/s)
kuaipai = 250;
% SNR=randi([-5 20]);
load noise92_new.mat;
jiequshijian = 0.3;
r = [0.5 2.72+1.5 1.8               %2米   0到180°
    0.5 2.8+1.5 1.8  
    0.5 2.88+1.5 1.8 
    0.5 2.96+1.5 1.8 
    0.5 3.04+1.5 1.8 
    0.5 3.12+1.5 1.8 
    0.5 3.20+1.5 1.8 
    0.5 3.28+1.5 1.8];                % Receiver position [x y z] (m)
% r = [1 2.72 1.2];
 micPos = r'; 
%
s1 = [0.5              0.5                        1.8             % Source position [x y z] (m)
0.5+4*cos(75*pi/180)   4.5-4*sin(75*pi/180)       1.8
0.5+4*cos(60*pi/180)   4.5-4*sin(60*pi/180)       1.8
0.5+4*cos(45*pi/180)   4.5-4*sin(45*pi/180)       1.8
0.5+4*cos(30*pi/180)   4.5-4*sin(30*pi/180)       1.8
0.5+4*cos(15*pi/180)   4.5-4*sin(15*pi/180)       1.8

4.5                  4.5                        1.8

0.5+4*cos(15*pi/180)   4.5+4*sin(15*pi/180)       1.8
0.5+4*cos(30*pi/180)   4.5+4*sin(30*pi/180)       1.8
0.5+4*cos(45*pi/180)   4.5+4*sin(45*pi/180)       1.8
0.5+4*cos(60*pi/180)   4.5+4*sin(60*pi/180)       1.8
0.5+4*cos(75*pi/180)   4.5+4*sin(75*pi/180)       1.8
0.5                    8.5                        1.8];
%%%
% plot(s1(:,1),s1(:,2),'o'),hold on 
% xlim([0 9]);ylim([0 9])
% plot(r(:,1),r(:,2),'o'),hold on
%%%%

L = [6 11 3.6];                % Room dimensions [x y z] (m)
beta2 = 0;     
beta = 0.7;                 % Reverberation time (s)
n = 10000;                   % Number of samples
mtype = 'omnidirectional';
order = -1;
order2 = 0;
load h3_nono_newnew.mat
load(hunxiangfile);

%%%%%%%%%数据
fileFolder = fullfile('D:\DOA\timit\signal\test_wav_16k'); %搜索目录
dirOutput = dir(fullfile(fileFolder,'*.wav'));%获取目录下所有“wav”格式音频文件信息
fileNames = {dirOutput.name};%获取音频文件的名字，放入数组fileNames中
filePath = {dirOutput.folder};%获取音频文件目录，存放入数组filePath中
total_withnoise=[];
fd=[];
hxtotal=[];
hxtotal_hunxiang=[];
hxtotal_nono=[];
hxtotal_noise=[];
meizushuliang= 250;
for weizhi =  2 : 12
    count_i = 0 ;
    for i = 1 : meizushuliang
          weizhi
          count_i = count_i + 1
          %%
          load(hunxiangfile)
          %%
          h3_hunxiang = h3-h3_nono;
          aa = fullfile(filePath(1,i),fileNames(1,i));
          wavepath = cell2mat(aa(1,1));
          [audio,fd(count_i)]=audioread(wavepath);
            for zhenyuan = 1 : 8     %没混响与多径
              hx_nono = conv(audio, h3_nono(:,zhenyuan,weizhi));
              hxtotal_nono = [hxtotal_nono,hx_nono];  
            end
            for zhenyuan = 1 : 8                      %%%有混响与多径
              hx = conv(audio, h3(:,zhenyuan,weizhi));
              hxtotal = [hxtotal,hx];  
            end
            for zhenyuan = 1 : 8                      %%%只有混响无多径
              hx_hunxiang = conv(audio, h3_hunxiang(:,zhenyuan,weizhi));
              hxtotal_hunxiang = [hxtotal_hunxiang,hx_hunxiang];    
            end
         
         hxkuaipai = hxtotal(1:jiequshijian*fs,:);            %有直径有多径   没噪声
        
         hxkuaipai_nono = hxtotal_nono(1:jiequshijian*fs,:);  %%仅有直径    没噪声
         hxkuaipai_hunxiang = hxtotal_hunxiang(1:jiequshijian*fs,:);  %%%仅有多径   没噪声 
        %%%%%%%%%%%%%%%%%%%%%%%%%5
          %%%%得到DRR
           changdu = length(totalnoise1);
           
           noiseSig=[];
           for channel = 1 : 8
              suijishu = randi([fs*jiequshijian/2+1 changdu-fs*jiequshijian/2-1]);
              noise = totalnoise1(suijishu-fs*jiequshijian/2:suijishu+fs*jiequshijian/2-1,1);
              noiseSig(channel,:) = noise';
              input = hxkuaipai';   %%%有混响有直达
              [output_noisy(channel,:),noiseSig(channel,:)]=add_noisem2(input(channel,:),noiseSig(channel,:),SNR);
              y = output_noisy';       
           end

    %%%%%%%%%%%
     noise_kuaipai =noiseSig';
    noise_contain_rever=noise_kuaipai + hxkuaipai_hunxiang;
            
            channel8_drr =[];
            channel8_spp =[];
            channel8_lps =[];
        %%%%%%%%%%%%%%%      

       for channel_new = 1:8
            direct_signal=hxkuaipai_nono(:,channel_new) ;
            rever_signal = hxkuaipai_hunxiang(:,channel_new);
            noise_signal = noise_kuaipai(:,channel_new);
            dir_re_signal = hxkuaipai(:,channel_new);
            y_signal=y(:,channel_new);
            noise_contain_rever_signal = noise_contain_rever(:,channel_new);
                
            dir_stft    = stft(direct_signal,16000,Window=hann(256,"periodic"),OverlapLength=128,FFTLength=256);
            rever_stft = stft(rever_signal,16000,Window=hann(256,"periodic"),OverlapLength=128,FFTLength=256);       
            noise_stft = stft(noise_signal,16000,Window=hann(256,"periodic"),OverlapLength=128,FFTLength=256);
            dir_re_stft = stft(dir_re_signal,16000,Window=hann(256,"periodic"),OverlapLength=128,FFTLength=256);
            y_stft = stft(y_signal,16000,Window=hann(256,"periodic"),OverlapLength=128,FFTLength=256);

%           [noisePowMat_clean, snrPost1_clean, SPP_clean, SPP2_clean] = noisepowproposed(direct_signal, fs);   %%%SPP
%             [noisePowMat,snrPost1,SPP] = noisepowproposed_new(y_signal, 16000, noise_signal);
            [noisePowMat,snrPost1,SPP] = noisepowproposed_new(y_signal, 16000, noise_contain_rever_signal); 

            dir_stft2=abs(dir_stft);
            rever_stft2=abs(rever_stft);
            noise_stft2=abs(noise_stft);
            dir_re_stft2=abs(dir_re_stft);
            y_stft2 = abs(y_stft);
         %%%
            dir_psd=dir_stft2.^2;
            rever_psd=rever_stft2.^2;
            noise_psd=noise_stft2.^2;
            dir_re_psd=dir_re_stft2.^2;
            y_psd = y_stft2.^2;
            
            log_LPS = log10(y_psd);   
            ratio_DRR_nolog = dir_psd./rever_psd;   %%%%
            ratio_SPP=SPP;  
            
         %%%%%%  DRR变成mask的样子方便拟合
            c_para=5;
            rou_para=-1.2;
            top=10^(c_para*rou_para/10);
            drr_mask = top./(top+ratio_DRR_nolog.^rou_para); %%%这里已经是做了sigmoid 之后的drr了


         %%%half
            duichen_DRR = drr_mask(128:256,:);
                        
            duichen_SPP = ratio_SPP;
                        
            duichen_LPS = log_LPS(128:256,:);
            
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%             channel8_drr = [channel8_drr;duichen_DRR];
            channel8_spp = [channel8_spp;duichen_SPP];
            channel8_lps = [channel8_lps;duichen_LPS];
         %%%%%%%%%%%%%%%       DRR and SPP get  and  ratio_new get
          end
            total_withnoise_lps(:,:,count_i,weizhi)=channel8_lps;
%             total_withoutnoise_drr(:,:,count_i,weizhi)=channel8_drr;
            total_spp(:,:,count_i,weizhi)=channel8_spp;

            total_withnoise_signal(:,:,count_i)=y;

            hxtotal=[];
            hxtotal_nono=[];
            hxtotal_noise=[];
            hxtotal_hunxiang=[];
         
    end

%%%切片并分别储存
[cun_1,cun_2,cun_3,cun_4] = size(total_withnoise_lps);
t_separate = 36;
dimen5 = floor(cun_2/t_separate);
numnum_new = 1;
for hh4 = cun_4 : cun_4
    for hh3 = 1 : cun_3
        for hh5 = 1 :dimen5
        qiepian_LPS  = total_withnoise_lps(:, (hh5-1)*t_separate+1:hh5*t_separate, hh3, hh4);
        x_LPS(:,:,numnum_new)  = qiepian_LPS;
        qiepian_spp  = total_spp(:, (hh5-1)*t_separate+1:hh5*t_separate, hh3, hh4);
        target_spp(:,:,numnum_new)  = qiepian_spp;
        y_label(numnum_new)   =  hh4-1;
        %%%增加DRR的切片
%         qiepian_DRR = total_withoutnoise_drr(:,(hh5-1)*t_separate+1:hh5*t_separate, hh3, hh4);
%         target_DRR(:,:,numnum_new) = qiepian_DRR;
        %%%%%%%%%%%%%
        numnum_new = numnum_new + 1;
        end
    end
end
%%%对应的SPP
%%%%%%%%%%%存下来
numnum_cun = numnum_new-1;

file_name = ['verifydata_' num2str(SNR) 'db_0' num2str(hunxiang) '_noise92_100_05s_ang' num2str(weizhi) '.h5']

h5create(file_name,'/x_verify',[1032 t_separate numnum_cun]);
h5create(file_name,'/y_verify',[1 numnum_cun],'Datatype','int32');
h5create(file_name,'/target_SPP',[1032 t_separate numnum_cun]);
% h5create(file_name,'/target_DRR',[1032 t_separate numnum_cun]);

h5write(file_name,'/x_verify',x_LPS);
h5write(file_name,'/y_verify',y_label);
h5write(file_name,'/target_SPP',target_spp);
% h5write(file_name,'/target_DRR',target_DRR);

file_name2 = ['y_' num2str(SNR) 'db_0' num2str(hunxiang) '_noise92_100_05s_ang' num2str(weizhi) '.mat'] 

save(file_name2,'total_withnoise_signal')



total_withoutnoise_LPS=[];
total_withnoise_LPS=[];
total_withoutnoise_STFT=[];
total_withnoise_STFT=[];
total_withoutnoise_SPP=[];
total_withnoise_SPP=[];
total_withnoise_signal=[];
total_DRR=[];

total_withnoise_lps=[];
total_spp=[];
total_withoutnoise_drr=[];
target_DRR=[];
target_spp=[];
y_label = [];
x_LPS = [];
qiepian_spp=[];
qiepian_DRR=[];
qiepian_LPS=[];

end
end
end
%%%%%%注意的是：STFT的图，最后一行是不对称的，而SPP是第一行不对称的

%改beta、SNR、路径名称（train与train）、最后存的名字、h3是否load、分频。


