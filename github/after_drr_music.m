% restoredefaultpath
% 这个是不分频率的
clc;clear;close all;
% load hxkuaipai_nono_new5_noise0.mat
% load y_new5_noise0.mat
% addpath 'D:\DOA'
SNR = -10
hunxiangjishu = 0
for hunxiangwenjianming = 40:10:60
hunxiangjishu = hunxiangjishu+1;
yongjuzhenbiaoda = 0;
for xuanzejiaodu = 2:12
    xuanzejiaodu
    yongjuzhenbiaoda = yongjuzhenbiaoda + 1;
        filename = ['G:\SPP_contain_rever\' 'y_' num2str(SNR) 'db_0' num2str(hunxiangwenjianming) '_noise92_100_05s_ang' num2str(xuanzejiaodu) '.mat'];
        load(filename);
        filename2=['G:\SPP_contain_rever\result_small_use14\' 'result_IRM_contain_rever_' num2str(SNR) 'db_0' num2str(hunxiangwenjianming) '_ang' num2str(xuanzejiaodu) '.h5'];
        

%         data_irm = h5read(filename2,'/result_IRM');  %
% %         data_spp = permute(data_spp,[3,2,1]);
%         data_drr = h5read(filename2,'/result_DRR');  %
% %         data_drr = permute(data_drr,[3,2,1]);
        data_product = h5read(filename2,'/result_SPP_tongdaobufenkai_total');  %%%取错名字了
        data_product = permute(data_product,[3,2,1]);

%         data_DRR = h5read(filename2,'/result_DRR');  %%%取名字取错了  result_DRR实际是product
%         data_IRM = h5read(filename2,'/result_IRM');  %%%取名字取错了  result_DRR实际是product
%         data_product=data_IRM.*data_DRR;
%         data_product = permute(data_product,[3,2,1]);


%         total_mask = data_irm.*data_drr;
%         

        
        %%%%%%复原STFT 从129*8到256*8
       total_mask_full=[];
       mask_256=[];
       %这里是SPP的那种回
        for i = 1:8
            mask_part = data_product(((i-1)*129+1):(i*129),:,:);
            
            mask_fli = flipud(mask_part);
            mask_top = mask_fli(2:128,:,:); 
%             y_top_new = conj(y_top);
            total_mask = [mask_top ;mask_part];
            total_mask_full =[total_mask_full;total_mask];
        end
        
        %%%%%%%%%%%%把y全部换成乘以了mask之后的

        new_total_withnoise_signal=[];
        for yshuliang = 1 : 250
            for tongdao = 1 : 8
                y_each = total_withnoise_signal(:,tongdao,yshuliang);
                y_stft = stft(y_each,16000,Window=hann(256,"periodic"),OverlapLength=128,FFTLength=256);
%                 a=y_stft(128,:)
%                 b=y_stft(256,:)
                mask_dangqian256 = total_mask_full(((tongdao-1)*256+1):(tongdao*256),:,yshuliang);
                new_y_stft = y_stft .* mask_dangqian256;
                
                iy_new = istft(new_y_stft,16000,'Window',hann(256,"periodic"),'OverlapLength',128,'FFTLength',256);
                %%%%%    看一眼前后图
                new_total_withnoise_signal(:,tongdao,yshuliang)=iy_new; 
            end
        end



%%%%%%%%%上面变成了optimal statistic estimation        
m = 0;     
error_total=[];
method = 'SRP-PHAT';
%     method = 'MUSIC';
    % method = 'SNR-MVDR';
for yshuliang = 1:250
    y = new_total_withnoise_signal(:,:,yshuliang);
    % load ix_total_new-5.mat
%     addpath(genpath('./../'));
    % addpath('./wav files');
    %% 音频文件和传声器位置坐标
    % fileName = 'example.wav';  
    % micPos = ... 
    % ...%  mic1	 mic2   mic3   mic4   mic5   mic6   mic7  mic8
    %     [ 0.037 -0.034 -0.056 -0.056 -0.037  0.034  0.056 0.056;  % x
    %       0.056  0.056  0.037 -0.034 -0.056 -0.056 -0.037 0.034;  % y
    %     -0.038   0.038 -0.038  0.038 -0.038  0.038 -0.038 0.038]; % z
    
    micPos =  [0.5 2.72+1.5 1.8               %2米   0到180°
    0.5 2.8+1.5 1.8  
    0.5 2.88+1.5 1.8 
    0.5 2.96+1.5 1.8 
    0.5 3.04+1.5 1.8 
    0.5 3.12+1.5 1.8 
    0.5 3.20+1.5 1.8 
    0.5 3.28+1.5 1.8]'; 
    azBound = [-90 90]; % 方位角搜索范围
    elBound = 0;   % 俯仰角搜索范围。若只有水平面：则elBound=0;
    gridRes = 0.1;          % 方位角/俯仰角的分辨率
    alphaRes = 0.1;          % 分辨率
     

    wlen = 256;
    window = hann(wlen);
    noverlap = 0.5*wlen;
    nfft = 256;
    nsrc = 1;               % 声源个数
    c = 340;                % 声速
    freqRange = [ ];         % 计算的频率范围 []为所有频率
    pooling = 'sum';        % 如何聚合各帧的结果：所有帧取最大或求和{'max' 'sum'}
    
    %% 读取音频文件(fix)
    % [x,fs] = audioread(fileName);
    % x_full = hxkuaipai_nono;
    
    x_full = y;
    % x_full = ix_total;
    qiepian = 4480;
    for qiepianshu = 1 : 1
        
    x = x_full((qiepianshu-1)*qiepian+1:qiepianshu*qiepian,:);
    fs=16e3;
        [nSample,nChannel] = size(x);
        if nChannel>nSample, error('ERROR:输入信号为nSample x nChannel'); end
        [~,nMic,~] = size(micPos);
        if nChannel~=nMic, error('ERROR:麦克风数应与信号通道数相等'); end
        %% 保存参数(fix)
        Param = pre_paramInit(c,window, noverlap, nfft,pooling,azBound,elBound,gridRes,alphaRes,fs,freqRange,micPos);
        %% 定位(fix)
        if strfind(method,'SRP')
            specGlobal = doa_srp(x,method, Param);
        elseif strfind(method,'SNR')
            specGlobal = doa_mvdr(x,method,Param);
        elseif strfind(method,'MUSIC')
            specGlobal = doa_music1(x,Param,nsrc);
        else 
        end
        
        %% 计算角度
        minAngle                   = 10;         % 搜索时两峰之间最小夹角
        specDisplay                = 0;          % 是否展示角度谱{1,0}
        % pfEstAngles = post_sslResult(specGlobal, nsrc, Param.azimuth, Param.elevation, minAngle);
        % 绘制角谱
        % [pfEstAngles,figHandle] = post_findPeaks(specGlobal, Param.azimuth, Param.elevation, Param.azimuthGrid, Param.elevationGrid, nsrc, minAngle, specDisplay);
        [pfEstAngles,figHandle] = post_findPeaks(specGlobal, Param.azimuth, Param.elevation, Param.azimuthGrid, Param.elevationGrid, nsrc, minAngle, specDisplay);
        
        azEst = pfEstAngles(:,1)';
        elEst = pfEstAngles(:,2)';
%         for i = 1:nsrc
%             fprintf('切片数为：%d \n 第 %d 个声源方位为: \n Azimuth (Theta): %.0f \t Elevation (Phi): %.0f \n\n',qiepianshu,i,azEst(i),elEst(i));
%         end
        close all;
                if xuanzejiaodu == 2
                    trueang=-75;
                elseif xuanzejiaodu ==3
                    trueang=-60;
                elseif xuanzejiaodu ==4
                    trueang=-45;
                elseif xuanzejiaodu ==5
                    trueang=-30;
                elseif xuanzejiaodu ==6
                    trueang=-15;
                elseif xuanzejiaodu ==7
                    trueang=0;
                elseif xuanzejiaodu ==8
                    trueang=15;
                elseif xuanzejiaodu ==9
                    trueang=30;
                elseif xuanzejiaodu ==10
                    trueang=45;
                elseif xuanzejiaodu ==11
                    trueang=60;
                elseif xuanzejiaodu ==12
                    trueang=75;
                end
        result(yshuliang) = azEst;
        error_total(yshuliang) = abs(azEst-trueang);
        if error_total(yshuliang)<=5;
            m = m+1;
        end
    end
end
error_final = mean(error_total)
aaa_juzhenbiaoda(1,yongjuzhenbiaoda)=error_final;
zhunquelv = m/250
aaa_juzhenbiaoda(2,yongjuzhenbiaoda)=zhunquelv;
end
a_sanwei_juzhenbiaoda(:,:,hunxiangjishu)=aaa_juzhenbiaoda;
end
%%%%%%注意的是：STFT的图，最后一行是不对称的，而SPP是第一行不对称的