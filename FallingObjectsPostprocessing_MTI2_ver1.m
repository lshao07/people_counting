clear;clc;clf
 close all%
% 实现动态显示下落抛物 以及 下落路线
% 在22组数据中，结果不好
eps=1e-12;
%%% data acquisition
% homeDir='E:\ZJL_Sync\Cascade_Radar_外场测试\Dropping Objects\FallingDetecting\';
%homeDir='C:\Users\Siyis\Nextcloud\Cascade_Radar_外场测试\Dropping Objects\FallingDetecting\';
%dataDir=[homeDir,'RadarData\'];
dataDir=['C:\Users\ZJ\Documents\matlab_workspace\radar_falling_detect\'];
%filename=[dataDir,'BBTX1BaiDuCanFallingSSY13.bin'];
filename=['BBTX1BaiDuCanFallingSSY16.bin'];
% addpath('E:\TarfProjectFunc')
rawdata = (readDCA1000(filename));  %gpuArray

%% Parameter Setup
LightSpeed=3e8;
lamda = 3e8/77e9;
d = lamda/2;
% API:ProfileConfig:
% [],Bandwidth, idle_time(100.00us), Adc_Start_Time(6.00us), Ramp_End_Time(150.00us), [], [], duty cycle?, [], Samples_per_Chirp, Fs(ksps), [], [], RX_GAIN(dB),[]
% 0,1481987413,10000,600,15000,0,0,207,0,1024,10000,0,0,30,0,

[Samples_per_Chirp, Chirp_Duration, ADC_Fs, Num_chirp_loops, FramPeriod, Num_Frames, slope, BandWidth, R_Maximum, R_resulo, V_Maximum, V_resulo] ...
    =  RadarParamsExtract(dataDir,filename);

%%% Profile Configuration
% ProfileConfig,0,1435384036,300,500,2300,0,0,145,0,256,15000,2,1,30,0,
Num_Lanes=4;

% Samples_per_Chirp=256;
% Num_chirp_loops=128;
% Num_Frames=256;
Num_SamplesPerChannel = Samples_per_Chirp * Num_chirp_loops * Num_Frames;

% ADC_Fs=8000e3;
% ka=40.024e12;
% Bandwidth=Samples_per_Chirp/ADC_Fs*ka;

% Idle_Time=20e-6;
% Ramp_End_Time=40e-6;
% Chirp_Duration=Ramp_End_Time + Idle_Time; % Chirp_Duration_us = Ramp_End_Time_us + Idle_Time_us

% T_chirp_ADC=Samples_per_Chirp/ADC_Fs;

nfft_d= Samples_per_Chirp;   %2^nextpow2  回波信号采样点数 解距离最少点数
nfft_v= Num_chirp_loops;  %2^nextpow2  回波信号脉冲数 解速度最少点数
nfft_f = Num_Frames;  %2^nextpow2  回波信号脉冲数 解速度频率最少点数
nfft_angle = 180;
% t=(1:Num_SamplesPerChannel)*1/ADC_Fs;
StartFreq=77e9;
% Fcarrier=StartFreq+BandWidth/2;
% Range Doppler Parameters
ChirpBeat=ADC_Fs/nfft_d*(0:nfft_d-1); % Up chirp beat frequency
Range = LightSpeed/2*ChirpBeat/slope; % range vector
RadicalSpeed=V_Maximum*(-1:2/nfft_v:1-2/nfft_v); % radical speed vector

%% 2D-FFT
loops = Num_Frames;
selected_R_F=cell(1,loops);

win = repmat(hamming(nfft_v).',nfft_d,1);  %hamming(K)是一个列向量
amp = (0.5*(cos((-0.5*nfft_v+1:0.5*nfft_v)/nfft_v*2*pi)+1)).^2+1e-4;
win_amp = repmat(amp,nfft_d,1);  %
a = 1;  b = 0.5*[1, -1];
stat=zeros(nfft_v,nfft_d); 
data_nci=sum(rawdata,1)/size(rawdata,1);
ValiFramN = 1;
angle_list = [];
range_list = [];
for Frame_ID = 80:loops
   
%     data_mat=reshape(data_nci, Samples_per_Chirp, Num_chirp_loops, Num_Frames); %%% Range FFT
%     DataRangeFft = fft(data_mat(:,:,Frame_ID),nfft_d,1)/(nfft_d/2); %X[k]为个点的频率分量（除X[0]外），幅值=模值(X[k])/(N/2)
%     %%% MTI
%     mti_out = filter(b, a, DataRangeFft,[],2);                                        % 三脉冲对消运算（同一距离点在3个相邻脉冲重复间隔内的数据做相减运算）
%     % 形成由10个脉冲组成的二维数据矩阵（慢时间维10个脉冲，快时间维1000个距离点）
%     %%% Doppler FFT
%     DataVelocityFftNoMTI = fftshift(fft(DataRangeFft,nfft_v,2),2)/(nfft_v/2);%
%     DataVelocityFft = fftshift(fft(mti_out.*win,nfft_v,2),2).*win_amp/(nfft_v/2);%
%     
    for ChannlNum = 1:4    % receive channel choice
        data_mat_angle=reshape(rawdata(ChannlNum,:), Samples_per_Chirp, Num_chirp_loops, Num_Frames); %%% Range FFT
        DataRangeFft_angle(:,:,ChannlNum) = fft(data_mat_angle(:,:,Frame_ID),nfft_d,1)/(nfft_d/2); %X[k]为个点的频率分量（除X[0]外），幅值=模值(X[k])/(N/2)
        %%% MTI
        mti_out_angle = filter(b, a, DataRangeFft_angle(:,:,ChannlNum),[],2);                                        % 三脉冲对消运算（同一距离点在3个相邻脉冲重复间隔内的数据做相减运算）
        % 形成由10个脉冲组成的二维数据矩阵（慢时间维10个脉冲，快时间维1000个距离点）
        %%% Doppler FFT
        DataVelocityFftNoMTI_angle(:,:,ChannlNum) = fftshift(fft(DataRangeFft_angle(:,:,ChannlNum),nfft_v,2),2)/(nfft_v/2);%
        DataVelocityFft_angle(:,:,ChannlNum) = fftshift(fft(mti_out_angle.*win,nfft_v,2),2).*win_amp/(nfft_v/2);%
    end
    DataVelocityFftNoMTI = sum(DataVelocityFftNoMTI_angle,3)/size(DataVelocityFftNoMTI_angle,3);
    DataVelocityFft = sum(DataVelocityFft_angle,3)/size(DataVelocityFft_angle,3);
    for n=1:nfft_d   %range
        for m=1:nfft_v  %chirp
            temp=DataVelocityFft_angle(n,m,:);    
            temp_fft=fftshift(fft(temp,nfft_angle));    %对2D FFT结果进行nfft_angle点FFT
            angle_profile(n,m,:)=temp_fft;  
        end
    end
%     RD_dB = 20*log10(abs(DataVelocityFftNoMTI)) ;
%     RD_dB(:,65) = RD_dB(:,2);
%     
%     RDMTI_dB = 20*log10(abs(DataVelocityFft)) ;
%     RDAll_dB = (RD_dB+RDMTI_dB);    
    % ==========================================================================
    
    %% CFAR
    amp_mat = abs(DataVelocityFft).'+abs(DataVelocityFftNoMTI).' ;%
%    amp_mat = abs(DataVelocityFft).' ;
    [M,N]=size(amp_mat);
    %%%%%%%%%二维筛选%%%%%%%%%
    N_ref_2D = M/4;  %参考区域大小
    clear DataVelocityFft DataVelocityFftNoMTI 
    %计算区域1的和
    sum1 = sum(amp_mat(1:N_ref_2D,1:N_ref_2D)); %按行相加
    sum_ref_2D(1,1) = sum(sum1(1,:));            %按列相加
    %计算区域2的和
    sum2 = sum(amp_mat(1:N_ref_2D,(N-N_ref_2D+1):N));
    sum_ref_2D(1,2) = sum(sum2(1,:));
    %计算区域3的和
    sum3 = sum(amp_mat((M-N_ref_2D+1):M,1:N_ref_2D));
    sum_ref_2D(1,3) = sum(sum3(1,:));
    %计算区域4的和
    sum4 = sum(amp_mat((M-N_ref_2D+1):M,(N-N_ref_2D+1):N));
    sum_ref_2D(1,4) = sum(sum4(1,:));
    
    %%%产生门限Threshold%%%
    SNR_Threshold = 10;   % 检测门限 dB
    Threshold = ((sum(sum_ref_2D(:,1:4)) - max(sum_ref_2D))/(N_ref_2D^2*3))*10^(SNR_Threshold/10);
    
    %%%寻找大于门限的点并记录其位置%%%
    location = zeros(M,N);  %填1表示该点大于门限
    for i=1:M
        for j=1:N
            if(amp_mat(i,j) >= (1*Threshold))
                location(i,j) = 1;
            else
                location(i,j) = 0;
            end
        end
    end
    %%%筛选结果显示%%%
    selected = zeros(M,N);
    for i=1:1:M
        for j=1:1:N
            if(location(i,j) == 1)
                selected(i,j) = amp_mat(i,j);
            else
                selected(i,j) = 0;
            end
        end
    end
    %%%%%%%%%%一维筛选%%%%%%%%%%
    %%%%距离维筛选%%%%
    amp_mat = selected;
    guardLen_R = 8; %保护单元
    referLen_R = 8;  %参考单元
    threshold_R = zeros(M,N); %记录距离维每个观察点的门限值
    
    selected_R = zeros(M,N);
    for i=1:M
        for j=1:N
            if(j <= (guardLen_R + referLen_R) )                     %观察点过于偏左，左侧点数不足
                threshold_R(i,j) = (sum(amp_mat(i,(j+guardLen_R+1):(j+guardLen_R+referLen_R))))/referLen_R; %观察点右侧参考点相加求平均
            elseif( j >= (N - guardLen_R - referLen_R + 1))         %观察点过于偏右，右侧点数不足
                threshold_R(i,j) = (sum(amp_mat(i,(j-guardLen_R-referLen_R):(j-guardLen_R-1))))/referLen_R; %观察点左侧参考点求平均
            else                                                    %观察点居中，左右点数均足够
                sum_R_left = sum(amp_mat(i,(j-guardLen_R-referLen_R):(j-guardLen_R-1)));
                sum_R_right = sum(amp_mat(i,(j+guardLen_R+1):(j+guardLen_R+referLen_R)));
                threshold_R(i,j) = max(sum_R_left,sum_R_right)/(referLen_R); %观察点两侧取大的一侧求均值
            end
            
            %利用距离维门限进行一维筛选
            if(amp_mat(i,j) >= (1*threshold_R(i,j)))
                selected_R(i,j) = amp_mat(i,j);
            else
                selected_R(i,j) = 0;
            end
        end
    end
    %%%%%%%%%%%距离维筛选结束%%%%%%%%%%%%%
    
    %%%%%频率维筛选%%%%%
    guardLen_F = 8;
    referLen_F = 8;
    threshold_F = zeros(M,N); %记录距离维每个观察点的门限值
    selected_R_F{Frame_ID} = zeros(M,N); 
    objN= 0;        posi =[];
    for j=1:N
        for i=1:M
            if selected_R(i,j)~=0
                if(i <= (guardLen_F + referLen_F) )                 %观察点过于偏左，左侧点数不足
                    threshold_F(i,j) = (sum(amp_mat((i+guardLen_F+1):(i+guardLen_F+referLen_F),j)))/referLen_F; %观察点右侧参考点相加求平均
                elseif( i >= (M - guardLen_F - referLen_F + 1))     %观察点过于偏右，右侧点数不足
                    threshold_F(i,j) = (sum(amp_mat((i-guardLen_F-referLen_F):(i-guardLen_F-1),j)))/referLen_F; %观察点左侧参考点求平均
                else                                                %观察点居中，左右点数均足够
                    sum_F_left = sum(amp_mat((i-guardLen_F-referLen_F):(i-guardLen_F-1),j));
                    sum_F_right = sum(amp_mat((i+guardLen_F+1):(i+guardLen_F+referLen_F),j));
                    threshold_F(i,j) = max(sum_F_left,sum_F_right)/(referLen_F); %观察点两侧取大的一侧求均值
                end
                %利用频率维门限进行一维筛选
              
                if(amp_mat(i,j) >= (1*threshold_F(i,j)))
                    selected_R_F{Frame_ID}(i,j) = amp_mat(i,j);
                   if abs(i-65)>3 
                       objN=objN+1;  
                      posi(objN,:) = [RadicalSpeed(i)-0.4 Range(j)-0.4 0.8 0.8 ,j , i];
 %                          posi(objN,:) = [-4 0 0.01 0.01]; 
                   end
                else
                    selected_R_F{Frame_ID}(i,j) = 0;
                end
            end
        end
    end  

 stat=stat+selected_R_F{Frame_ID};

 %% figure
%imagesc(RadicalSpeed,Range,RDAll_dB );
imagesc(RadicalSpeed,Range,20*log10(selected_R_F{Frame_ID}).' );  
if objN>0
    if ValiFramN==1
        [mv,ml] = min(posi(:,1)-posi(:,2)) ;%选取速度最小，高度最高的作为初始点
        if posi(ml,1)<0
            %valiPosi(ValiFramN,:) = posi(ml,:);
            %rectangle('Position',valiPosi(ValiFramN,:),'EdgeColor','r','LineWidth',1);
            valiPosi(ValiFramN,:) = posi(ml,:);
            rectangle('Position',valiPosi(ValiFramN,1:4),'EdgeColor','r','LineWidth',1);
            [val_angle,idx_angle] = max(angle_profile(valiPosi(ValiFramN,5),valiPosi(ValiFramN,6),:));
            fw = (idx_angle-nfft_angle/2-1)/nfft_angle;         %空间频率
            theta = asin(fw*lamda/d);  %角度公式
            angle = theta*180/pi;
            angle_list = cat(1,angle_list,angle);
            range_list = cat(1,range_list,valiPosi(ValiFramN,2));
            fprintf('目标角度： %f°\n',angle);
            ValiFramN = ValiFramN+1;%
        end
    else
        if objN == 1
            if posi(ml,2)<=valiPosi(ValiFramN-1,2)&&abs(posi(ml,1)-valiPosi(ValiFramN-1,1))<0.5
                valiPosi(ValiFramN,:) = posi(ml,:);
                rectangle('Position',valiPosi(ValiFramN,1:4),'EdgeColor','r','LineWidth',1);
                [val_angle,idx_angle] = max(angle_profile(valiPosi(ValiFramN,5),valiPosi(ValiFramN,6),:));
                fw = (idx_angle-nfft_angle/2-1)/nfft_angle;         %空间频率
                theta = asin(fw*lamda/d);  %角度公式
                angle = theta*180/pi;
                angle_list = cat(1,angle_list,angle);
                range_list = cat(1,range_list,valiPosi(ValiFramN,2));
                fprintf('目标角度： %f°\n',angle);
                ValiFramN = ValiFramN+1;%     d
            end
        else
            %            for nn = 1:objN
            [mv,ml]= min(norm(posi(:,1:2)-valiPosi(ValiFramN-1,1:2))) ;
            if posi(ml,2)<=valiPosi(ValiFramN-1,2)
                valiPosi(ValiFramN,:) = posi(ml,:);
                rectangle('Position',valiPosi(ValiFramN,1:4),'EdgeColor','r','LineWidth',1);
                [val_angle,idx_angle] = max(angle_profile(valiPosi(ValiFramN,5),valiPosi(ValiFramN,6),:));
                fw = (idx_angle-nfft_angle/2-1)/nfft_angle;         %空间频率
                theta = asin(fw*lamda/d);  %角度公式
                angle = theta*180/pi;
                angle_list = cat(1,angle_list,angle);
                range_list = cat(1,range_list,valiPosi(ValiFramN,2));
                fprintf('目标角度： %f°\n',angle);
                ValiFramN = ValiFramN+1;%
            end       
        end
    end
end
% delete(h)
set(gca,'YDir','normal')  % 显示图像时，纵坐标本来就是翻转的
title('R-D Diagram')
xlabel('Speed (m/s)')
ylabel('Range (m)')
colorbar
% caxis([-250 130])
title (['Frame' num2str(Frame_ID)] )
drawnow
Trajectory= getframe;

end

 %% ========================================================================
figure
imagesc(RadicalSpeed,Range,20*log10(stat.'))
set(gca,'YDir','normal')
% title('Trajectory with MTI')
xlabel('Speed (m/s)')
ylabel('Range (m)')
colorbar
caxis([-250 130])
% hold on
% plot(valiPosi(:,1)+0.4,valiPosi(:,2)+0.4,'r*--')
figure 
c = linspace(1,5,length(angle_list));
sz = 15;
scatter(angle_list,range_list,sz,c)
axis([-90 90 0 30]);
xlabel('Angle (degree)')
ylabel('Range (m)')

figure
plot(angle_list)
xlabel('frame contains target')
ylabel('Angle (degree)')

idx_SSY=strfind(filename,'SSY');
idx_bin=strfind(filename,'.bin');
file_id=str2double(filename(idx_SSY+3:idx_bin-1));
figfilename=([dataDir,'Trajectory with MTI\',sprintf('data_%d.eps',file_id)]);
% saveas(gcf,figfilename,'psc2')
% 生成视频，并保存
writerObj =VideoWriter('Trajectory.avi'); % 生成一个avi动画
open(writerObj); % 打开avi动画
writeVideo(writerObj,Trajectory); % 将保存的动画写入到视频文件中
close(writerObj); % 关闭动画

