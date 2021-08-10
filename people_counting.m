clear all;
close all;

%%% data acquisition
dataDir=['D:\matlab_workspace\people_counting\iwr6843_08_03_test_data\'];

filename=[dataDir,'adc_data_0.bin'];

% rawdata = readDCA1000_6843(filename); %gpuArray
load rawdata.mat

%% Parameter Setup
LightSpeed=3e8;
lambda = 3e8/60.75e9;
d = lambda/2;
azi_bin = -70:0.75:70;%nu
ele_bin = -20:0.75:20;%mu
azi_init = -sind(70);
ele_init = -sind(20);
azi_step = -2*azi_init/length(azi_bin);
ele_step = -2*ele_init/length(ele_bin);
rfFreqScaleFactor = 2.7;
Num_Chirp_config = 3;
Num_TX = 3;
Num_RX = 4;
m_ind = [0, 0,-1,-1,-2,-2,-3,-3,-2,-2,-3,-3];
n_ind = [0,-1,-1, 0, 0,-1,-1, 0,-2,-3,-3,-2]; 
phaseRot = [1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1];
nu = sind(azi_bin);
for i = 1:length(nu)
    tmp_az = exp(1i*pi.*m_ind*nu(i));
    sv_az(i,:) = phaseRot.*tmp_az;
end
mu = sind(ele_bin);
for i = 1:length(mu)
    tmp_el = exp(1i*pi.*n_ind*mu(i));
    sv_el(i,:) = tmp_el;
end
sv_az = sv_az.';
sv_el = sv_el.';
%sv_ele = sv_ele .* [1;-1;1;-1];

[Samples_per_Chirp, Chirp_Duration, ADC_Fs, Num_chirp_loops, FramPeriod, Num_Frames, slope, BandWidth, R_Maximum, R_resulo, V_Maximum, V_resulo] ...
    =  RadarParamsExtract_6843(dataDir,filename,rfFreqScaleFactor,Num_Chirp_config,Num_TX);



Num_Chirp_config = 3;
Num_SamplesPerChannel = Samples_per_Chirp * Num_chirp_loops * Num_Frames;
nfft_d = Samples_per_Chirp;   %2^nextpow2  回波信号采样点数 解距离最少点数
nfft_v = Num_chirp_loops;
ChirpBeat=ADC_Fs/nfft_d*(0:nfft_d-1); % Up chirp beat frequency
Range = R_Maximum*(0:nfft_d-1)/nfft_d;%LightSpeed/2*ChirpBeat/slope; % range vector
RadicalSpeed=V_Maximum*(-1:2/nfft_v:1-2/nfft_v); 

Discard_R_left = 4;
Discard_R_right = 4;
Discard_AZ_left = 2;
Discard_AZ_right = 2;
Ref_size_R = 8;
Ref_size_AZ = 12;
Guard_size_R = 4;
Guard_size_AZ = 8;
range_th = 5;
az_th = 8/Ref_size_AZ;
sidelobe_th = 0.4;
Firstpass_select1_4_5_8 = true;
Secondpass_select8_7_11_12 = true;

if Firstpass_select1_4_5_8
    sv_az_final = cat(1,sv_az(1,:),sv_az(4,:),sv_az(5,:),sv_az(8,:));
else
    sv_az_final = cat(1,sv_az(2,:),sv_az(3,:),sv_az(6,:),sv_az(7,:));
end

if Secondpass_select8_7_11_12
    sv_el_final = cat(1,sv_el(8,:),sv_el(7,:),sv_el(11,:),sv_el(12,:));
else
    sv_el_final = cat(1,sv_el(5,:),sv_el(6,:),sv_el(9,:),sv_el(10,:));
end
az_select_chirp_cfg = [1,2];
el_select_chirp_cfg = [2,3];

detected_points_total = [];
%%
for Frame_ID = 147:147%Num_Frames
    
    detected_points_info_per_frame = [];
 for ChannlNum = 1:4    % receive channel choice
        IQ_data(:,:,:,:,ChannlNum)=reshape(rawdata(ChannlNum,:), Samples_per_Chirp, Num_Chirp_config, Num_chirp_loops,Num_Frames); 
        
        %%1D FFT
        DataRangeFft(:,:,:,ChannlNum) = fft(IQ_data(:,:,:,Frame_ID,ChannlNum),nfft_d,1)/(nfft_d/2); %X[k]为个点的频率分量（除X[0]外），幅值=模值(X[k])/(N/2)
        
        %%Static Clutter Removal
        Average_samples = Static_Removal(DataRangeFft(:,:,:,ChannlNum));
        DataRangeFft_Static_removed(:,:,:,ChannlNum) = DataRangeFft(:,:,:,ChannlNum) - Average_samples;
%         figure
%         plot(abs(Average_samples(:,1)))
%         figure
%         plot(abs(DataRangeFft_Static_removed(:,1,1,ChannlNum)))
 end
        %%Capon
        %Sptial cov matrix
        
%         for ChirpNum = 1:Num_chirp_loops
%                 for ChirpConfig = 1:Num_TX
%                     Xn(:,:,ChirpNum) = squeeze(DataRangeFft(:,ChirpConfig,ChirpNum,:)).';%%暂时使用TX2 后续TX2-TX1-TX0
%                     Xn_sq(:,:,ChirpConfig,ChirpNum) = Xn(:,:,ChirpNum)*Xn(:,:,ChirpNum)' ./ 96;
%                 end
%         end
%         R_xx(:,:,Frame_ID) = mean(mean(Xn_sq,3),4);
%         
%         DLF = 0.5*std(diag(R_xx(:,:,Frame_ID)));
%         R_xx(:,:,Frame_ID) = R_xx(:,:,Frame_ID) + DLF*(trace(R_xx(:,:,Frame_ID))/Num_RX).*eye(4);
%         R_AZ_heatmap(:) = diag(1./(sv'*(R_xx(:,:,Frame_ID)\sv)));
%         weight_matrix(:,:) = R_xx(:,:,Frame_ID)\sv .* R_AZ_heatmap;


            
        for SampleNum = 1:Samples_per_Chirp
            for ChirpNum = 1:Num_chirp_loops
                %TX1,TX0 MIMO阵列
                Xn(:,:,ChirpNum) = cat(1,squeeze(DataRangeFft_Static_removed(SampleNum,az_select_chirp_cfg(1),ChirpNum,1)),squeeze(DataRangeFft_Static_removed(SampleNum,az_select_chirp_cfg(1),ChirpNum,4)),squeeze(DataRangeFft_Static_removed(SampleNum,az_select_chirp_cfg(2),ChirpNum,1)),squeeze(DataRangeFft_Static_removed(SampleNum,az_select_chirp_cfg(2),ChirpNum,4)));
                Xn_sq(:,:,ChirpNum) = Xn(:,:,ChirpNum)*Xn(:,:,ChirpNum)';
            end
            R_xx(:,:,SampleNum,Frame_ID) = mean(Xn_sq,3);
             %diagonal loading is applied to the R matrix to ensure stability
            %DLF = 0.2*std(diag(R_xx(:,:,SampleNum,Frame_ID)));
            DLF = 0.001;
            R_xx(:,:,SampleNum,Frame_ID) = R_xx(:,:,SampleNum,Frame_ID) + DLF*(trace(R_xx(:,:,SampleNum,Frame_ID))/Num_RX).*eye(4);
            %R_xx(:,:,Frame_ID) = Spatial_cov(IQ_data,FrameID);
            
            %Range_Azimuth_Heatmap
            R_xx_inv = inv(R_xx(:,:,SampleNum,Frame_ID));
            R_AZ_heatmap(SampleNum,:) = diag(1./(sv_az_final'*R_xx_inv*sv_az_final));
            %weight_matrix(SampleNum,:,:) = (R_xx_inv * sv_az_final) .* R_AZ_heatmap(SampleNum,:);
        end
        
       R_AZ_pwr_heatmap = abs(R_AZ_heatmap);
       figure(1);
       mesh(abs(R_AZ_heatmap))
       %weight = squeeze(mean(weight_matrix,1));% not sure
       
       %%2 pass CFAR
       [M,N] = size(R_AZ_pwr_heatmap);
       % range
       %单侧guard 和 ref cell 长度
       Guard_length_R = Guard_size_R/2;
       Ref_length_R = Ref_size_R/2;
       Guard_length_AZ = Guard_size_AZ/2;
       Ref_length_AZ = Ref_size_AZ/2;
       
       detection_map = zeros(M,N);
       detect_az_list = [];
       detect_r_list = [];
       SNR_list = [];
       %CASO First pass

       for i=1+Discard_AZ_left:N-Discard_AZ_right
           for j=1+Discard_R_left:M-Discard_R_right
                if j <= 1 + Discard_R_left + Guard_length_R %太靠左，左侧无Ref cell 
                    R_mean_right = mean(R_AZ_pwr_heatmap(j + Guard_length_R + 1 :2: j + Guard_length_R + Ref_length_R,i),1);
                    R_mean_left = R_mean_right;
                    %sum_left = sum_right - R_AZ_heatmap(j + Guard_length + 1,i) + R_AZ_heatmap(j + Guard_length + Ref_length + 1,i);
                elseif j > 1 + Discard_R_left + Guard_length_R && j < 1 + Discard_R_left + Guard_length_R + Ref_length_R %太靠左，左侧Ref cell个数不足
                    R_mean_right = mean(R_AZ_pwr_heatmap(j + Guard_length_R + 1 : j + Guard_length_R + Ref_length_R,i),1);
                    R_mean_left = R_mean_right;%mean(R_AZ_pwr_heatmap(Discard_R_left + 1:j - Guard_length_R - 1,i),1);
                    %mean_left = (mean_right*Ref_length + sum(R_AZ_pwr_heatmap(Discard_R_left + 1:j - Guard_length - 1,i)))/(Ref_length + j - 1 - Discard_R_left - Guard_length);
                elseif j >= M - Discard_R_right - Guard_length_R % 太靠右，右侧无Ref cell
                    R_mean_left = mean(R_AZ_pwr_heatmap(j - Guard_length_R - Ref_length_R : j - Guard_length_R - 1,i),1);
                    R_mean_right = R_mean_left;
                elseif j < M - Discard_R_right - Guard_length_R && j > M - Discard_R_right - Guard_length_R - Ref_length_R % 太靠右，右侧Ref cell个数不足
                    R_mean_left = mean(R_AZ_pwr_heatmap(j - Guard_length_R - Ref_length_R : j - Guard_length_R - 1,i),1);
                    R_mean_right = R_mean_left;%mean(R_AZ_pwr_heatmap(j + Guard_length_R + 1:M - Discard_R_right,i),1);         
                else %居中，两侧ref + guard cell 个数充足
                    R_mean_left = mean(R_AZ_pwr_heatmap(j - Guard_length_R - Ref_length_R : j - Guard_length_R - 1,i),1);
                    R_mean_right = mean(R_AZ_pwr_heatmap(j + Guard_length_R + 1 : j + Guard_length_R + Ref_length_R,i),1);      
                end
                noise_R = min(R_mean_left,R_mean_right);
                if R_AZ_pwr_heatmap(j,i) > range_th * noise_R
                    %second pass CFAR
                    
                    if i <= 1 + Discard_AZ_left + Guard_length_AZ %太靠左，左侧无Ref cell 
                        AZ_mean_right = mean(R_AZ_pwr_heatmap(j,i + Guard_length_AZ + 1 : i + Guard_length_AZ + Ref_length_AZ),2);
                        AZ_mean_left = AZ_mean_right;
                    elseif i > 1 + Discard_AZ_left + Guard_length_AZ && i < 1 + Discard_AZ_left + Guard_length_AZ + Ref_length_AZ %太靠左，左侧Ref cell个数不足
                        AZ_mean_right = mean(R_AZ_pwr_heatmap(j,i + Guard_length_AZ + 1 : i + Guard_length_AZ + Ref_length_AZ),2);
                        AZ_mean_left = mean(R_AZ_pwr_heatmap(j,Discard_AZ_left + 1:i - Guard_length_AZ - 1),2);
                    elseif i >= N - Discard_AZ_right - Guard_length_AZ % 太靠右，右侧无Ref cell
                        AZ_mean_left = mean(R_AZ_pwr_heatmap(j,i - Guard_length_AZ - Ref_length_AZ : i - Guard_length_AZ - 1),2);
                        AZ_mean_right = AZ_mean_left;
                    elseif i < N - Discard_AZ_right - Guard_length_AZ && i > N - Discard_AZ_right - Guard_length_AZ - Ref_length_AZ % 太靠右，右侧Ref cell个数不足
                        AZ_mean_left = mean(R_AZ_pwr_heatmap(j,i - Guard_length_AZ - Ref_length_AZ : i - Guard_length_AZ - 1),2);
                        AZ_mean_right = mean(R_AZ_pwr_heatmap(j,i + Guard_length_AZ + 1:N - Discard_AZ_right),2);         
                    else %居中，两侧ref + guard cell 个数充足
                        AZ_mean_left = mean(R_AZ_pwr_heatmap(j,i - Guard_length_AZ - Ref_length_AZ : i - Guard_length_AZ - 1),2);
                        AZ_mean_right = mean(R_AZ_pwr_heatmap(j,i + Guard_length_AZ + 1 : i + Guard_length_AZ + Ref_length_AZ),2);      
                    end
                    maxpwr = max(abs(R_AZ_pwr_heatmap(j,:)));
                    noise_AZ = min(AZ_mean_left,AZ_mean_right);
                    if abs(R_AZ_pwr_heatmap(j,i)) > az_th * noise_AZ
                        detection_map(j,i) = 1;
                        detect_r_list = cat(1,detect_r_list,j);
                        detect_az_list = cat(1,detect_az_list,i);
                        SNR_list = cat(1,SNR_list,abs(R_AZ_pwr_heatmap(j,i))/noise_R);
                    else
                        if R_AZ_pwr_heatmap(j,i) > sidelobe_th * maxpwr && R_AZ_pwr_heatmap(j,i) > R_AZ_pwr_heatmap(j,i-1) && R_AZ_pwr_heatmap(j,i) > R_AZ_pwr_heatmap(j,i+1)  
                            detection_map(j,i) = 1;
                            detect_r_list = cat(1,detect_r_list,j);
                            detect_az_list = cat(1,detect_az_list,i);
                            SNR_list = cat(1,SNR_list,abs(R_AZ_pwr_heatmap(j,i))/noise_R);
                        end
                    end 
                end
            end
       end
%        figure(2)
%        mesh(detection_map)
       

       %% elevation est with capon BF
       ele_spectrum = zeros(Samples_per_Chirp,length(ele_bin));
       detection_map_3D = zeros(M,N,length(ele_bin));
       BF_weight = zeros(12,96);% 一个方向4个天线， 共96个距离点
       for range_idx = 1:length(detect_r_list)
           range_val = Range(detect_r_list(range_idx));
           for ChirpNum = 1:Num_chirp_loops
                %TX2,TX1 TX0 MIMO阵列 all 12 antennas
                Xk(:,:,ChirpNum) = cat(1,squeeze(DataRangeFft_Static_removed(detect_r_list(range_idx),1,ChirpNum,1:4)),squeeze(DataRangeFft_Static_removed(detect_r_list(range_idx),2,ChirpNum,1:4)),squeeze(DataRangeFft_Static_removed(detect_r_list(range_idx),3,ChirpNum,1:4)));
                Xk_sq(:,:,ChirpNum) = Xk(:,:,ChirpNum)*Xk(:,:,ChirpNum)';
           end
                R_xx_k(:,:,range_idx,Frame_ID) = mean(Xk_sq,3);
                 %diagonal loading is applied to the R matrix to ensure stability
                %DLF_k = 0.2*std(diag(R_xx_k(:,:,range_idx,Frame_ID)));
                DLF_k = 0.03;
                R_xx_k(:,:,range_idx,Frame_ID) = R_xx_k(:,:,range_idx,Frame_ID) + DLF_k*(trace(R_xx_k(:,:,range_idx,Frame_ID))/(Num_RX*Num_TX)).*eye(12);

            %1D_elevation_spectrum
            sv_az_el = sv_el .* sv_az(:,detect_az_list(range_idx));
            R_xx_k_inv = inv(R_xx_k(:,:,range_idx,Frame_ID));
            ele_spectrum = diag(1./(sv_az_el'*R_xx_k_inv*sv_az_el));
            %[ele_peak_tmp,ele_peak_idx] = findpeaks(abs(ele_spectrum));
            [ele_peak_tmp,ele_peak_idx] = max(abs(ele_spectrum));
            ele_left = ele_peak_idx - 1;
            if ele_left < 1
                ele_left = 1;
            end
            ele_right = ele_peak_idx + 1;
            if ele_right > length(ele_bin)
                ele_right = length(ele_bin);
            end
            maxIdxInterp = (ele_left * abs(ele_spectrum(ele_left)) + ele_peak_idx * ele_peak_tmp + ele_right * abs(ele_spectrum(ele_right)))/(abs(ele_spectrum(ele_left)) + ele_peak_tmp + abs(ele_spectrum(ele_right)));
            ele_Est = asind(ele_init + ele_step * maxIdxInterp); 
            tempAZ = (azi_init + detect_az_list(range_idx) * azi_step)/cosd(ele_Est);
            if abs(tempAZ) < 1
                azi_Est = asind(tempAZ);
            else
                continue
            end
%             figure(5)
%             plot(abs(ele_spectrum))
            %apply weights
            if length(ele_peak_tmp) == 1
                sv_az_el_max = sv_az_el(:,ele_peak_idx);
                BF_weight(:,detect_r_list(range_idx),detect_az_list(range_idx)) = R_xx_k_inv * sv_az_el_max;
                for chirp_cfg_idx = 1:Num_Chirp_config
                    weighted_cfg(:,chirp_cfg_idx,:) = squeeze(DataRangeFft_Static_removed(detect_r_list(range_idx),chirp_cfg_idx,:,1:4)) .* BF_weight(1 + 4*(chirp_cfg_idx-1):4 + 4*(chirp_cfg_idx-1),detect_r_list(range_idx),detect_az_list(range_idx)).';
                end
                for ch_idx = 1:4
                    for chirp_cfg_idx = 1:Num_Chirp_config
                        DataVelocityFft(:,detect_r_list(range_idx),detect_az_list(range_idx),4*(chirp_cfg_idx-1)+ch_idx) = fftshift(fft(squeeze(weighted_cfg(:,chirp_cfg_idx,ch_idx)),nfft_v,1))/(nfft_v/2);
                    end
                end
                DataVelocityFft_avg = mean(DataVelocityFft,4);
                [~,velocity_idx] = max(DataVelocityFft_avg(:,detect_r_list(range_idx),detect_az_list(range_idx)));
                velocity = RadicalSpeed(velocity_idx);
                SNR = SNR_list(range_idx);
                detected_points_info_per_frame = cat(1,detected_points_info_per_frame,[range_val,azi_Est,ele_Est,velocity,SNR,Frame_ID] );
%                 figure(6)
%                 plot(abs(DataVelocityFft(:,detect_r_list(range_idx),detect_az_list(range_idx),1)));
            end  
       end
       detected_points_info_per_frame = sort(detected_points_info_per_frame,1);
       detected_points_total = cat(1,detected_points_total,detected_points_info_per_frame);
end

       
       