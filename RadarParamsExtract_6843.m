function [ADCSample, ChirpPeriod, ADCFs, nchirp_loops, FramPeriod, FramNum, slope, BandWidth, R_Maximum, R_resulo, V_Maximum, V_resulo] =  RadarParamsExtract_6843(addr,rawfileName,rfFreqScaleFactor,Num_Chirp_config,Num_TX)
%%  radar raw data choose
% addr = 'E:\ZJL_Sync\Cascade_Radar_Õ‚≥°≤‚ ‘\Dropping Objects\FallingDetecting\RadarData';   % File address
cd(addr)
% rawfileName = 'BBTX1BaiDuCanFallingSSY1.bin'; %Fc400BPSKb200t05M0_NAT1_V08W0S2 Change only this line   rawData3D_simple2D

%% Parameter load
% addpath('E:\TarfProjectFunc')   %addpath
paraFileName = [rawfileName(1:end-4) '_LogFile.txt'] ;
delimiterIn = ' ';
para = importdata(paraFileName,delimiterIn) ;
s1 = strfind(para,'API:ProfileConfig');
idx1 = find(cellfun(@(x)~isempty(x),s1,'UniformOutput',true));
s2 = strfind(para,'API:FrameConfig');
idx2 = find(cellfun(@(x)~isempty(x),s2,'UniformOutput',true));
pc = strsplit(para{idx1(end),1},',');
pf = strsplit(para{idx2(end),1},',');

ADCSample = str2double(pc{1,11}) ;%number
ChirpPeriod = (str2double(pc{1,4})+str2double(pc{1,6}))/100 ;  %us
ADCFs = str2double(pc{1,12})*1000  ; % sps
nchirp_loops = str2double(pf{1,5}) ;
FramPeriod = str2double(pf{1,6})/1e2/2;  % us
FramNum  =   str2double(pf{1,4})     ;
%slope =  (str2double(pc{1,9})+2)*0.048*1e12   ;  %
slope = ((str2double(pc{1,9})*rfFreqScaleFactor*1e3*900)/(2^26))*1e12   ;  %

fs = 1/ (ChirpPeriod*Num_Chirp_config*1e-6);

BandWidth = ADCSample/ADCFs*slope;
lamda = 3e8/60.75e9;

R_resulo =3e8/2/BandWidth; % range resolution
V_resulo =  lamda/(2*ChirpPeriod*Num_Chirp_config*1e-6*nchirp_loops*Num_Chirp_config)  ;

%% Display
R_Maximum = ADCFs*0.9*3e8/(2*slope);
% R_Maximum = ADCFs/slope*3e8/2;
V_Maximum = lamda/4*fs;
sprintf('Maximum Range= %.2f (m)\nRange Resolution= %.3f (m)\nMaximum Velocity= %.2f (m)\nVelocity Resolution= %.3f (m)',...
    R_Maximum,R_resulo,V_Maximum,V_resulo)

end


