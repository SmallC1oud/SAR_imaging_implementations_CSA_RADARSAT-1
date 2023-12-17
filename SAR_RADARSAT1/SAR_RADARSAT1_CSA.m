%% 合成孔径雷达 星载SAR数据成像---RADARSAT-1数据处理之CSA成像
% 成像数据在 CD_run_params.mat 中有储存
% 使用CDdata1.mat 的数据进行成像
% 成像流程：
% 1、原始数据(绘图)
% 2、方位向FFT，变换到距离多普勒域(RD域)，进行"补余RCMC"
% 3、距离向FFT，变换到二维频域，进行"距离压缩"、"SRC"、"一致RCMC"
% 4、距离向IFFT，变换到多普勒频域，进行"方位压缩"和"附加相位校正"
% 5、方位向IFFT，回到图像域，绘制成像结果
clc;clear;close all
% ---------------------------------------------------------------------
%% 参数设置
% 雷达参数
C = physconst('LightSpeed');                                                 % 光速
% C = 3e8;
f_0 = 5.3e9;                                                                 % 雷达中心频率
lambda = C/f_0;                                                              % 波长
V = 7062;                                                                    % 平台速度--等效雷达速度
% 距离向
K_r = -0.72135e12;                                                           % 距离调频率
T_p = 30.111e-6;                                                             % 脉冲宽度
T_r = 41.75e-6;                                                              % 脉冲持续时间
F_r = 32.317e6;                                                              % 距离向采样率
d_tau = 1/F_r;                                                               % 距离向采样时间间隔
% 方位向
K_a = 1733;                                                                  % 方位向调频率
F_a = 1256.98;                                                               % 方位向采样率--脉冲重复频率
f_eta_c = -6900;                                                             % 多普勒中心频率
d_eta = 1/F_a;                                                               % 距离向采样时间间隔
% 其他参数
t_0 = 6.5956e-3;                                                             % 获取数据时间--数据窗开始时间
R_0 = t_0*C/2;                                                               % 最短斜距
%% 获取原始数据
% load('CDdata1.mat')
SAR_data = importdata('S:\MATLAB_code_2023\SAR_2023\SAR_exercise_3\RadarsatData\CD_Data_Processing\CDdata1.mat');
% SAR_data = importdata('S:\MATLAB_code_2023\SAR_2023\SAR_exercise_3\RadarsatData\CD_Data_Processing\CDdata2.mat');
%% 绘制原始回波数据
figure(1)
imagesc(abs(double(SAR_data)));
% title('原始数据'); %未补零

SAR_data=double(SAR_data);
[N_a,N_r] = size(SAR_data);

%% 参考数据设置
R_ref = (t_0+N_r/2*d_tau)*C/2;                                               % 参考距离
V_ref = V;                                                                   % 参考速度
f_ref = f_eta_c;                                                             % 参考频率

%% 数据末尾补零
Z_a = 800+1760;                                                              % 方位向补零数
Z_r = ceil(T_r/d_tau)+698;                                                   % 距离向补零数
SAR_data = cat(2,SAR_data,zeros(N_a,Z_r));                                   % 距离向补零
SAR_data = cat(1,zeros(800,N_r+Z_r),SAR_data);
SAR_data = cat(1,SAR_data,zeros(1760,N_r+Z_r));                              % 方位向补零
N_a = N_a+Z_a;%+800=2336
N_r = N_r+Z_r;%+1350=3398

figure(2),imagesc(abs(SAR_data));axis image;set(gcf,'Color','w');
% title('时域---补零后的原始信号幅度');
xlabel('距离向（采样点）');ylabel('方位向（采样点）');
% imwrite(uint8(abs(SAR_data)),'absData.bmp','bmp');
%% 时间轴、频率轴设置
tau = t_0+(0:N_r-1)*d_tau;                                                   % 距离向时间轴
f_tau = ((0:N_r-1)-N_r/2)/N_r*F_r;                                           % 距离向频率轴
eta = ((0:N_a-1)-N_a/2)*d_eta;                                               % 方位向时间轴
f_eta = f_eta_c+((0:N_a-1)-N_a/2)/N_a*F_a;                                   % 方位向频率轴

%% 中间变量设置
D_feta_V = sqrt(1-C^2*f_eta.^2/(4*V^2*f_0^2));                               % D(feta,V)，式7.17
D_feta_Vref = sqrt(1-C^2*f_eta.^2/(4*V_ref^2*f_0^2));                        % 双曲距离方程：D(feta,Vref)，式7.17
D_fref_V = sqrt(1-C^2*f_ref^2/(4*V^2*f_0^2));                                % D(fref,V)，式7.17
D_fref_Vref = sqrt(1-C^2*f_ref^2/(4*V_ref^2*f_0^2));                         % D(fref,Vref)，式7.17

K_src = (2*V^2*f_0^3*D_feta_V.^3)./(C*R_0*f_eta.^2);                         % SRC滤波器的调频率，式6.22        
K_m = K_r./(1-K_r./K_src);                                                   % 距离多普勒域：改变后的距离向调频率K_m，式6.21
RCM_bulk_feta = (1./D_feta_Vref-1/D_fref_Vref)*R_ref;                        % 一致RCM，，式7.22
alpha = D_fref_Vref./D_feta_Vref-1;                                          % 频偏参数，式7.28

%% 对齐到方位向多普勒中心
S0 = SAR_data.*exp(-1i*2*pi*f_eta_c*(eta'*ones(1,N_r)));                     % 搬移至多普勒频率中心

%% 变换到距离多普勒域，实现变标相乘
Srd = fftshift(fft(fftshift(S0,1),[],1),1);                                  % 原信号做方位向FFT

tt = 2/C*(R_0/D_fref_V+RCM_bulk_feta)-2*R_ref./(C*D_feta_Vref);              % P205 (7.26) (7.27)
Ssc = exp(1i*pi*K_m.*alpha.*tt.^2);                                          % 线性调频变标方程变标方程 P207 (7.30)
S1 = Srd.*(Ssc'*ones(1,N_r));                                                % 变标相乘，，式7.31

%% 变换到二维频域，实现RC、SRC和一致RCMC
S2 = fftshift(fft(fftshift(S1,2),[],2),2);                                   % 信号变换到二维频域 

WindowR = ones(N_a,1)*kaiser(N_r,2.5)';                                      % 距离窗
WindowA = kaiser(N_a,2.5)*ones(1,N_r);                                       % 方位窗
S2 = S2.*WindowR.*WindowA;                                                   % 加窗，式7.32
Hm = exp(1i*pi./((K_m'*ones(1,N_r)).*(1+alpha'*ones(1,N_r))).*(ones(N_a,1)*f_tau).^2); % 合并的距离压缩和一致RCMC滤波器，，式7.32,1+alpha
Hrcm = exp(1i*4*pi/C*(RCM_bulk_feta'*ones(1,N_r)).*(ones(N_a,1)*f_tau));             % 一致RCM，式7.32，式7.22
S3 = S2.*Hm.*Hrcm;                                                           % 相位相乘，式7.33

%% 变换到距离多普勒域，实现方位压缩和附加相位校正
S4 = ifftshift(ifft(ifftshift(S3,2),[],2),2);                                % 距离向IFFT，式7.34

figure(3),imagesc(abs(S4));axis image;set(gcf,'Color','w');
% title('RD域：RC、SRC和一致RCMC后的信号幅度');
xlabel('距离向（采样点）');ylabel('方位向（采样点）');

Hac = exp(-1i*4*pi*R_0*f_0*D_feta_V/C);                                       % 方位匹配滤波器，式7.34第一个指数项
Hpc = exp(1i*(4*pi*K_m/C^2).*(1-D_feta_Vref/D_fref_Vref).*(R_0./D_feta_V-R_ref./D_feta_V).^2); % 相位校正滤波器，式7.34第二个指数项
S5 = S4.*(Hac'*ones(1,N_r)).*(Hpc'*ones(1,N_r));                              % 相位相乘，实现滤波，或者取H复共轭相消

%% 变换到图像域
Stt = ifftshift(ifft(ifftshift(S5,1),[],1),1);                                % 信号的方位向IFFT
%% 调整图像
Img = flipud(abs(Stt));                                                       % 上下翻转图像
Ntau = round((T_r/d_tau/2-1)/2);                                              % 计算距离向弃置区长度
Nta = round((800/2-1)/2);
Img = Img(Z_a-Nta+1:N_a-Nta,Ntau+1:Ntau+N_r-Z_r);                             % 裁剪图像有效区域
Img = Img/max(max(Img));
Img = 20*log10(Img+eps);
Img(Img<-60) = -60;
Img = uint8((Img+60)/50*255);
figure,imagesc(Img);axis image;set(gcf,'Color','w');
figure,imshow(Img);
%===========================================================================================================================
%===========================================================================================================================
% 9704*8192
%% 获取原始数据
% load('CDdata1.mat')
SAR_data = importdata('S:\MATLAB_code_2023\SAR_2023\SAR_exercise_3\RadarsatData\CD_Data_Processing\CDdata1_9704_8192.mat');
% SAR_data = importdata('S:\MATLAB_code_2023\SAR_2023\SAR_exercise_3\RadarsatData\CD_Data_Processing\CDdata2_9704_8192.mat');
%% 绘制原始回波数据
figure
imagesc(abs(double(SAR_data)));
% title('原始数据'); %未补零

SAR_data=double(SAR_data);
[N_a,N_r] = size(SAR_data);

%% 参考数据设置
R_ref = (t_0+N_r/2*d_tau)*C/2;                                               % 参考距离
V_ref = V;                                                                   % 参考速度
f_ref = f_eta_c;                                                             % 参考频率
Z_a = 0;                                                                     % 方位向补零数
Z_r = 0;                                                                     % 距离向补零数

figure,imagesc(abs(SAR_data));axis image;set(gcf,'Color','w');
% title('时域---补零后的原始信号幅度');
xlabel('距离向（采样点）');ylabel('方位向（采样点）');
%% 时间轴、频率轴设置
tau = t_0+(0:N_r-1)*d_tau;                                                   % 距离向时间轴
f_tau = ((0:N_r-1)-N_r/2)/N_r*F_r;                                           % 距离向频率轴
eta = ((0:N_a-1)-N_a/2)*d_eta;                                               % 方位向时间轴
f_eta = f_eta_c+((0:N_a-1)-N_a/2)/N_a*F_a;                                   % 方位向频率轴

%% 中间变量设置
D_feta_V = sqrt(1-C^2*f_eta.^2/(4*V^2*f_0^2));                               % D(feta,V)，式7.17
D_feta_Vref = sqrt(1-C^2*f_eta.^2/(4*V_ref^2*f_0^2));                        % 双曲距离方程：D(feta,Vref)，式7.17
D_fref_V = sqrt(1-C^2*f_ref^2/(4*V^2*f_0^2));                                % D(fref,V)，式7.17
D_fref_Vref = sqrt(1-C^2*f_ref^2/(4*V_ref^2*f_0^2));                         % D(fref,Vref)，式7.17

K_src = (2*V^2*f_0^3*D_feta_V.^3)./(C*R_0*f_eta.^2);                         % SRC滤波器的调频率，式6.22        
K_m = K_r./(1-K_r./K_src);                                                   % 距离多普勒域：改变后的距离向调频率K_m，式6.21
RCM_bulk_feta = (1./D_feta_Vref-1/D_fref_Vref)*R_ref;                        % 一致RCM，，式7.22
alpha = D_fref_Vref./D_feta_Vref-1;                                          % 频偏参数，式7.28

%% 对齐到方位向多普勒中心
S0 = SAR_data.*exp(-1i*2*pi*f_eta_c*(eta'*ones(1,N_r)));                     % 搬移至多普勒频率中心

%% 变换到距离多普勒域，实现变标相乘
Srd = fftshift(fft(fftshift(S0,1),[],1),1);                                  % 原信号做方位向FFT

tt = 2/C*(R_0/D_fref_V+RCM_bulk_feta)-2*R_ref./(C*D_feta_Vref);              % P205 (7.26) (7.27)
Ssc = exp(1i*pi*K_m.*alpha.*tt.^2);                                          % 线性调频变标方程变标方程 P207 (7.30)
S1 = Srd.*(Ssc'*ones(1,N_r));                                                % 变标相乘，，式7.31

%% 变换到二维频域，实现RC、SRC和一致RCMC
S2 = fftshift(fft(fftshift(S1,2),[],2),2);                                   % 信号变换到二维频域 

WindowR = ones(N_a,1)*kaiser(N_r,2.5)';                                      % 距离窗
WindowA = kaiser(N_a,2.5)*ones(1,N_r);                                       % 方位窗
S2 = S2.*WindowR.*WindowA;                                                   % 加窗，式7.32
Hm = exp(1i*pi./((K_m'*ones(1,N_r)).*(1+alpha'*ones(1,N_r))).*(ones(N_a,1)*f_tau).^2); % 合并的距离压缩和一致RCMC滤波器，，式7.32,1+alpha
Hrcm = exp(1i*4*pi/C*(RCM_bulk_feta'*ones(1,N_r)).*(ones(N_a,1)*f_tau));     % 一致RCM，式7.32，式7.22
S3 = S2.*Hm.*Hrcm;                                                           % 相位相乘，式7.33

%% 变换到距离多普勒域，实现方位压缩和附加相位校正
S4 = ifftshift(ifft(ifftshift(S3,2),[],2),2);                                % 距离向IFFT，式7.34
figure,imagesc(abs(S4));axis image;set(gcf,'Color','w');
% title('RD域：RC、SRC和一致RCMC后的信号幅度');
xlabel('距离向（采样点）');ylabel('方位向（采样点）');

Hac = exp(-1i*4*pi*R_0*f_0*D_feta_V/C);                                       % 方位匹配滤波器，式7.34第一个指数项
Hpc = exp(1i*(4*pi*K_m/C^2).*(1-D_feta_Vref/D_fref_Vref).*(R_0./D_feta_V-R_ref./D_feta_V).^2); % 相位校正滤波器，式7.34第二个指数项
S5 = S4.*(Hac'*ones(1,N_r)).*(Hpc'*ones(1,N_r));                              % 相位相乘，实现滤波，或者取H复共轭相消
%% 变换到图像域
% Stt = ifftshift(ifft(ifftshift(S5,1),[],1),1);                              
Stt = ifft(ifftshift(S5,1),[],1);                                             % 信号的方位向IFFT
%% 调整图像
Img = flipud(abs(Stt));                                                       % 上下翻转图像
Img = Img/max(max(Img));
Img = 20*log10(Img+eps);
Img(Img<-60) = -60;
Img = uint8((Img+60)/50*255);
figure,imagesc(Img);axis image;set(gcf,'Color','w');
figure,imshow(Img);
%==========================================================================