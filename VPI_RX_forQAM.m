%VPI_RX_forQAM   Receiver side DSP for 16QAM fiber optical transmission
%system. Phase estimation is realized by blind phase search algorithm.
%
%   Version:            1.0
%   Author:             Yiming Wen
%   Contact:            yimingwen@hust.edu.cn
%   Last modified:      26/11/2025

%% Import Data
addpath(genpath('C:\Users\HP\Desktop\example\getstarted\MatlabCode\OptDSP_lite-master'));

close all
clear
y_Rx1=load('C:\Users\HP\Desktop\example\getstarted\MatlabCode\Outputs\sim001.mat');  
y_Rx=y_Rx1.y_Rx_X;

%% Parameter Settings
Fb=32e9;
Nss=2;
Fs=Nss*Fb;
N_length=51200;
DispSpectrum('Original signal',Fs,y_Rx)
M=64;
disp('--------------------------------------------------------');


%% Chromatic Dispersion Compensation
Dispersion = 16e-6;
fiberlength = 100e3;
SymbolRate = Fs;
N = 2;
y_Rx = CDC(y_Rx,Dispersion,fiberlength,SymbolRate,N);
y_Rx = y_Rx.';
figure;plot(real(y_Rx(:)),imag(y_Rx(:)),'b.');axis equal
title('receive after static equalization ')

%% DC Cancellation
y_Rx = y_Rx-mean(y_Rx);    %mean求平均
figure;plot(real(y_Rx(:)),imag(y_Rx(:)),'b.');axis equal
title('receive after DC cancellation ')

DispSpectrum('DC Cancellation',Fs,y_Rx);

%% IQ Imbanlance Compensation
[y_Rx]=GSOP_X(y_Rx);
figure;plot(real(y_Rx(:)),imag(y_Rx(:)),'b.');axis equal



%% Clock Recovery
% Not neccessary for VPI simulation, for there is only one frame of signal
% transmitted.

%% Match Filtering
rolloff = 0.2;   %滚降系数
span = 10;       %跨度
sps = Nss;       %每符号采样数
b = rcosdesign(rolloff,span,sps,"sqrt");     %设计升余弦或根升余弦滤波器,"sqrt"表示根升余弦
y_Rx = upfirdn(y_Rx,b);                      %MATLAB函数，实现上采样、滤波、下采样
y_Rx = y_Rx(span+1 : end-span);              %去除滤波器引入的延迟
figure;plot(real(y_Rx(:)),imag(y_Rx(:)),'b.');axis equal
title('receive after rcosine ')
grid on 
box on  

%% DownSampling
y_Rx = decimate3(y_Rx,Nss);
figure;plot(real(y_Rx(:)),imag(y_Rx(:)),'b.')
title('DownSampling Normalization');axis equal
grid on
box on   

%% Dynamic Equalization 
Taps =   31;
StepSize= 5e-6;
Convergence_Length =8000;
DN=0;
[y_Rx,indm] = RDE_Equalizer_blind(y_Rx,Taps,Convergence_Length,StepSize,Nss,DN,M);

y_Rx = QAM_normalization(y_Rx,M);
title('receive after dynamic equalization ')

%% Frequency Offset Estimation
Fd = 0.5*Fs;
[y_Rx_to,Frequency_offsetout]=FFT_Frequencyoffset_estimation(y_Rx,Fd);  %频偏补偿后的信号和估计出的频偏值
y_Rx=y_Rx_to;
y_Rx = QAM_normalization(y_Rx,M);

%% Phase Estimation
y_Rx=y_Rx(50:end);
B = 24;
N= 20;
[y_Rx,phase_est] = PhaseEstimate_BPS_QAM(y_Rx.',B,N,M)   ; 
figure; plot(-phase_est);
figure('color','w');plot(real(y_Rx(:)),imag(y_Rx(:)),'b.');axis equal
grid on 
box on   

%% DD-LMS
Taps = 13;
Length = 14000;
StepSize1 = 20e-5;
StepSize2 = 20e-5;
[RI,RQ] = IQ_DD_LMS_Equalizer(real(y_Rx),imag(y_Rx),Taps,Length,StepSize1,StepSize2,M);

y_Rx = RI+1i*RQ;
y_Rx = QAM_normalization(y_Rx,M);



scatterPlot(y_Rx);
title('IQ Constellation of the Transmitted Signal with Oversampling')

%% ====== PAS( CCDM ) + LDPC Soft Demap & Decode ======
% tx_dump.mat 已经由发端保存，且包含 txDump 结构体
dumpPath = 'C:\Users\HP\Desktop\example\getstarted\MatlabCode\Outputs\tx_dump.mat';  
load(dumpPath, 'txDump');

% --- 基本检查 ---
if txDump.M ~= M
    warning('TX M=%d but RX M=%d, please check.', txDump.M, M);
end

% --- EVM / SNR（用于估计 N0） ---
[EVM, SNR] = EVM_Measurement(M, y_Rx);
SNRlin = 10^(SNR/10);

% 注意：你的 QAM_normalization 不是 UnitAveragePower=1，而是拉到 sqrt(42) 尺度
% 所以用当前符号能量 Es 来保证 N0 尺度一致
Es = mean(abs(y_Rx).^2);
N0 = Es / SNRlin;

fprintf('[RX] EVM=%.3f%%, SNR=%.2f dB, Es=%.3f, N0=%.3e\n', EVM*100, SNR, Es, N0);

% --- 生成 PAS 的 2D 星座先验概率 symProb2D ---
% 优先从 PCS 里取 1D 幅度概率
C = txDump.C;  % qammod(0:63,64,'gray').'
if isfield(txDump.PCS, 'symProb1D') && isfield(txDump.PCS, 'C_I')
    symProb1D = txDump.PCS.symProb1D(:);
    C_I = txDump.PCS.C_I(:);
    symProb2D = pas_makeSymProb2D_64QAM(C, C_I, symProb1D);
else
    warning('PCS.symProb1D / PCS.C_I not found, using uniform priors.');
    symProb2D = [];
end

% --- 计算 MAP LLR（Gray + 可选非均匀先验）---
% 输出 LLRs 为 [nSyms*log2(M) x 1]，按“每符号连续6比特”排列
LLRs = LLR_eval_gray_ps(y_Rx, N0, C, symProb2D);

% --- 准备LDPC解码参数 ---
PCM_FEC = txDump.PCM_FEC;          % 或者 dvbs2ldpc(txDump.R_FEC)
idx_FEC = txDump.PCS.FEC.idx;      % 关键索引
nCols = size(PCM_FEC, 2);
nIter_FEC = 50;                    % 可调

% --- 可选：FEC 同步（强烈建议你保持开启，因为你前面丢了 y_Rx(50:end) 等操作）---
doSync = true;

if doSync
    % 同步参考比特：优先用 PCS.txBits_trunc（你在 Tx_PS_CCDM_PAS 截断后保存了）
    if isfield(txDump.PCS, 'txBits_trunc')
        txBits_ref = txDump.PCS.txBits_trunc;
    else
        % 兜底：用输入比特（可能比实际短/长，仿真不一定完全对齐）
        warning('PCS.txBits_trunc not found, using txDump.txBits_in as sync reference (may be imperfect).');
        txBits_ref = txDump.txBits_in;
    end

    % FEC_syncBits 典型输出：对齐后的 LLR、对齐后的 idx_FEC、FEC block 数等
    % 如果你的 FEC_syncBits 输出参数名字不同，按报错改一下即可
[LLRs_sync, idx_FEC_sync, nFEC, bestShift] = FEC_syncBits_simple(txBits_ref, LLRs, idx_FEC, nCols);
fprintf('[RX] Sync bestShift=%d bits\n', bestShift);

    LLRs_use = LLRs_sync; 
    idx_use  = idx_FEC_sync;

    fprintf('[RX] FEC sync done. nFEC=%d\n', nFEC);
else
    LLRs_use = LLRs.';   % LDPC_decoder 可能希望行向量
    idx_use  = idx_FEC;
end

% --- LDPC 解码 ---
% 你的 LDPC_decoder 期望 LLRs 的形状可能是列向量或行向量，这里统一给列向量更稳
N = size(PCM_FEC,2);

LLRvec = LLRs_use(:);
nBlocks = floor(numel(LLRvec)/N);
LLRvec = LLRvec(1:nBlocks*N);

% 如果你同步后通常只剩 1 个码字，就用单码字 idx（最稳）
idx_FEC_one = unique(mod(txDump.PCS.FEC.idx(:)-1, N) + 1, 'stable');

[rxBits_info, nBlocks_out] = LDPC_decoder(LLRvec(1:N), PCM_FEC, idx_FEC_one, nIter_FEC);
fprintf('[RX] LDPC decoded blocks = %d (used 1 block)\n', nBlocks_out);

% --- BER 统计（与 Tx 的参考信息比特对齐）---
% 最正确：用 txBits_trunc 的前 numel(rxBits_info) 位作为参考
if isfield(txDump.PCS, 'txBits_trunc')
    txBits_true = txDump.PCS.txBits_trunc(:);
else
    txBits_true = txDump.txBits_in(:);
end

disp('[RX] LDPC decode finished. (Q/BER print disabled)');
disp('--------------------------------------------------------');

fprintf('bestShift mod nCols = %d\n', mod(bestShift, nCols));
