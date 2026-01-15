function [txI,txQ] = VPI_QAM_TX(y_in)
%VPI_QAM_TX   Prepare data and modulate into QAM format for VPI simulation, working with Tx CoSimInterface.vtms

%% Transmission Parameters
close all

% Load Libraries:
addpath(genpath('C:\Users\HP\Desktop\example\getstarted\MatlabCode\OptDSP_lite-master'));

RGB = fancyColors();

SIG.M = 64;                 % QAM constellation size
SIG.symRate = 32e9;         % total symbol-rate of the signal
SIG.bitRate_net = 100e9;    % net bit-rate [Gbps]
SIG.modulation = 'QAM';     % modulation type [QAM/PAM]
SIG.rollOff = 0.2;          % roll-off factor
SIG.nPol = 1;               % number of polarizations
% SIG.nSyms = 25600;          % total number of simulated symbols
nSpS = 2;                   % number of samples per symbol
SIG.nSyms = floor(length(y_in.band.E) / nSpS);
R_FEC = 5/6;             % FEC rate


pilotRate = 1;              % Pilot rate
useCPE2 = false;            % flag to decide whether to use 2nd CPE or not


%nBpS_net = SIG.bitRate_net / (SIG.nPol*SIG.symRate*FEC_rate*pilotRate);
nBpS_net = SIG.bitRate_net / (SIG.nPol*SIG.symRate*pilotRate);
TX.SIG = setSignalParams('symRate',SIG.symRate,'M',SIG.M,...
    'nPol',SIG.nPol,'nBpS',nBpS_net,'nSyms',SIG.nSyms,...
    'roll-off',SIG.rollOff,'modulation',SIG.modulation);

% Modulation Parameters:
TX.QAM = QAM_config(TX.SIG);

% Bit Parameters:
TX.BIT.source = 'randi';
TX.BIT.seed = 1;

% Pulse Shaping Filter Parameters:
TX.PS.type = 'RRC';
TX.PS.rollOff = TX.SIG.rollOff;
TX.PS.nTaps = 128;

% DAC Parameters:
TX.DAC.RESAMP.sampRate = nSpS*TX.SIG.symRate;

% DSP Pilots Parameters:
TX.PILOTS.active = true;
TX.PILOTS.rate = pilotRate;
TX.PILOTS.option = 'outerQPSK';
% TX.PILOTS.scaleFactor = 1;

% FEC Parameters:
% TX.FEC.active = false;
% TX.FEC.rate = FEC_rate;
% TX.FEC.nIter = 50;

% PCS Parameters:
TX.PCS.method = 'CCDM';

% Generate Tx Bits
TX.BIT.txBits = Tx_generateBits(SIG.nSyms,TX.QAM.M,TX.QAM.nPol,TX.BIT);

%% Generate Transmitted Symbols
% [S.tx,txSyms,TX.QAM] = Tx_ProbShaping(TX.BIT.txBits,TX.QAM,TX.SIG,TX.FEC.rate);
C = qammod(0:63,64,'gray').';
H_DM = 3.5;
TX.PCS.blockLength = 1800;

[S.tx,txSyms,TX.PCS] = Tx_PS_CCDM_PAS(C,SIG.M,H_DM,R_FEC,TX.PCS,TX.BIT.txBits);   %Iutput H

%% ===== 保存给接收端用的关键数据（结构体 txDump）=====
savePath = 'C:\Users\HP\Desktop\example\getstarted\MatlabCode\Outputs\tx_dump.mat';  

txDump = struct();

% ===== 基本体制参数 =====
txDump.M       = SIG.M;
txDump.nBpS    = log2(SIG.M);
txDump.symRate = SIG.symRate;
txDump.rollOff = SIG.rollOff;
txDump.nSpS    = nSpS;

% ===== PAS / CCDM 参数 =====
txDump.H_DM    = H_DM;
txDump.R_FEC   = R_FEC;
txDump.C       = C;              % qammod(0:63,64,'gray').'
txDump.PCS     = TX.PCS;         % 含 PCS.FEC.idx、symProb1D/C_I/txBits_trunc/shapedBits 等

% ===== LDPC 相关 =====
txDump.PCM_FEC = dvbs2ldpc(R_FEC);

% ===== 参考比特（输入到 PAS 的原始比特；截断后的比特建议在 PCS.txBits_trunc 里）=====
txDump.txBits_in = TX.BIT.txBits;

% ===== 对齐/排错信息 =====
txDump.VPI_LEN  = length(y_in.band.E);
txDump.nSyms    = SIG.nSyms;
txDump.scaleIQ  = 0.1;           % 重要：你 TX 输出给 VPI 的缩放

% 保存发送符号索引/符号（不大，排错有用）
txDump.txSyms   = txSyms;
txDump.Stx      = S.tx;          % 1sps 符号（脉冲成形前）

% ===== 版本信息=====
txDump.timestamp = datestr(now, 30);
txDump.note = 'PAS(CCDM)+DVBS2-LDPC TX dump for VPI co-sim';

save(savePath, 'txDump', '-v7.3');
fprintf('[TX] Saved txDump to: %s\n', savePath);


%% Plot Transmitted Constellation Symbols:
plotProbShaping_PDF_const('const',TX.QAM.IQmap,'symbols',txSyms(:),'color',RGB.itred);

% Plot Transmitted Constellation Symbols:
scatterPlot(S.tx,'markersize:sig',20);
title('IQ Constellation of the Transmitted Symbols','interp','latex')

% Pulse Shaping
[S.txSC,TX.PS] = pulseShaper(S.tx,nSpS,TX.PS);

% Plot Pulse Shaping Taps:
figure();
plot(TX.PS.W);
title('Pulse Shaping Taps')

% Plot Signal Spectrum:
figure();
pwelch(S.txSC(1,:),1e4,[],[],TX.DAC.RESAMP.sampRate,'centered')
% Plot Time-Domain Signal:
figure();
ni = 50;
nf = 100;
plot(real(S.txSC(1,ni:nf)));
hold on;
plot(imag(S.txSC(1,ni:nf)));
hold off;
legend('I','Q')
title('Time-Domain Signal')
% Plot Time-Domain Signal Constellation with Oversampling:
scatterPlot(S.txSC(1,:));
title('IQ Constellation of the Transmitted Signal with Oversampling')

% Plot Time-Domain Signal Constellation after Downsampling
scatterPlot(S.txSC(1,1:nSpS:end),'markersize:sig',20);
title('IQ Constellation of the Transmitted Signal after Downsampling')


%% Assign to VPI
I_full = real(S.txSC(1,:)) * 0.1;
Q_full = imag(S.txSC(1,:)) * 0.1;

VPI_LEN = length(y_in.band.E);   % = 43200

% ===== 强制对齐 VPI 帧长 =====
if length(I_full) >= VPI_LEN
    I = I_full(1:VPI_LEN);
    Q = Q_full(1:VPI_LEN);
else
    I = [I_full zeros(1,VPI_LEN-length(I_full))];
    Q = [Q_full zeros(1,VPI_LEN-length(Q_full))];
end



% ---------------------- 新增3：传给VPI的最终I/Q星座图 ----------------------
figure('Name','传给VPI的最终I/Q星座');
plot(I,Q,'g.','MarkerSize',3);
grid on; axis equal;
xlabel('Tx I 分量'); ylabel('Tx Q 分量');
title('VPI发送端最终I/Q星座（幅度×0.1）');
% ----------------------------------------------------------------------------
txI = y_in;
txQ = y_in;

for i = 1:length(y_in.band.E) 
    txI.band.E(i) = I(i);
    txQ.band.E(i) = Q(i);
end

fprintf('PAS symbols = %d\n', length(txSyms));
fprintf('VPI symbols = %d\n', length(y_in.band.E)/nSpS);

end


