% 绘制地震动反应谱及相应时程图
% Created on Tus Mar 8 15:00:00 2022
% @author: Vincent NT, ylpxdy@live.com
clc;clear;
% 示例：

%% 输入
recordFolder = '.\input';  % ### 文件夹
recordName = 'Oak_2_50_26_FN.acc';  % ### 文件名
% recordFolder = 'D:\Wen\Research\MAS\Duration\Chandramohan-Baker Database\1974 Lima, Peru - Arequipa\WF';  % 文件夹
% recordName = 'peru.c.PER01.058';  % 文件名

fmtType = 1;  % 0 = PEER; 1 = 1974 Lima, Peru from CESMD;
%% 读取地震动波形数据

if fmtType == 0
    %%%% PEER格式
    [wave, dt, NPTS, rsn] = getAmpDtPEER(recordFolder,recordName);
%     wave = wave./max(abs(wave));  % standardization of the waveform
    
elseif fmtType == 1
    %%%% 通用格式 for 1974 Lima, Peru from CESMD
    %%% 1）开头特定行（共HeaderLines行）记录地震和地震记录信息；
    %%% 2）数据部分以每行特定列记录（共DataCol列），且数据分隔符统一，数据位置统一；
    %%% 3）文件末尾仅以数据结束，并无其他结尾标记说明或符号。
    headerLines = 7;   % ### 表头行数
    dataCol = 5;    % ### 数据列数
    formatString = '%f %f %f %f %f';
%     formatString = '%*fS %f %f %f %f %f %f %f %f %f %fA %*f';  %%% 每行格式，参考网址：
    %%% https://ww2.mathworks.cn/help/matlab/ref/textscan.html?s_tid=doc_ta
    samplePoints = 2000;   % ### 采样的总点数
    dt = 0.02;  % ### 采样间隔时间，或 1/单位时间采样点数量
    scalarUnit = 1; % ### 单位变换因子，此处为CM/SEC**2（或gal）变为M/SEC**2
    
    wave = getAmpGeneral(recordFolder, recordName, ...
        formatString, headerLines, dataCol, samplePoints);
    
    wave(isnan(wave)) = [];
    wave = wave.*scalarUnit;
end

%% 计算反应谱及其他烈度指标

% 反应谱输出参数设置
abs_psd = 1; % ### Output abs (1), psd (0) or both (other values)
load("periodPEER.mat")
T = periodPEER'; % ### Natural period Interval
% fig = 0; % Default: plot the figure
% variable = 'A'; % Default: acceleration response spectrum
% normalize = 0; % Default: normalization ON

%%% Explore the impact of damping rate
kesi = 0.01:0.01:0.1;  % ### 阻尼比
% kesi = 0.05;

Sa = zeros(size(kesi,2),size(T,2));

for i = 1:1:size(kesi,2)
    % 计算反应谱
    [T, peak_abs, peak_psd] = responseSpectrum(...
        wave,dt,kesi(i),abs_psd,T,0,'A',0); 
    Sa(i,:) = peak_abs;
end

% 计算其他烈度指标（表），包括kesi=0.05的反应谱
units = 'g';  % ### 输入wave的单位
T1 = 0.60;    % ### 目标结构的第一自振周期
imTable = intensityCalculate(recordName, wave, dt, units, 0.05, T1, 0.05, T);

%% 绘制反应谱
% 绘图
figure
% 长持时波反应谱
p = plot(T, Sa, '-', 'LineWidth', 2.);

% 图设置
showPeriodStrat = 0;  % 显示区域的起始周期（s）
showPeriodEnd = 6;  % 显示区域的终止周期（s）
set(gca,'XLim',[showPeriodStrat showPeriodEnd])  % x轴显示范围
xlabel('\itT \rm(s)');  % x轴名
ylabel('\itSa \rm( \itT \rm{, 5%) (g)}');  % y轴名
xlim([0.01, 10])          % show range of T
set(gca,'XScale','log')   % show in log(T)
set(gca,'YScale','log')   % show in log(T)
grid on
lgd = legend(string(kesi),"Location","best");  % 图例
title(lgd, 'Damping ratio')

%% 绘制地震波图
% 绘制波形图的设置
selectPoint = [28.58, 24.46];  % 标记时间点（s）
gmFig = 1;   % 是否画地震动加速度，速度，位移图：1是，0否
SaFig = 0;   % 是否画反应谱：1是，0否
% 时间，加速度，速度，位移，谱加速度

% 绘制波形图
[time,~,vel,dsp,~] = plotGMandSa(wave,dt,selectPoint,gmFig,SaFig);

disp('Finish!')
