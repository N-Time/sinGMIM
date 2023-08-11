% ���Ƶ��𶯷�Ӧ�׼���Ӧʱ��ͼ
% Created on Tus Mar 8 15:00:00 2022
% @author: Vincent NT, ylpxdy@live.com
clc;clear;
% ʾ����

%% ����
recordFolder = '.\input';  % ### �ļ���
recordName = 'Oak_2_50_26_FN.acc';  % ### �ļ���
% recordFolder = 'D:\Wen\Research\MAS\Duration\Chandramohan-Baker Database\1974 Lima, Peru - Arequipa\WF';  % �ļ���
% recordName = 'peru.c.PER01.058';  % �ļ���

fmtType = 1;  % 0 = PEER; 1 = 1974 Lima, Peru from CESMD;
%% ��ȡ���𶯲�������

if fmtType == 0
    %%%% PEER��ʽ
    [wave, dt, NPTS, rsn] = getAmpDtPEER(recordFolder,recordName);
%     wave = wave./max(abs(wave));  % standardization of the waveform
    
elseif fmtType == 1
    %%%% ͨ�ø�ʽ for 1974 Lima, Peru from CESMD
    %%% 1����ͷ�ض��У���HeaderLines�У���¼����͵����¼��Ϣ��
    %%% 2�����ݲ�����ÿ���ض��м�¼����DataCol�У��������ݷָ���ͳһ������λ��ͳһ��
    %%% 3���ļ�ĩβ�������ݽ���������������β���˵������š�
    headerLines = 7;   % ### ��ͷ����
    dataCol = 5;    % ### ��������
    formatString = '%f %f %f %f %f';
%     formatString = '%*fS %f %f %f %f %f %f %f %f %f %fA %*f';  %%% ÿ�и�ʽ���ο���ַ��
    %%% https://ww2.mathworks.cn/help/matlab/ref/textscan.html?s_tid=doc_ta
    samplePoints = 2000;   % ### �������ܵ���
    dt = 0.02;  % ### �������ʱ�䣬�� 1/��λʱ�����������
    scalarUnit = 1; % ### ��λ�任���ӣ��˴�ΪCM/SEC**2����gal����ΪM/SEC**2
    
    wave = getAmpGeneral(recordFolder, recordName, ...
        formatString, headerLines, dataCol, samplePoints);
    
    wave(isnan(wave)) = [];
    wave = wave.*scalarUnit;
end

%% ���㷴Ӧ�׼������Ҷ�ָ��

% ��Ӧ�������������
abs_psd = 1; % ### Output abs (1), psd (0) or both (other values)
load("periodPEER.mat")
T = periodPEER'; % ### Natural period Interval
% fig = 0; % Default: plot the figure
% variable = 'A'; % Default: acceleration response spectrum
% normalize = 0; % Default: normalization ON

%%% Explore the impact of damping rate
kesi = 0.01:0.01:0.1;  % ### �����
% kesi = 0.05;

Sa = zeros(size(kesi,2),size(T,2));

for i = 1:1:size(kesi,2)
    % ���㷴Ӧ��
    [T, peak_abs, peak_psd] = responseSpectrum(...
        wave,dt,kesi(i),abs_psd,T,0,'A',0); 
    Sa(i,:) = peak_abs;
end

% ���������Ҷ�ָ�꣨��������kesi=0.05�ķ�Ӧ��
units = 'g';  % ### ����wave�ĵ�λ
T1 = 0.60;    % ### Ŀ��ṹ�ĵ�һ��������
imTable = intensityCalculate(recordName, wave, dt, units, 0.05, T1, 0.05, T);

%% ���Ʒ�Ӧ��
% ��ͼ
figure
% ����ʱ����Ӧ��
p = plot(T, Sa, '-', 'LineWidth', 2.);

% ͼ����
showPeriodStrat = 0;  % ��ʾ�������ʼ���ڣ�s��
showPeriodEnd = 6;  % ��ʾ�������ֹ���ڣ�s��
set(gca,'XLim',[showPeriodStrat showPeriodEnd])  % x����ʾ��Χ
xlabel('\itT \rm(s)');  % x����
ylabel('\itSa \rm( \itT \rm{, 5%) (g)}');  % y����
xlim([0.01, 10])          % show range of T
set(gca,'XScale','log')   % show in log(T)
set(gca,'YScale','log')   % show in log(T)
grid on
lgd = legend(string(kesi),"Location","best");  % ͼ��
title(lgd, 'Damping ratio')

%% ���Ƶ���ͼ
% ���Ʋ���ͼ������
selectPoint = [28.58, 24.46];  % ���ʱ��㣨s��
gmFig = 1;   % �Ƿ񻭵��𶯼��ٶȣ��ٶȣ�λ��ͼ��1�ǣ�0��
SaFig = 0;   % �Ƿ񻭷�Ӧ�ף�1�ǣ�0��
% ʱ�䣬���ٶȣ��ٶȣ�λ�ƣ��׼��ٶ�

% ���Ʋ���ͼ
[time,~,vel,dsp,~] = plotGMandSa(wave,dt,selectPoint,gmFig,SaFig);

disp('Finish!')
