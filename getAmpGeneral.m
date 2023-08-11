% Func: Get a GM Wave amplication from the data file in a general format.
% Created on Tus Mar 8 15:00:00 2022
% @author: Vincent NT, ylpxdy@live.com


% % Input:
% % folder name
%
% % How to build a correct formatSpec in MATLAB?
% % https://ww2.mathworks.cn/help/matlab/ref/textscan.html?s_tid=doc_ta
%
% % 1) Data table without row variable or row note:
% recordFolder = 'D:\Wen\Research\MAS\Duration\Chandramohan-Baker Database\1985 Valparaiso, Chile - Llolleo\WF';
% recordName = '1985.c.062w47ll.o0a';  % file name
% %%%% FORMAT: '%f %f %f %f %f %f %f %f %f %f'
% dataCol = 10;  % data column
% formatString = repmat('%f ',1,dataCol);  % read format in a line
% samplePoints = 23279;   % sampling points
% headerLines = 7;   % header lines
%
% % 2) Data table with row variable (at the beginning column) or note (at last
% % column):
% recordFolder = 'D:\Wen\Research\MAS\Duration\Chandramohan-Baker Database\1974 Lima, Peru - Arequipa\WF';
% recordName = 'peru.c.PER01.058';
% %%%% FORMAT: '%*fS %f %f %f %f %f %f %f %f %f %fA %*f'
% dataCol = 10;  % data columns
% formatStringAmp = repmat('%f ',1,dataCol-1);  % read format in a line
% formatStringFirst = '%*fS';
% formatStringLast = '%fA %*f';
% formatString = [formatStringFirst,' ',formatStringAmp, formatStringLast]; % build format
% samplePoints = 1972;   % sampling points
% headerLines = 4;   % header lines
% 
% % e.g.
% data = getAmpGeneral(recordFolder,recordName,...
%     formatString,headerLines,dataCol,samplePoints);

function dataVector = getAmpGeneral(recordFolder, recordName, ...
    formatString, headerLines, dataCol, samplePoints)

    fid = fopen([recordFolder '\' recordName], 'r');  % open data file
    
    c = textscan(fid,formatString,(ceil(samplePoints/dataCol)),...
        'headerlines',headerLines,'EmptyValue',NaN);
    
    % supply the empty data (no space) at the final row with NaN
    colNum = zeros(1,size(c,2));   % find the max column number
    for i = 1:1:size(c,2)
        colNum(i) = size(c{i},1);
    end
    colNumMax = max(colNum);
    colNumMin = min(colNum);
    
    if colNumMax ~= colNumMin   % judge the inconsistent row length
        cNew = cell(1,size(c,2));   % to store the supplied nan data 
        for i = 1:1:size(c,2)
            cNew{1,i} = [c{i}; NaN(colNumMax-size(c{i},1))];  % supply NaN
        end
    else
        cNew = c;
    end

    data = cell2mat(cNew);  % read data

    % trasform data matrix to a wave vector in n*1
    dataVector = data';
    dataVector = dataVector(:);

    fclose(fid);

end