function [xcoords,ycoords] = loadChanMap(chanMap_folder,label)
% waveformMean: cluster *  tSampleN * nChan (384)
% load channel map and scale it to a matrix of size row *col 
% based on chanMap
% wfIm2 : cluster * row * col * tSampleN
% col: max(xc); row: max(yc)
switch label
    case 'NPUHD2'
        chanMap = 'NPUHD2_bank0_ref0.mat';        
    case 'NP10'
        chanMap = 'neuropixPhase3B1_kilosortChanMap.mat';
    case 'NP20'
        chanMap = 'NPtype24_doubleLengthStripe_botRow0_ref0.mat';
end
load(fullfile(chanMap_folder,chanMap));