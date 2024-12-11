%% Add paths 
githubDir = 'C:\Users\Steinmetz lab\Documents\git';
addpath(genpath(fullfile(githubDir, 'npy-matlab')))
addpath(genpath(fullfile(githubDir, 'spikes')))
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\cbrewer2'))
%% specify folder that has *.bin raw data file and probe type
labels{1} = 'NP20';
labels{2} = 'NP10';
labels{3} = 'NPUHD2';
%% pre-processing to get mean waveform from raw data *.bin files
% ksfolders{1} = '\\sahale.biostr.washington.edu\data\Subjects\ZYE_0041\2021-09-03\2\p1_g0\p1_g0_imec0';
% ksfolders{2} = '\\sahale.biostr.washington.edu\data\Subjects\LK_0003\2021-07-06\1\p0_g0';
% ksfolders{3} = '\\sahale.biostr.washington.edu\data\Subjects\JRS_0033\2024-02-27\4\p3_g0\p3_g0_imec1';
% for i = 1:3
%     ksSubFolder = ksSubFolder{i};
%     waveformMean = getMeanWaveform(ksSubFolder);
%     save([label '_meanWaveform.mat'],'waveformMean','unitIDs');
% end
%% load meanWavform matrix 
% waveformMean: ncluster * nchan (384) *tSampleN (82)
id = 2;
label = labels{id};
load(['meanWaveform_' label '.mat']);
%% load channel map
chanMap_folder = 'C:\Users\Steinmetz lab\Documents\git\Kilosort2_J\configFiles';
[xcoords,ycoords] = loadChanMap(chanMap_folder,label);
%% subtract baseline for mean waveform
waveformMean = waveformMean*2.34;                                          % multiply by gain from .bin data
wfBaseline = mean(waveformMean(:,:,1:10),3);                                        % subtract baseline for each chan
wfBaseline2 = repmat(wfBaseline,[1,1,82]);
waveformMean = waveformMean-wfBaseline2;
%% get footprints
footprint = nan(size(waveformMean,1),1);
shank_spacing = 250;                                                       % specify shank spacing (use 250, if single shank)
for i = 1:size(waveformMean,1)
    thisWF = squeeze(waveformMean(i,:,:));
    footprint(i,1) = getFootprint(thisWF,xcoords,ycoords,shank_spacing);
end
%% plot waveforms
% yscale = 0.2; 
% siteSz = 12;  
% siteN = 36;
yscale = 0.1;                                                              % scale peak voltage, to avoid overlapping waveforms
siteSz = 6;                                                                % scale waveform length
siteN = 48;                                                                % total sites to plot (max 96 for 4 shank NP20)
indx = [1:5];
h1 = plotWaveform(waveformMean(indx,:,:),footprint(indx),xcoords,ycoords,...
    siteN,siteSz,yscale,shank_spacing);                                   % plot example waveform
% save_folder = fullfile('C:\Users\Steinmetz lab\Documents\MATLAB\footprint\main2');
% print(h1,fullfile(save_folder,['footprint_' label '_2.pdf']),...
%     '-dpdf', '-bestfit', '-painters');