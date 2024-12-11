function waveformMean = getMeanWaveform(ksSubFolder)
% extract mean waveform from *.bin file

% ksfolders{1} = '\\sahale.biostr.washington.edu\data\Subjects\ZYE_0041\2021-09-03\2\p1_g0\p1_g0_imec0';
% ksfolders{2} = '\\sahale.biostr.washington.edu\data\Subjects\LK_0003\2021-07-06\1\p0_g0';
% ksfolders{3} = '\\sahale.biostr.washington.edu\data\Subjects\JRS_0033\2024-02-27\4\p3_g0\p3_g0_imec1';
nspikeThresh = 100;
ksSubFolder  = ksfolder;
[wf1,spikeTimesAll,spikeCluAll] = filterCluster(ksSubFolder,nspikeThresh);
unitIDs = wf1.unitIDs;
waveformMean = wf1.waveFormsMean;
% save(fullfile(save_folder, [label '_meanWaveform.mat']),'waveformMean','unitIDs');