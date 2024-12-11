%% preprocess data 
% save mean waveform from *.bin files, by averaging waveforms from raw data
% using kisosort timestamps 
ksfolders{1} = '\\sahale.biostr.washington.edu\data\Subjects\ZYE_0041\2021-09-03\2\p1_g0\p1_g0_imec0';
ksfolders{2} = '\\sahale.biostr.washington.edu\data\Subjects\LK_0003\2021-07-06\1\p0_g0';
ksfolders{3} = '\\sahale.biostr.washington.edu\data\Subjects\JRS_0033\2024-02-27\4\p3_g0\p3_g0_imec1';
labels{1} = 'NP20';
labels{2} = 'NP10';
labels{3} = 'NPUHD2';
%%
for i = 1:3
    ksSubFolder = ksfolders{i};
    label = labels{i};
    waveformMean = data_preprocessing(ksSubFolder,label);
    %%
    waveformMean = waveformMean*2.34;                                          % multiply by gain from .bin data
    % subtract baseline for mean waveform
    wfBaseline = mean(waveformMean(:,:,1:10),3);                               % subtract baseline for each chan
    wfBaseline2 = repmat(wfBaseline,[1,1,82]);
    waveformMean = waveformMean-wfBaseline2;
    save(['meanWaveform_' label '_2.mat'],'waveformMean');
    writeNPY(waveformMean,['meanWaveform_' label '.npy']);
end
