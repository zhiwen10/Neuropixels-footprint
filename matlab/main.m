%% Add paths 
githubDir = 'C:\Users\Steinmetz lab\Documents\git';
addpath(genpath(fullfile(githubDir, 'Neuropixels-footprint')))
%% specify probe type
labels{1} = 'NP20';
labels{2} = 'NP10';
labels{3} = 'NPUHD2';
labels{4} = 'NPUHD';
id = 2;
label = labels{id};
%% save mean waveform from *.bin files and kilosort results (optional)
% data_preprocessing;

%% load meanWavform matrix: ncluster * nchan (384) *tSampleN (82)
data_folder = 'data';
load(fullfile(data_folder,['meanWaveform_' label '.mat']));

%% load channel map
chanMap_folder = 'chanMaps';
[xcoords,ycoords] = loadChanMap(chanMap_folder,label);

%% get footprints
footprint = nan(size(waveformMean,1),1);
shank_spacing = 250;                                                       % specify shank spacing (use 250, if single shank)
for i = 1:size(waveformMean,1)
    thisWF = squeeze(waveformMean(i,:,:));
    footprint(i,1) = getFootprint(thisWF,xcoords,ycoords,shank_spacing);
end
% ft = table(footprint);
% writetable(ft,[label '_footprint_matlab.csv']);

%% plot waveforms
yscale = 0.1;                                                              % scale peak voltage, to avoid overlapping waveforms
siteSz = 6;                                                                % scale waveform length
siteN = 48;                                                                % total sites to plot (max 96 for 4 shank NP20)
indx = [1:5];
h1 = plotWaveform2(waveformMean(indx,:,:),footprint(indx),...
    xcoords,ycoords,siteN,siteSz,yscale,shank_spacing);                    % plot example waveform

% print(h1,['footprint_' label '_2.pdf'],...
%     '-dpdf', '-bestfit', '-painters');