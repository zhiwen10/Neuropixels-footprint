function [wf1,spikeTimesAll,spikeCluAll] = filterCluster(ksSubFolder,nspikeThresh)
sp = loadKSdir(ksSubFolder);
%% keep cluster if they are good by KS default or mabually curated in phy2
% if manually curated in phy2
% gcluster = sp.cids(sp.cgs==2);
%if not manually curated in phy, and use rez.good automated selected from KS2
load([ksSubFolder '\rez.mat']);
gcluster = find(rez.good==1)-1;
% gcluster = 0:40;
gcluster = gcluster';
if ~isempty(gcluster)
    %% filter by spike numbers (by default >50)
    for m = 1:length(gcluster)
        UnitID = gcluster(m);
        SpikeTimes = sp.st(sp.clu==UnitID);
        nSpikes(m) = size(SpikeTimes,1);
    end
    gcluster = gcluster(find(nSpikes>nspikeThresh));
    %%
    [lic,loc] = ismember(sp.clu,gcluster);
    spikeTimesAll = sp.st(lic);
    spikeCluAll = sp.clu(lic);
    %%
    gwfparams.dataDir = ksSubFolder;
    dir1 = dir(fullfile(ksSubFolder, '*.ap.bin'));
    % gwfparams.fileName = fullfile(dir1.folder, dir1.name);
    gwfparams.fileName =  dir1.name;
    gwfparams.dataType = sp.dtype;
    gwfparams.nCh = sp.n_channels_dat;
    gwfparams.wfWin = [-40 41];
    gwfparams.nWf = 100; % randomly choose 100 spikes from each cluster
    % gwfparams.spikeTimes = sp.st;
    % gwfparams.spikeClusters = sp.clu;
    gwfparams.spikeTimes = spikeTimesAll;
    gwfparams.spikeClusters = spikeCluAll;
    gwfparams.Fs = sp.sample_rate;
    % In wf, wf.waveformsMean: nClu x 256 (sites) x 82 (nSamp)
    wf1 = getWaveForms2(gwfparams);
else
    wf1 = []; spikeTimesAll = []; spikeCluAll = [];
end
end
