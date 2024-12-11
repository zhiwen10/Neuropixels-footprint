function h1 = plotWaveform2(wfs,radius,xcoords,ycoords,siteN,siteSz,yscale,shank_spacing)
%% load data
% wfs: cluster*nChan*tSampleN
tRange_old = [-40 40];
tRange_new = 25:80;
wfs = wfs(:,:,tRange_new);
%%
xc =sort(unique(xcoords)); yc = sort(unique(ycoords));
tRangeN = size(wfs,3);           
nchan =  size(wfs,2);
%%
h1 = figure('Renderer', 'painters', 'Position', [50 50 900 500]);
ncluster = size(wfs,1);
for kk = 1:ncluster
    thisWF = squeeze(wfs(kk,:,:));
    %% get amp map
    troughs_min = min(thisWF(:));
    peaks_max = max(thisWF(:));
    sign = 1;                                                                  % defalut is absolute trough is larger than peak
    if abs(troughs_min) < abs(peaks_max)
        sign = 0;
        thisWF = -thisWF;
    end
    amps = zeros(1,nchan);
    for chan = 1:nchan
        clear trough_h peak_h
        thisChanWF = squeeze(thisWF(chan,:));
        [~,trough_idx] = min(thisChanWF);
        thisChanWF2 = thisChanWF(trough_idx(1):end);
        peak_idx = find(thisChanWF2 == max(thisChanWF2),1)+trough_idx-1;

        trough_h = thisChanWF(trough_idx);
        peak_h = thisChanWF(peak_idx);
        amps(chan) = peak_h+abs(trough_h);
    end
    %% select 32 sites closest to the peak amp site
    amps1 = sort(amps);
    amps = amps-mean(amps1(1:10));
    [~,max_chan] = max(amps);
    max_pos = [xcoords(max_chan), ycoords(max_chan)];
    shank_pos = round(xcoords/shank_spacing)+1;
    max_shank_pos = round(max_pos(1)/shank_spacing)+1;
    allSitesOnShank_index = (shank_pos==max_shank_pos);
    allSitesOnShank_I = find(allSitesOnShank_index);
    if iscolumn(xcoords)
        xcoords = xcoords';
        ycoords = ycoords';
    end
    allSitesPos = [xcoords;ycoords]';
    allSitesPos = allSitesPos(allSitesOnShank_index,:);
    point_diff = allSitesPos-max_pos;
    distance = vecnorm(point_diff,2,2);
    [distance1,I] = sort(distance);
    I2 = allSitesOnShank_I(I);
    sites = I2(1:siteN);
    %% plot current cluster
    ax(kk) = subplot(1,ncluster,kk);
    xMin(kk) = min(xcoords(sites));
    xMax(kk) = max(xcoords(sites));
    for i = 1:numel(sites)
        site = sites(i);
        wv_chan = squeeze(thisWF(site,:));
        trH(i) = plot(xcoords(site)+linspace(1, siteSz, tRangeN), ...
            ycoords(site)+ wv_chan*yscale, 'color',[0.8,0.8,0.8], 'LineWidth', 1);
        hold on;
    end
    text(max_pos(1),max_pos(2),'*','fontSize',12);
    %% hihlight channels within radius
    threshold = find(distance-radius(kk)<=0);
    I3 = allSitesOnShank_I(threshold);
    xcoords3 = xcoords(I3);
    ycoords3 = ycoords(I3);
    distance3 = distance(threshold);
    thisWF3 = thisWF(I3,:);
    [~,I4] = sort(distance3);
    thisWF3_sort = thisWF3(I4,:);
    xcoords3_sort = xcoords3(I4);
    ycoords3_sort = ycoords3(I4);
    nsites = size(thisWF3_sort,1);
    % color1 = flipud(hot(nsites));
    color1 = flipud(cbrewer2('seq','reds',nsites));
    for i = 1:nsites
        hold on;
        wv_chan = squeeze(thisWF3_sort(i,:));
        trH(i) = plot(xcoords3_sort(i)+linspace(1, siteSz, tRangeN), ...
            ycoords3_sort(i)+ wv_chan*yscale, 'color',color1(i,:), 'LineWidth', 1);
    end
    title(['radius=' num2str(radius(kk)) ' \mum']);
    yscaleBar = 30;
    plot([xMin(kk),xMin(kk)],[max_pos(2),max_pos(2)+yscaleBar*yscale],'r','lineWidth',1);
end
%
xLength = max(xMax-xMin);
for kk = 1:ncluster
    xlim(ax(kk),[xMin(kk),xMin(kk)+xLength+siteSz]);
end
