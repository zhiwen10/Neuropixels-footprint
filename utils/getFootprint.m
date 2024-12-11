function footprint = getFootprint(thisWF,xcoords,ycoords, shank_spacing)
%% get amp map by suming trough and peak after trough 
troughs_min = min(thisWF(:));
peaks_max = max(thisWF(:));
% if abs(peak)>abs(trough), i.e., inverted waveform, then flip the waveform
% sign. Default is absolute trough larger than peak
sign = 1;                                                                  
if abs(troughs_min) < abs(peaks_max)
    sign = 0;
    thisWF = -thisWF;
end
amps = zeros(1,384);
for chan = 1:384
    clear trough_h peak_h
    thisChanWF = squeeze(thisWF(chan,:));
    [~,trough_idx] = min(thisChanWF);
    thisChanWF2 = thisChanWF(trough_idx(1):end);
    peak_idx = find(thisChanWF2 == max(thisChanWF2),1)+trough_idx-1;

    trough_h = thisChanWF(trough_idx);
    peak_h = thisChanWF(peak_idx);
    amps(chan) = peak_h+abs(trough_h);
end
amps1 = sort(amps);
amps = amps-mean(amps1(1:10));
%% if 4 shank probe, only use shank where peak site is at 
[amps_max,I_max] = max(amps);

maxAmpSite_shank = round(xcoords(I_max)./shank_spacing)+1;
alSites_shank = round(xcoords./shank_spacing)+1;
siteOnShank_index = (alSites_shank == maxAmpSite_shank);
xcoords = xcoords(siteOnShank_index);
ycoords = ycoords(siteOnShank_index);
amps = amps(siteOnShank_index);
thisWF = thisWF(siteOnShank_index,:);
%% find out row spacing
% if all sites are single column, with trunk staggering like double length
% NP2.0, then collapse into single column. Otherwise, fillmissing 2 mess up
xc1 = unique(xcoords);
ysort = sort(ycoords(xcoords==xc1(1)));
ysort_diff = diff(ysort);
[C,ia,ic] = unique(ysort_diff);
a_counts = accumarray(ic,1);
if isrow(C)
    C = C';
end
value_counts = [C, a_counts];
value_counts_sort = sortrows(value_counts,2,'descend');
row_spacing = value_counts_sort(1,1);
%% if too much staggering, then fillmissing2 fails
% therefore, we would collapse into a single column for NP20 double length
% situation
length_expected = round(numel(xcoords)./numel(unique(xcoords)))*row_spacing;
length_true = max(ycoords)-min(ycoords);
if length_true > length_expected*1.2 
    xcoords(:) = min(xcoords);
end
%% define a high resolution matrix using min and max of site position in chanMap
% at 1 um resolution, fill in the 384 sites within this matrix
% then intepolate the rest missing entries, using "fillmissing2"
xmin = min(xcoords); xmax = max(xcoords);
ymin = min(ycoords); ymax = max(ycoords);
width = ceil(xmax-xmin+1);
height = ceil(ymax-ymin+1);
amp2 = nan(height,width);
xcoords2 = round(xcoords-xmin+1);
ycoords2 = round(ycoords-ymin+1);
for i = 1:numel(xcoords)
    amp2(ycoords2(i),xcoords2(i)) = amps(i);
end
if numel(unique(xcoords))==1
    amp3 = fillmissing(amp2,"linear");
else
    amp3 = fillmissing2(amp2,"cubic");
end
    
[amp_max,I] = max(amp3(:));
[I1,I2] = ind2sub(size(amp3),I);
%% sample amp at different radius from the peak site
% first design the sampling position matrix
rho = 0:400;
theta = zeros(size(rho));
count1 = 1;
for beta = -pi:pi/6:pi
    [xq(count1,:),yq(count1,:)] = pol2cart(theta+beta,rho);
    count1 = count1+1;
end
xq2 = round(xq)+I2;
yq2 = round(yq)+I1;
%% read out the values at position matrix
vq = nan(size(xq2));
for i = 1:size(xq2,1)
    for j = 1:size(xq2,2)
        if yq2(i,j)>=1 & yq2(i,j)<=size(amp3,1) & xq2(i,j)>=1 & xq2(i,j)<=size(amp3,2)
            vq(i,j) = amp3(yq2(i,j),xq2(i,j));
        end
    end
end
%% average across samples at different radius
% find footprint as the radius at which voltage drops to 30uV
vq_mean = mean(vq,1,"omitnan");
footprint = find(vq_mean<=30,1);
if isempty(footprint)
    footprint = nan;
end