function row_spacing = getRowSpacing(xcoords,ycoords)
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