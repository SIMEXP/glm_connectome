perc_disc = glmc_build_perc_disc(vol_perc,mask)
% SYNTAX: perc_disc = glmc_build_perc_disc(vol_perc,mask)
%
% VOL_PERC (3D volume) the map of percentage of discovery
% MASK (3D volume) the networks that were used to generate the map of percentage
%
% PERC_DISC (scalar) the overall percentage of discovery (average of percentage 
%   of discovery across all networks).

list_net = unique(mask);
list_net = list_net (list_net~=0);

perc_disc = 0;
for nn = 1:length(list_net)
    val = unique(vol_perc(mask==list_net(nn));
    if length(val)>1
        error('More than one percentage of discovery was found in network %i. This should not be the case',list_net(nn));
    end
    perc_disc = perc_disc + val(1);
end
perc_disc = perc_disc / length(list_net);