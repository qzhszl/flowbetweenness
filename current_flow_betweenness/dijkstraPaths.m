function [dist, sigma, P, S] = dijkstraPaths(adj, s)
n = size(adj,1);
dist = inf(n,1);
sigma = zeros(n,1);
P = cell(n,1);
S = [];

dist(s) = 0;
sigma(s) = 1;

visited = false(n,1);
Q = containers.Map('KeyType','int32','ValueType','double');
Q(s) = 0;

while ~isempty(Q)
    % extract min dist
    [~,idx] = min(cell2mat(values(Q)));
    keysQ = cell2mat(keys(Q));
    v = keysQ(idx);
    dv = Q(v);
    remove(Q,v);
    
    if visited(v), continue; end
    visited(v) = true;
    
    S(end+1) = v;
    for w = find(adj(v,:)>0)
        vw_dist = dv + adj(v,w);
        if vw_dist < dist(w)
            dist(w) = vw_dist;
            sigma(w) = sigma(v);
            P{w} = v;
            Q(w) = dist(w);
        elseif abs(vw_dist - dist(w))<1e-12
            sigma(w) = sigma(w) + sigma(v);
            P{w}(end+1) = v;
        end
    end
end
end


