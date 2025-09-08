function edgeBC = edgeBetweennessCentrality(G)
% EDGEBETWEENNESSCENTRALITY Compute edge betweenness centrality
%
%   edgeBC = edgeBetweennessCentrality(G)
%
% Input:
%   G - MATLAB graph object (undirected, unweighted or weighted)
%
% Output:
%   edgeBC - m-by-1 vector of edge betweenness centrality values
%            (in the same order as G.Edges)
%
% Note:
%   - This implementation is based on Brandes (2001) algorithm.
%   - If G has edge weights, they are interpreted as distances
%     (shortest-path length). For unweighted graphs all weights = 1.
%   - Complexity: O(n*m) for unweighted, O(n*m + n^2 log n) for weighted.

n = numnodes(G);
m = numedges(G);
edgeBC = zeros(m,1);

if ~isempty(G.Edges.Weight)
    weights = G.Edges.Weight;
else
    weights = ones(m,1);
end

adj = adjacency(G, 'weighted');

% --- Brandes algorithm for edges ---
for s = 1:n
    % Initialization
    S = [];                     % stack
    P = cell(n,1);              % predecessors
    sigma = zeros(n,1);         % number of shortest paths
    dist = inf(n,1);            % distance
    
    sigma(s) = 1;
    dist(s) = 0;
    
    % Use Dijkstra if weighted, BFS if unweighted
    if all(weights==1)
        % --- BFS (unweighted) ---
        Q = java.util.LinkedList();
        Q.add(s);
        while ~Q.isEmpty()
            v = Q.remove();
            S(end+1) = v;
            neigh = find(adj(v,:));
            for w = neigh
                % Path discovery
                if dist(w) == inf
                    Q.add(w);
                    dist(w) = dist(v) + 1;
                end
                % Path counting
                if dist(w) == dist(v) + 1
                    sigma(w) = sigma(w) + sigma(v);
                    P{w}(end+1) = v;
                end
            end
        end
    else
        % --- Dijkstra (weighted) ---
        [dist, sigma, P, S] = dijkstraPaths(adj, s);
    end
    
    % Dependency accumulation
    delta = zeros(n,1);
    while ~isempty(S)
        w = S(end); S(end) = [];
        for v = P{w}
            c = (sigma(v)/sigma(w)) * (1 + delta(w));
            % Find edge (v,w)
            eid = findedge(G, v, w);
            edgeBC(eid) = edgeBC(eid) + c;
            delta(v) = delta(v) + c;
        end
    end
end
edgeBC = edgeBC / 2;
% % Undirected: divide by 2
% if ~isdirected(G)
%     edgeBC = edgeBC / 2;
% end

end

%% --- helper: Dijkstra shortest paths (for weighted graphs) ---

