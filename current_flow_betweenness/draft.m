clear,clc
% unweighted path graph 4 nodes
A = GenerateERfast(100,0.2,1);
G = graph(A);
plot(G,'EdgeLabel',G.Edges.Weight,'NodeColor',[0.8500 0.3250 0.0980], ...
'EdgeAlpha',0.5,'LineWidth',1,'MarkerSize',7,'EdgeLabelColor',[0 0.4470 0.7410],'NodeFontSize',10);


[nodeCB, edgeCB, edges] = current_flow_betweenness(A, 'normalize', false);
disp('edges (i j w) and edge normalized CF betweenness:');
disp([edges, edgeCB]);
disp('node normalized CF betweenness:'); disp(nodeCB);



[nodeflowbetweeness] = FlowNodeBetweeness(G)/4;
diff = find(abs(nodeflowbetweeness-nodeCB)>0.0000001)

[linkendpoint, flowlinkbetweeness] = FlowLinkBetweeness(G);

linkbetweenness2 = [linkendpoint,flowlinkbetweeness];

edgeEndpoints = sort(edges(:,1:2), 2);      % 无向边时排序
linkEndpoints = sort(linkendpoint, 2);

% 找出匹配关系
[tf, loc] = ismember(edgeEndpoints, linkEndpoints, 'rows');

% 匹配的索引
matchedEdges = find(tf);
matchedLinks = loc(tf);

% 拿出对应的值
edgeValues = edgeCB(matchedEdges);
linkValues = flowlinkbetweeness(matchedLinks);

% 用 find 找不一样的边
diffIdx = find(abs(edgeValues- linkValues)>0.000001);

% 输出不一样的边信息
badEdges   = matchedEdges(diffIdx);   % edges 中的行号
badLinks   = matchedLinks(diffIdx);   % linkendpoint 中的行号
badEdgeVal = edgeValues(diffIdx);
badLinkVal = linkValues(diffIdx);

% 打包成表格
T = table(badEdges, badLinks, badEdgeVal, badLinkVal)




function [nodeCB, edgeCB, edges] = current_flow_betweenness(A, varargin)
% CURRENT_FLOW_BETWEENNESS  Compute current-flow betweenness for nodes and edges
%
%   [nodeCB, edgeCB, edges] = current_flow_betweenness(A)
%   [nodeCB, edgeCB, edges] = current_flow_betweenness(A,'normalize',true)
%
% Inputs:
%   A         - n-by-n symmetric adjacency matrix (weighted). A(i,j)=0 means no edge.
%               Should be nonnegative and symmetric for undirected graph.
% Optional name-value:
%   'normalize' - logical (default false). If true:
%                   - edgeCB is divided by number of unordered pairs n*(n-1)/2
%                   - nodeCB is divided by number of unordered pairs not involving the node: (n-1)*(n-2)/2
%
% Outputs:
%   nodeCB    - n-by-1 vector of node current-flow betweenness (unnormalized or normalized)
%   edgeCB    - m-by-1 vector of edge current-flow betweenness (corresponding to edges list)
%   edges     - m-by-3 array: [i, j, w] listing edges (i<j)
%
% Notes:
%   - Complexity: roughly O(m * n log n) due to per-edge sorting; building L^+ costs O(n^3).
%   - This implementation follows the electrical-flow (current-flow) definition:
%       edge contribution for pair (s,t) is |I_ij^{st}| = w_ij * | x_s - x_t |,
%     where x = L^+_{i,:} - L^+_{j,:}.
%
% Example:
%   % small test (unweighted path of 4 nodes)
%   A = diag(ones(3,1),1) + diag(ones(3,1),-1);
%   [nc, ec, E] = current_flow_betweenness(A, 'normalize', true)

% Parse input
p = inputParser;
addRequired(p,'A',@(x) isnumeric(x) && ismatrix(x) );
addParameter(p,'normalize',false,@(x) islogical(x) || x==0 || x==1);
parse(p,A,varargin{:});
normalize = logical(p.Results.normalize);

% Validate and preprocess
A = double(A);
if any(A(:)<0)
    error('Adjacency matrix must have nonnegative weights.');
end
if ~isequal(A,A.')
    warning('Adjacency matrix not symmetric: symmetrizing by (A + A.'')/2.');
    A = (A + A.')/2;
end

n = size(A,1);
if n ~= size(A,2)
    error('Adjacency matrix must be square.');
end

% Build Laplacian L = D - A
deg = sum(A,2);
L = diag(deg) - A;

% Compute Moore-Penrose pseudoinverse Lplus of L via spectral decomposition
% (more stable than pinv for Laplacian singularity)
opts.disp = 0;
[V, D] = eig(L);
d = diag(D);
% numerical threshold
tol = max(n * eps(max(d)), 1e-12);
nonzero = d > tol;
if sum(nonzero) == 0
    error('Graph Laplacian is zero (no edges).');
end
% Construct pseudoinverse: invert nonzero eigenvalues
Dplus = zeros(n);
Dplus(nonzero,nonzero) = diag(1./d(nonzero));
Lplus = V * Dplus * V';

% Build edge list (upper triangle)
[i_idx, j_idx, w_vals] = find(triu(A,1));
m = numel(w_vals);
edges = [i_idx, j_idx, w_vals];

edgeCB = zeros(m,1);
% For computing node CB we will accumulate two things:
%   sumIncidentEdgeCB(v) = sum of edgeCB for edges incident to v
%   and sumIncidentEndpointPairs(v) = sum over incident edges of sum_{t != v} w * | x_v - x_t |
sumIncidentEdgeCB = zeros(n,1);
sumIncidentEndpointPairs = zeros(n,1);

% Pre-allocate temp vectors
x = zeros(n,1);

for e = 1:m
    i = i_idx(e);
    j = j_idx(e);
    w = w_vals(e);
    % x = Lplus(i,:) - Lplus(j,:)  (row vector)
    x = (Lplus(i,:) - Lplus(j,:)).';    % column vector length n

    % compute sum_{s<t} |x_s - x_t| efficiently via sorting
    xs = sort(x); % ascending
    k = (1:n).';
    % sum_{s<t} (x_t - x_s) where sorted => sum_{s<t}|...|
    sdiff = sum( xs .* (2*k - n - 1) );  % scalar

    edge_contrib = w * sdiff;
    edgeCB(e) = edge_contrib;

    % accumulate for endpoints
    sumIncidentEdgeCB(i) = sumIncidentEdgeCB(i) + edge_contrib;
    sumIncidentEdgeCB(j) = sumIncidentEdgeCB(j) + edge_contrib;

    % Now compute contributions for pairs where one endpoint is endpoint node,
    % i.e., sum_{t != i} w * | x_i - x_t |  and sum_{t != j} w * | x_j - x_t |
    % We'll compute both by direct vectorized absolute diffs (O(n)).
    xi = x(i);
    xj = x(j);
    % sum over t (including t=i gives zero) -> same as excluding i
    sum_abs_i = w * sum( abs(xi - x) );
    sum_abs_j = w * sum( abs(xj - x) );

    sumIncidentEndpointPairs(i) = sumIncidentEndpointPairs(i) + sum_abs_i;
    sumIncidentEndpointPairs(j) = sumIncidentEndpointPairs(j) + sum_abs_j;
end

% Node betweenness:
% For node v, nodeCB_v = 0.5 * sum_{e incident v} (edgeCB_e - sum_{t != v} w_e * |x_v - x_t|)
% We have sumIncidentEdgeCB(v) = sum_{e inc v} edgeCB_e
% and sumIncidentEndpointPairs(v) = sum_{e inc v} sum_{t != v} w_e * |x_v - x_t|
nodeCB = 0.5 * ( sumIncidentEdgeCB - sumIncidentEndpointPairs );

% Normalization (optional)
if normalize
    % number of unordered pairs (s,t): Npairs = n*(n-1)/2
    Npairs = n*(n-1)/2;
    edgeCB = edgeCB / Npairs;
    % For node normalization: number of unordered pairs s<t with s,t != v is (n-1)*(n-2)/2
    denom_nodes = (n-1)*(n-2)/2;
    if denom_nodes > 0
        nodeCB = nodeCB / denom_nodes;
    else
        nodeCB = nodeCB; % small n case, leave as is
    end
end

end
