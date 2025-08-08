function [y,l_sg] = flowsubgraph(G,node_source,node_sink)
    % l_sg: the length(number of links) of the flow subgraph
    % y: links on the flow subgraph
    if isinf(distances(G, node_source,node_sink))
        l_sg = 0;
        y=[];
    else
    %     Linknumber = numedges(G);
        N = numnodes(G);              % node number
        Linkweight = G.Edges.Weight;
        Lmatrix = diag(Linkweight);  % when link weight w_l = 1/r
        % Lmatrix = inv(diag(Linkweight));  % when link weight w_l = r
        B = -full(incidence(G));
        Q=B*Lmatrix*B.';
        PQ = pinv(Q);
        C = Lmatrix*B.'*PQ.';
        e_i = zeros(N,1);
        e_i(node_source)=1;
        e_j = zeros(N,1);
        e_j(node_sink)=1;
        y = C*(e_i-e_j);   % the current on eack link
        %     u =PQ*e_i-PQ*e_j;   % the voltage of each node
        l_sg = length(find((abs(y))>0.0000001));
    end
end

