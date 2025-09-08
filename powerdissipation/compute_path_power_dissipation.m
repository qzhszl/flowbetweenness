function [power_givenst,EdgeEnergy_SP,connected_flag] = compute_path_power_dissipation(G,node_source, node_destination)
    % force the current flow following a path, then show what the power
    % dissipation.
    % node s and t must connected with each other
    m = numedges(G);
    EdgeEnergy_SP = zeros(m,1);
    [shortest_path, shortest_pathlength] = shortestpath(G,node_source,node_destination);
    if isfinite(shortest_pathlength)
        connected_flag = 1;
        power_givenst = shortest_pathlength;
        for idx = 1:(numel(shortest_path)-1)
            u = shortest_path(idx);
            v = shortest_path(idx+1);
            eID = findedge(G,u,v);
            r = G.Edges.Weight(eID);
            EdgeEnergy_SP(eID) = r;
        end
    else
        power_givenst = 0;
        connected_flag = 0;
    end
end