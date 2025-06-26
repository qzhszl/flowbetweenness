clear,clc
% read data from clu, sum all different external simulation and save it in
% a new file

N_vec = [200];
p_start_vec = zeros(length(N_vec),1);
countN = 1; 
for N = N_vec
    p_start_vec(countN) = round(log(N)/N,4);
    countN = countN+1;
end

countN = 1;
external_simutimes = 20;
for N = N_vec
    p_vec = linspace(p_start_vec(countN), 1, 15);
    p_vec = round(p_vec,4);
    countp = 1;
    for p= p_vec
        filename = sprintf("D:\\data\\flow betweenness\\IERP\\IERP_N%dERp%.4f.txt",N,p);
        final_result = [];
        for ex_simu_time = 1:external_simutimes
            filename_each_simu = sprintf("D:\\data\\flow betweenness\\IERP\\IERP_N%dERp%.4fSimu%d.txt",N,p,ex_simu_time);
            results = readmatrix(filename_each_simu);
            final_result = [final_result;results];
            delete(filename_each_simu);
        end
        writematrix(final_result,filename)
    end
    countN = countN+1;
end

