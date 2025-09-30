clear,clc

for N = 1000
    avg = [8,10,20,50,100,200,500];
    p_vec = avg./(N-1)
    for p=p_vec
        simu_flowbetweenness_distribution(N,p,1)
    end
end
