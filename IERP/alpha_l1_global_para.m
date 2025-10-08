function alpha = alpha_l1_global_para(A, D)
    r = D(A>0)./A(A>0);
    w = A(A>0)./D(A>0);
    
    % Sort and find weighted median
    [rsort, idx] = sort(r(:));
    w = w(idx);
    cw = cumsum(w)/sum(w);
    alpha = rsort(find(cw>=0.5, 1));
end
