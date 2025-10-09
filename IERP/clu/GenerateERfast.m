function A = GenerateERfast(n,p,weighted)
    A = rand(n,n) < p;
    A = triu(A,1);
    if weighted == 0
     
    elseif weighted == 1
        linkweight_matrix = rand(n,n);
        A = A.*linkweight_matrix;
    else
        linkweight_matrix = randi(weighted,n,n);
        A = A.*linkweight_matrix;
    end
    
    A = A + A';
end
