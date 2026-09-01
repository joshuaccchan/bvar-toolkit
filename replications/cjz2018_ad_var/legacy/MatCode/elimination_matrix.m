function [m0] = elimination_matrix(n)
    m0 = zeros(0.5*n*(n+1), n^2);
    src = 1;
    tgt = 1;
    for j = 1:n
        for i = 1:n           
            if i >= j
                m0(tgt, src) = 1;
                tgt = tgt + 1;
            end
            src = src + 1;
        end
    end
