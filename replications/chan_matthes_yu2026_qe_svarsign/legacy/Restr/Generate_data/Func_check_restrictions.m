function condition = Func_check_restrictions(R)

[n, m] = size(R);

    c = [];
    for i = 1:m-1
        i
        for j = i+1:m
            j
            Ri = R(:,i); Rj = R(:,j); Rij = Ri.*Rj;
            if length(find(Rij==1)) ~= 0 && length(find(Rij==-1)) ~= 0
                c = [c, 1];
            else
                c = [c, 0];
            end
            c
        end
    end
    condition = (length(find(c==0)) == 0);

end