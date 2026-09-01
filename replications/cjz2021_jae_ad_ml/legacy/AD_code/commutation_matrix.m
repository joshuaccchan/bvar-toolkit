


    function m0=commutation_matrix(r,c)
    i=1:r;
    j=1:c;
    src=reshape(r*repmat((j-1),r,1)+repmat(i',1,c),r*c,1);
    tgt=reshape((c*repmat((i-1),c,1)+repmat(j',1,r))',r*c,1);
    m0=sparse(tgt,src, ones(r*c,1));