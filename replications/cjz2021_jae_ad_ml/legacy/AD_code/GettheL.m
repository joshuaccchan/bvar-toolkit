    function L=GettheL(S,n)
        L.v=chol(S.inv)';
        L.d=d_cholasky(L.v,S.invd);
        L.dt=d_trans(L.v,L.d);
   L.inv=(L.v\speye(n));
   L.tinv=(L.v'\speye(n));
    end 