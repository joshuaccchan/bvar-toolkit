function S=dmatrix(S,n,k)
  S.inv=S.v\speye(n);
  S.invd=d_minverse(S.inv,S.d);
  S.det=det(S.v);
  S.ddet=d_det(S.det,S.v,S.d,n,k);