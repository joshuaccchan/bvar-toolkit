function  d_Q=d_normal(Sig,Sigd,Q,mu,mud,forecast_period, Sample)

    S=kron(Sig,ones(forecast_period,1));
    S_d=kron(Sigd,ones(forecast_period,1));
    
    Q=repmat(Q,Sample,1);
    f=normpdf(Q,mu, S);
    minusdF=repmat(f.*(Q-mu)./(2*S),1,3).*S_d+repmat(f,1,3).*mud;
  
    minusdF=squeeze(mean(reshape(minusdF,forecast_period,Sample,3),2));
    f=repmat(mean(reshape(f,forecast_period,Sample),2),1,3);
    d_Q=minusdF./f;
