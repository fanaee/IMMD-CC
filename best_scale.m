function ns=best_scale(x,m,min_scale,max_scale)

k=0;
for s=min_scale:2:max_scale
    k=k+1;
    ZXX{k}=matrix_disc(x,s,1);
    mv=median(1:s);
    ns(k,1)=length(find(ZXX{k}==1)); ns(k,2)=length(find(ZXX{k}==mv)); ns(k,3)=length(find(ZXX{k}==s));
    ntt=sum(ns(k,1:3));
    ns(k,4)=abs(ns(k,2)-ns(k,1))/ntt+abs(ns(k,2)-ns(k,3))/ntt+abs(ns(k,1)-ns(k,3))/ntt;
    ns(k,5)=s;
end
