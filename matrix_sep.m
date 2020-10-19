function [xs,xs2]=matrix_sep(x,ZC,scale)

middle_value=median(1:scale);

xs{1}=NaN(size(x));
xs{middle_value}=NaN(size(x));
xs{scale}=NaN(size(x));

idx1=ZC==1;
xs{1}(idx1)=x(idx1);

idx2=ZC==middle_value;
xs{middle_value}(idx2)=x(idx2);

idx3=ZC==scale;
xs{scale}(idx3)=x(idx3);

idxa=logical(mod(idx1+idx2+idx3,2));

xs2=NaN(size(x));
xs2(idxa)=1;