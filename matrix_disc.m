function [ZC]=matrix_disc(x,scale,rc)

if nargin<3
    rc=1;
end    
[n,m]=size(x);
XN=datascaling(x,rc);
ZC=round(XN*(scale-1))+1; 
