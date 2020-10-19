function [xn,mx,mn]=datascaling(x,m)

if m<3
    mx=max(x,[],m);
    mn=min(x,[],m);
else
    mx=max(x(:));
    mn=min(x(:));
end

xn = bsxfun(@minus, x, mn);
xn = bsxfun(@rdivide, xn, mx-mn);
