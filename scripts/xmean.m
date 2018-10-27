function x1 = xmean(xsim)
[n, ~] = size(xsim);
xm = mean(xsim(n - 300:n));
x1 = zeros(n,1);
x1 = x1 + xm;
end
