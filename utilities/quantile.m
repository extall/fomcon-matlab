function y = quantile(x)
%QUANTILE Compute quantile function Phi^-1(x)
    y = sqrt(2)*erfinv(2*x-1);
end

