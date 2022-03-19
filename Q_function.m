function [ Q ] = Q_function(x)
s=sqrt(2);
Q=0.5*erfc(x/s);
end

