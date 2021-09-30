function [testfunction]=harmfit(coefin,t)

% [output]=harmfit(coefin, t)
% coefin are the estimated amplitudes and phases of the constituents
% coefin(1)=mean
% coefin(2:Nharm+1) is Cn, Nharm the number of harmonics 
% coefin(Nharm+2:2*Nharm+1) is Dn
% t is the time vector
% output = mean + sum_k=1^k=Nharm Cn*sin(wn(k)*t) + Dn*cos(wn(k)*t)

global wn

a0=coefin(1);
Nharm=(length(coefin)-1)/2;
Cn=coefin(2:Nharm+1);
Dn=coefin(Nharm+2:length(coefin));
testfunction(1:length(t))=a0;

for k=1:length(Cn)
    testfunction=testfunction+Cn(k)*sin(wn(k)*t) + Dn(k)*cos(wn(k)*t);
end

end