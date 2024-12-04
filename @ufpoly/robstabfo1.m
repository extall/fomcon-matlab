function is_robuststable = robstabfo1(P)
%ROBSTABFO1 - Robust stability analysis of fractional-order polynomials.
%
% Syntax:
%   is_robuststable = robstabfo1(P)
%
% Inputs:
%   P - A fractional-order polynomial with uncertainty intevals
%
% Outputs:
%   is_robuststable - Logical result indicating robust stability.
%                     A message is also displayed indicating whether the system is robustly stable.
%
% Description:
%   This function evaluates the robust stability of fractional-order polynomials by generating
%   all combinations of coefficients within their bounds and analyzing the system's stability
%   over a defined frequency range. It determines whether or not the origin is included in 
% the value set of the inteval polynomial.
% Example:
%   P = ufpoly('[1,2]s^1.8 + s + s^.65 + 10');
%   robstabfo1(P)
lowerbounds = fliplr(P.a(:,1)');
upperbounds = fliplr(P.a(:,2)');
alpha =fliplr(P.na(:,1)');
had=vertcat(lowerbounds,upperbounds);
n=numel(alpha);
r=2^n;
c=n;
h=nan(r,c);
inx=[];


for i1=1:c
    st=2^(i1);
    id=ones(st/2,r/(st/2));
    id(:,2:2:end)=2;
    idx=reshape(id,[r,1]);
    inx=horzcat(inx,idx);
end

h=[];

for j1=1:r
    for k=1:c
        h(j1,k)=had(inx(j1,k),k);
    end
end
hh= sum(max(abs(lowerbounds(1:end-1)),abs(upperbounds(1:end-1))))/(min(abs(lowerbounds(end)),abs(upperbounds(end))));
hhh= min(hh^(1/(alpha(end)-alpha(end-1))),500);
w= 0:.01:hhh;
s= j*w;
laplas=[];
for a=1:numel(alpha)
    laplas(a,:)=s.^alpha(a);
end

results=[];
sumcoef=[];
for b=1:r
    for d=1:numel(alpha)
        sumcoef(d,1:numel(w))=h(b,d).*laplas(d,:);
    end
    results(b,1:numel(w))=sum(sumcoef,1);
end
%TrEnq = []
for d0 = 1:r
    for d1 = 1:r

        TrEn(d1,:) =  abs(results(d0,:)) + abs(results(d1,:)) - abs(results(d0,:) - results(d1,:));

    end
    TrEnq(d0,:) =  min(TrEn(d1,:));
end

RC = min(TrEnq,[],'all');


if RC<.00001
    display("the system is not robustly stable.")
else
    display ("the system is robustly stable.")

end