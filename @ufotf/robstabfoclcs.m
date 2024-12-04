function is_robustlystable = robstabfoclcs(Controller, Plant,N)
% ROBSTABFOCLCS - Analyze the robust stability of an interval fractional-order closed-loop control system.
%
% Syntax:
%   is_robustlystable = robstabfoclcs(Controller, Plant, N)
%
% Inputs:
%   Controller - A structure representing the controller.
%   Plant      - A structure representing the plant.
%   N          - (Optional) Number of delay partitions for analysis. Default: N = 3.
%
% Outputs:
%   is_robustlystable - Logical flag indicating whether the system is robustly stable.
%
% Description:
%   This function evaluates the robust stability of a fractional-order system with
%   interval uncertainty in its parameters and delay. The robust stability is determined
%   based on bounds of coefficients and fractional orders for the plant and controller.
%   If the system fails the robust stability criterion, the function suggests increasing
%   the number of delay partitions (N) or declares the system unstable.
%
% Example 1:
%   % Define a controller and plant (with an interval time delay)
%Plant= ufotf('[1.3,1.7]s^.3 + [1.4,1.6]', '[15,25]s^1.6 + [2.5,3.5]s^.3+ [1.5,2.5]',[.01 .15])
 %Controller = ufotf('2s^.3 +.5','s^.3')
 %robstabfoclcs(Controller, Plant)
%
% Example 2:
%   % Define a controller and plant (without an interval time delay)
%Plant=ufotf('[1.3,1.7]s^.3 + [1.4,1.6]', '[15,25]s^1.6 + [2.5,3.5]s^.3+ [1.5,2.5]',0)
 %Controller = ufotf('2s^.3 +.5','s^.3')
 %robstabfoclcs(Controller, Plant)


if nargin < 3
    N=3;
end
lowerboundsN = fliplr(Plant.b.a(:,1)');
upperboundsN = fliplr(Plant.b.a(:,2)');
alphaN =fliplr(Plant.b.na(:,1)');
alphaNC =fliplr(Controller.b.na(:,1)');
lowerboundsD = fliplr(Plant.a.a(:,1)');
upperboundsD = fliplr(Plant.a.a(:,2)');
alphaD =fliplr(Plant.a.na(:,1)');
alphaDC =fliplr(Controller.a.na(:,1)');
CoNumCon = Controller.b.a(:,1)';
hCN= sum(abs(CoNumCon(1:end)));
CoDeCon = Controller.a.a(:,1)';
hCDd= sum(abs(CoDeCon(1:end)));
hCD= sum(abs(CoDeCon(1:end-1)));
if numel(alphaDC)<2
    alphax =0;
else
    alphax = alphaDC(end-1);
end
wm = (sum(max(abs(lowerboundsN(1:end)),abs(upperboundsN(1:end))))*hCN+ hCDd*sum(max(abs(lowerboundsD(1:end-1)),abs(upperboundsD(1:end-1)))) + hCD*max(abs(lowerboundsD(end)),abs(upperboundsD(end))))/(min(abs(lowerboundsD(end)),abs(upperboundsD(end)))*hCDd(1));
alph = alphaD(end)+ alphaDC(end) - max([alphaD(end)+ alphax alphaD(end-1)+ alphaDC(end) alphaN(end)+ alphaNC(end)]);
wmax =min(wm^(1/alph),500);
delay = Plant.ioDelay;
l1 = Plant.ioDelay(1,1);
l2 = Plant.ioDelay(1,2);
hhh=2*pi/(l2-l1);
if hhh>wmax
    hhh= wmax;
else
    hhh = 2*pi/(l2-l1);
end
had1=vertcat(lowerboundsN,upperboundsN);
n1=numel(alphaN);
r1=2^n1;
c1=n1;
h1=nan(r1,c1);
inx1=[];
for i1=1:c1
    st1=2^(i1);
    id1=ones(st1/2,r1/(st1/2));
    id1(:,2:2:end)=2;
    idx1=reshape(id1,[r1,1]);
    inx1=horzcat(inx1,idx1);
end

h1=[];

for j1=1:r1
    for k1=1:c1
        h1(j1,k1)=had1(inx1(j1,k1),k1);
    end
end
w1= 0:.01:hhh;
s1= j*w1;
laplas1=[];
for a1=1:numel(alphaN)
    laplas1(a1,:)=s1.^alphaN(a1);
end
results1=[];
sumcoef1=[];
syms('s');
CoN = ufpoly2str(Controller.b,'*');
%CoN = eval(CoN);
if  isscalar(CoN)
    CoN = CoN;
    NofCo = CoN;
else
    CoN = eval(CoN);
    CN = matlabFunction(CoN);
    NofCo = CN(s1);
end
for b1=1:r1
    for d1=1:numel(alphaN)
        sumcoef1(d1,1:numel(w1))=NofCo.*h1(b1,d1).*laplas1(d1,:);
    end
    results1(b1,1:numel(w1))=sum(sumcoef1,1);
end
%% second part
had2=vertcat(lowerboundsD,upperboundsD);
n2=numel(alphaD);
r2=2^n2;
c2=n2;
h2=nan(r2,c2);
inx2=[];
for i2=1:c2
    st2=2^(i2);
    id2=ones(st2/2,r2/(st2/2));
    id2(:,2:2:end)=2;
    idx2=reshape(id2,[r2,1]);
    inx2=horzcat(inx2,idx2);
end

h2=[];

for j2=1:r2
    for k2=1:c2
        h2(j2,k2)=had2(inx2(j2,k2),k2);
    end
end
w2= 0:.01:hhh;
s2= j*w2;

laplas2=[];
for a2=1:numel(alphaD)
    laplas2(a2,:)=s2.^alphaD(a2);
end
results2=[];
sumcoef2=[];
syms('s');
CoD = ufpoly2str(Controller.a,'*');
if  isscalar(CoD)
    CoD=CoD;
    DofCo = CoD;
else
    CoD = eval(CoD);
    CD = matlabFunction(CoD);
    DofCo = CD(s2);
end
for b2=1:r2
    for d2=1:numel(alphaD)
        sumcoef2(d2,1:numel(w2))=DofCo.*h2(b2,d2).*laplas2(d2,:);
    end
    results2(b2,1:numel(w2))=sum(sumcoef2,1);
end
%% third part
nu=1;
nx=1;
for k=1:2:2*N-1
    %:2:2*N-1
    ple = exp(-s1*(l1 + ((k-1)*((l2 - l1)/(2*N)))));
    plo = (1./cos(w1*((l2-l1)/(2*N)))).*exp(-s1*(l1+(k*((l2-l1)/(2*N)))));
    plee = exp(-s1*(l1 + ((k+1)*((l2 - l1)/(2*N)))));
    pse  = exp(-s1.*(l1+(k-1)*((l2-l1)/(2*N))));
    pso  = exp(-s1*(l1+k*((l2-l1)/(2*N))));
    psee = exp(-s1*(l1+(k+1)*((l2-l1)/(2*N))));
    for d2 = 1:r2
        for d1 = 1:r1
            if d2 < r2
                TrEn1(d1,:) =  abs(ple.*results1(d1,:) + results2(d2,:)) + abs(ple.*results1(d1,:) + results2(d2+1,:)) - abs(results2(d2,:) - results2(d2+1,:));
                TrEn2(d1,:) =  abs(plo.*results1(d1,:) + results2(d2,:)) + abs(plo.*results1(d1,:) + results2(d2+1,:)) - abs(results2(d2,:) - results2(d2+1,:));
                TrEn3(d1,:) =  abs(plee.*results1(d1,:) + results2(d2,:)) + abs(plee.*results1(d1,:) + results2(d2+1,:)) - abs(results2(d2,:) - results2(d2+1,:));
                TrEn4(d1,:) =  abs(pse.*results1(d1,:) + results2(d2,:)) + abs(pse.*results1(d1,:) + results2(d2+1,:)) - abs(results2(d2,:) - results2(d2+1,:));
                TrEn5(d1,:) =  abs(pso.*results1(d1,:) + results2(d2,:)) + abs(pso.*results1(d1,:) + results2(d2+1,:)) - abs(results2(d2,:) - results2(d2+1,:));
                TrEn6(d1,:) =  abs(psee.*results1(d1,:) + results2(d2,:)) + abs(psee.*results1(d1,:) + results2(d2+1,:)) - abs(results2(d2,:) - results2(d2+1,:));
            else
                TrEn1(d1,:) =  abs(ple.*results1(d1,:) + results2(1,:)) + abs(ple.*results1(d1,:) + results2(r2,:)) - abs(results2(1,:) - results2(r2,:));
                TrEn2(d1,:) =  abs(plo.*results1(d1,:) + results2(1,:)) + abs(plo.*results1(d1,:) + results2(r2,:)) - abs(results2(1,:) - results2(r2,:));
                TrEn3(d1,:) =  abs(plee.*results1(d1,:) + results2(1,:)) + abs(plee.*results1(d1,:) + results2(r2,:)) - abs(results2(1,:) - results2(r2,:));
                TrEn4(d1,:) =  abs(pse.*results1(d1,:) + results2(1,:)) + abs(pse.*results1(d1,:) + results2(r2,:)) - abs(results2(1,:) - results2(r2,:));
                TrEn5(d1,:) =  abs(pso.*results1(d1,:) + results2(1,:)) + abs(pso.*results1(d1,:) + results2(r2,:)) - abs(results2(1,:) - results2(r2,:));
                TrEn6(d1,:) =  abs(psee.*results1(d1,:) + results2(1,:)) + abs(psee.*results1(d1,:) + results2(r2,:)) - abs(results2(1,:) - results2(r2,:));
            end
            TR1(nu) = min(TrEn1,[],'all');
            TR2(nu) = min(TrEn2,[],'all');
            TR3(nu) = min(TrEn3,[],'all');
            TR4(nu) = min(TrEn4,[],'all');
            TR5(nu) = min(TrEn5,[],'all');
            TR6(nu)  = min(TrEn6,[],'all');
        end
        nu=nu+1;
    end

    for d1 = 1:r1
        for d2 = 1:r2

            TreEn1(d2,:) =  abs(ple.*results1(d1,:) + results2(d2,:)) + abs(plo.*results1(d1,:) + results2(d2,:)) - abs(ple.*results1(d1,:) - plo.*results1(d1,:));
            TreEn2(d2,:) =  abs(plo.*results1(d1,:) + results2(d2,:)) + abs(plee.*results1(d1,:) + results2(d2,:)) - abs(plo.*results1(d1,:) - plee.*results1(d1,:));
            TreEn3(d2,:) =  abs(pse.*results1(d1,:) + results2(d2,:)) + abs(pso.*results1(d1,:) + results2(d2,:)) - abs(pso.*results1(d1,:) - pse.*results1(d1,:));
            TreEn4(d2,:) =  abs(pso.*results1(d1,:) + results2(d2,:)) + abs(psee.*results1(d1,:) + results2(d2,:)) - abs(pso.*results1(d1,:) - psee.*results1(d1,:));
            TeR1(nx) = min(TreEn1,[],'all');
            TeR2(nx) = min(TreEn2,[],'all');
            TeR3(nx) = min(TreEn3,[],'all');
            TeR4(nx) = min(TreEn4,[],'all');
        end
        nx=nx+1;
    end
end
minimumval1= min ([min(TR1) min(TR2) min(TR3) min(TR4) min(TR5) min(TR6) min(TeR1) min(TeR2) min(TeR3) min(TeR4)]);
%% forth part
nuu=1;
nxx=1;
for k=1:2:2*N-1
    %:2:2*N-1
    gle = exp(s1*(l1 + ((k-1)*((l2 - l1)/(2*N)))));
    glo = (1./cos(w1*((l2-l1)/(2*N)))).*exp(s1*(l1+(k*((l2-l1)/(2*N)))));
    glee = exp(s1*(l1 + ((k+1)*((l2 - l1)/(2*N)))));
    gse  = exp(s1.*(l1+(k-1)*((l2-l1)/(2*N))));
    gso  = exp(s1*(l1+k*((l2-l1)/(2*N))));
    gsee = exp(s1*(l1+(k+1)*((l2-l1)/(2*N))));
    for d1 = 1:r1
        for d2 = 1:r2
            if d1 < r1
                GrEn1(d2,:) =  abs(gle.*results2(d2,:) + results1(d1,:)) + abs(gle.*results2(d2,:) + results1(d1+1,:)) - abs(results1(d1,:) - results1(d1+1,:));
                GrEn2(d2,:) =  abs(glo.*results2(d2,:) + results1(d1,:)) + abs(glo.*results2(d2,:) + results1(d1+1,:)) - abs(results1(d1,:) - results1(d1+1,:));
                GrEn3(d2,:) =  abs(glee.*results2(d2,:) + results1(d1,:)) + abs(glee.*results2(d2,:) + results1(d1+1,:)) - abs(results1(d1,:) - results1(d1+1,:));
                GrEn4(d2,:) =  abs(gse.*results2(d2,:) + results1(d1,:)) + abs(gse.*results2(d2,:) + results1(d1+1,:)) - abs(results1(d1,:) - results1(d1+1,:));
                GrEn5(d2,:) =  abs(gso.*results2(d2,:) + results1(d1,:)) + abs(gso.*results2(d2,:) + results1(d1+1,:)) - abs(results1(d1,:) - results1(d1+1,:));
                GrEn6(d2,:) =  abs(gsee.*results2(d2,:) + results1(d1,:)) + abs(gsee.*results2(d2,:) + results1(d1+1,:)) - abs(results1(d1,:) - results1(d1+1,:));
            else
                GrEn1(d2,:) =  abs(gle.*results2(d2,:) + results1(1,:)) + abs(gle.*results2(d2,:) + results1(r1,:)) - abs(results1(1,:) - results1(r1,:));
                GrEn2(d2,:) =  abs(glo.*results2(d2,:) + results1(1,:)) + abs(glo.*results2(d2,:) + results1(r1,:)) - abs(results1(1,:) - results1(r1,:));
                GrEn3(d2,:) =  abs(glee.*results2(d2,:) + results1(1,:)) + abs(glee.*results2(d2,:) + results1(r1,:)) - abs(results1(1,:) - results1(r1,:));
                GrEn4(d2,:) =  abs(gse.*results2(d2,:) + results1(1,:)) + abs(gse.*results2(d2,:) + results1(r1,:)) - abs(results1(1,:) - results1(r1,:));
                GrEn5(d2,:) =  abs(gso.*results2(d2,:) + results1(1,:)) + abs(gso.*results2(d2,:) + results1(r1,:)) - abs(results1(1,:) - results1(r1,:));
                GrEn6(d2,:) =  abs(gsee.*results2(d2,:) + results1(1,:)) + abs(gsee.*results2(d2,:) + results1(r1,:)) - abs(results1(1,:) - results1(r1,:));
            end
            GR1(nuu) = min(GrEn1,[],'all');
            GR2(nuu) = min(GrEn2,[],'all');
            GR3(nuu) = min(GrEn3,[],'all');
            GR4(nuu) = min(GrEn4,[],'all');
            GR5(nuu) = min(GrEn5,[],'all');
            GR6(nuu)  = min(GrEn6,[],'all');
        end
        nuu=nuu+1;
    end

    for d1 = 1:r1
        for d2 = 1:r2

            GreEn1(d2,:) =  abs(gle.*results2(d2,:) + results1(d1,:)) + abs(glo.*results2(d2,:) + results1(d1,:)) - abs(gle.*results2(d2,:) -glo.*results2(d2,:));
            GreEn2(d2,:) =  abs(glo.*results2(d2,:) + results1(d1,:)) + abs(glee.*results2(d2,:) + results1(d1,:)) - abs(glo.*results2(d2,:) - glee.*results2(d2,:));
            GreEn3(d2,:) =  abs(gse.*results2(d2,:) + results1(d1,:)) + abs(gso.*results2(d2,:) + results1(d1,:)) - abs(gso.*results2(d2,:) - gse.*results2(d2,:));
            GreEn4(d2,:) =  abs(gso.*results2(d2,:) + results1(d1,:)) + abs(gsee.*results2(d2,:) + results1(d1,:)) - abs(gso.*results2(d2,:) - gsee.*results2(d2,:));
            GeR1(nxx) = min(GreEn1,[],'all');
            GeR2(nxx) = min(GreEn2,[],'all');
            GeR3(nxx) = min(GreEn3,[],'all');
            GeR4(nxx) = min(GreEn4,[],'all');
        end
        nxx=nxx+1;
    end
end
minimumval2= min ([min(GR1) min(GR2) min(GR3) min(GR4) min(GR5) min(GR6) min(GeR1) min(TeR2) min(GeR3) min(GeR4)]);
Auxi1 = min([minimumval1 minimumval2]);
% hh= sum(max(abs(lowerbounds(1:end-1)),abs(upperbounds(1:end-1))))/(min(abs(lowerbounds(end)),abs(upperbounds(end))));
% hhh=hh^(1/(alpha(end)-alpha(end-1)));
if 2*pi/(l2-l1)> wmax
    if Auxi1>.0001
        display("the system is robustly stable")
    else
        if  l2==0
            display("the system is not robustly stable")
        elseif N<15
            display ("increse N")
        elseif  N>=15
            display("the system is not robustly stable")
        end
    end
else
    w3= 2*pi/(l2-l1):.01:wmax;
    s3= j*w3;
    laplas3=[];
    for a3=1:numel(alphaN)
        laplas3(a3,:)=s3.^alphaN(a3);
    end
    results3=[];
    sumcoef3=[];
    syms('s');
    CoN3 = ufpoly2str(Controller.b,'*');
    if isscalar(CoN3)
        Con3=CoN3;
        NofCo3 = CoN3;
    else
        CoN3 = eval(CoN3);
        CN3 = matlabFunction(CoN3);
        NofCo3 = CN3(s3);
    end
    for b1=1:r1
        for d1=1:numel(alphaN)
            sumcoef3(d1,1:numel(w3))=NofCo3.*h1(b1,d1).*laplas3(d1,:);
        end
        results3(b1,1:numel(w3))=sum(sumcoef3,1);
    end
    s3= j*w3;
    laplas4=[];
    for a2=1:numel(alphaD)
        laplas4(a2,:)=s3.^alphaD(a2);
    end
    results4=[];
    sumcoef4=[];
    syms('s');
    CoD4 = ufpoly2str(Controller.a,'*');
    if isscalar(CoD4)
        CoD4 = CoD4;
        DofCo4 = CoD4;
    else
    CoD4 = eval(CoD4);
    CD4 = matlabFunction(CoD4);
    DofCo4 = CD4(s3);
    end
    for b2=1:r2
        for d2=1:numel(alphaD)
            sumcoef4(d2,1:numel(w3))=DofCo4.*h2(b2,d2).*laplas4(d2,:);
        end
        results4(b2,1:numel(w3))=sum(sumcoef4,1);
    end
    for i1 =1:length(w3)
        lan = 0:.01:1;

        for d0 = 1:r2
            for d1 = 1:r2

                E(d1,:) =  abs(lan* results4(d0,i1)+(1-lan)*results4(d1,i1));

            end
            ED(d0,i1) =  min(E(d1,:));
        end

        EN(:,i1) = max([results3(:,i1)]);
        Auxi2(:,i1) =  min(ED(:,i1)) - EN(:,i1);
    end

    Auxi = min ([Auxi1 Auxi2]);
    if Auxi>.0001
        display("the system is robustly stable")
    else
        if l2==0
            display("the system is not robustly stable")
        elseif N<15
            display ("increse N")
        elseif  N>=15
            display("the system is not robustly stable")
        end
    end
end
