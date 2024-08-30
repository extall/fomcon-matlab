function is_robuststable = robstabfo2(P1,P2)
lowerbounds1 = fliplr(P1.a(:,1)');
upperbounds1 = fliplr(P1.a(:,2)');
alpha1 =fliplr(P1.na(:,1)');
lowerbounds2 = fliplr(P2.a(:,1)');
upperbounds2 = fliplr(P2.a(:,2)');
alpha2 =fliplr(P2.na(:,1)');

had1=vertcat(lowerbounds1,upperbounds1);
n1=numel(alpha1);
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
B1 = max(max(abs(lowerbounds1),abs(upperbounds1)));
B2 = max(max(abs(lowerbounds2),abs(upperbounds2)));
A1 = min(min(abs(lowerbounds1),abs(upperbounds1)));
A2 = min(min(abs(lowerbounds2),abs(upperbounds2)));
hh = max(B1,B2)/(min(A1,A2));
alpha = union(alpha1,alpha2); 
hhh=hh^(1/(alpha(end)-alpha(end-1)));
w1= 0:.01:hhh;
s1= j*w1;
laplas1=[];
for a1=1:numel(alpha1)
    laplas1(a1,:)=s1.^alpha1(a1);
end

results1=[];
sumcoef1=[];
for b1=1:r1
    for d1=1:numel(alpha1)
        sumcoef1(d1,1:numel(w1))=h1(b1,d1).*laplas1(d1,:);
    end
    results1(b1,1:numel(w1))=sum(sumcoef1,1);
end

%% second part

had2=vertcat(lowerbounds2,upperbounds2);
n2=numel(alpha2);
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
for a2=1:numel(alpha2)
    laplas2(a2,:)=s2.^alpha2(a2);
end

results2=[];
sumcoef2=[];
for b2=1:r2
    for d2=1:numel(alpha2)
        sumcoef2(d2,1:numel(w2))=h2(b2,d2).*laplas2(d2,:);
    end
    results2(b2,1:numel(w2))=sum(sumcoef2,1);
end

%% robust stability testing function
for d1 = 1:r1
    for d2 = 1:r2
        if d1 < r1
            TrEn(d2,:) =  abs(results1(d1,:) + results2(d2,:)) + abs(results1(d1+1,:) + results2(d2,:)) - abs(results1(d1,:) - results1(d1+1,:));
        else
            TrEn(d2,:) =  abs(results1(1,:) + results2(d2,:)) + abs(results1(r1,:) + results2(d2,:)) - abs(results1(1,:) - results1(r1,:));
        end
    end
    TrEnq1(d1,:) =  min(TrEn(d2,:));
end

RC1 = min(TrEnq1,[],'all');

for d2 = 1:r2
    for d1 = 1:r1
        if d2 < r2
            TrE(d1,:) =  abs(results1(d1,:) + results2(d2,:)) + abs(results1(d1,:) + results2(d2+1,:)) - abs(results2(d2,:) - results2(d2+1,:));
        else
            TrE(d1,:) =  abs(results1(d1,:) + results2(1,:)) + abs(results1(d1,:) + results2(r2,:)) - abs(results2(1,:) - results2(r2,:));
        end
    end
    TrEnq2(d2,:) =  min(TrE(d1,:));
end

RC2 = min(TrEnq2,[],'all');

RC = min(RC1,RC2);

if RC<.0001
    display("not robust stable")
else
    display ("robustly stable")
end