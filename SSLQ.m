function [T,SS,QC,Ueqn,algorithm] = SSLQ(A,B,Qb,AlphaGain)
% [T,SS,QC,Ueqn,algorithm] = SSLQ(A,B,Qb,Gain)
% SingularLQ: 1) Calculates transformation T from matrices A and B,  to a Controllable Canonical Form:
% z = Tx,
% 2) Depending on Weithing Matrix Qb, calculates the Sliding Surface and its
% r-1 time derivatives [s \dot{s} \ddot{s} ...]^{T} = SS*z, 
% 3) Generates Quasi-Continuous Algorithm of order r with the specified Gain,
% 4) Calculates nominal equivalent control ueqn = Ueqn*z, so the complete control law is u = ueqn
% + QC


%% Tests
n = TestSLQ(A,B,Qb,AlphaGain);

%% Transformation

[T,Ab,Bb] = TransCC(A,B,n);

%% Orden de Singularidad

r = OrderSingular(Qb,n);

%% Particiones y K

K = PartK(Ab,Qb,r,n);

%% Construcción de SS

SS = Surfaces(K,r,n);

%% Control equivalente

Ueqn = -Ab(n,:) - [zeros(1,r) K];

%% Quasi-Continuous

[QC,S] = QuasiC(r);

%% Algorithm

algorithm = Algthm(T,SS,QC,Ueqn,r,AlphaGain,S);

end

function n = TestSLQ(A,B,Qb,AlphaGain)
%% Tests

[mA,nA] = size(A);
[mB,nB] = size(B);
n = 0;
if(mA~=nA)
    error('A matrix must be square');
else
   n = nA;
end

if(mB~=n)
    error('The A matrix must be square with as many rows as B');
end

[mQb,nQb] = size(Qb);

if(mQb~=nQb)
    error('Qb matrix must be square');
end

if(mQb~=n)
    error('The Qb matrix must be square with as many rows as A');
end

if(AlphaGain<0)
    error('The Alpha Gain must be positive');
end

end

function [T,Ab,Bb] = TransCC(A,B,n)

C = ctrb(A,B);
rnk = rank(C);

if(n~=rnk)
    error('Controllability Error: Pair (A,B) must be controllable.');
end

Cinv = eye(rnk)/C;
q = Cinv(rnk,:);
clear T;
for i=1:rnk
T(i,:) = q*A^(i-1);
end

Ab = T*A/T;
Bb = T*B;

end

function r = OrderSingular(Qb,n)


clear r;
for i=1:n

    r = (i-1) +1;
    if(  Qb(n-(i-1),n-(i-1))  )~=0
       break    
    end
         
end

end

function K = PartK(Ab,Qb,r,n)


if(r<n)
Q11 = Qb(1:n-r,1:n-r);
Q22 = Qb(n-(r-1),n-(r-1));
Q12 = Qb(1:n-r,n-(r-1));

A11 = Ab(1:n-r,1:n-r);
A22 = Ab(n-(r-1),n-(r-1));
A12 = Ab(1:n-r,n-(r-1));

Ah = A11 + A12*Q22^(-1)*Q12';
Bh = A12;

Qh = Q11 -Q12*Q22^(-1)*Q12';
Rh = Q22;

Db = cholcov(Qh);
O = obsv(Ah,Db);

if((n-r)~=rank(O))
    error('Observability Error: Pair (Ah,Db), Db^(T)DB=Qh, must be observable.');
end

P = are(Ah,Bh*inv(Rh)*Bh',Qh);
K = inv(Q22)*(A12'*P + Q12');


elseif (n==r)
    
    K = [];
    
else
    
end

end

function SS = Surfaces(K,r,n)

SS = zeros(r,n);

for i=1:r

    SS(i,:) = [ zeros(1,i-1) K 1 zeros(1,r-i) ];
    
end

end

function [QC,S] = QuasiC(r)
%% Quasi-Continuous

beta = [1   0 0 0 0 0 0 0 0 0;
        1   0 0 0 0 0 0 0 0 0;
        1   2 0 0 0 0 0 0 0 0;
        0.5 1 3 0 0 0 0 0 0 0;
        0.4 1 1.8 4 0 0 0 0 0 0;
        0.2 0.8 1.1 3.2 5 0 0 0 0 0; %Probados Hasta 6
        0.4 3 7 11 14 17 0 0 0 0; 
        0.3 4 8 13 17 21 15 0 0 0;
        0.2 5 9 15 19 24 28 31 0 0;
        0.2 6 10 16 20 25 29 34 40 0];
    
syms s sd sdd sddd sdddd sddddd;
S = [s sd sdd sddd sdddd sddddd];
syms Alpha

phi(1) = s;
N(1) = abs(s);
PHI(1) = phi(1)/N(1);

if r>1 
    for i = 2:r
        phi(i) = S(i) + beta(r,i-1)*N(i-1)^(-1/(r-(i-1)+1))*phi(i-1);
        N(i) = abs(S(i)) + beta(r,i-1)*N(i-1)^(-1/(r-(i-1)+1))*abs(phi(i-1));
        PHI(i) = phi(i)/N(i);
    end
end

QC = -Alpha* PHI(r);



end

function algorithm = Algthm(T,SS,QC,Ueqn,r,AlphaGain,S)
%% Algorithm

ssurf = '';
for i=1:r

ssurf = [ssurf, char(S(i)) ,'=ss(',num2str(i),');',char(10)];
    
end

algorithm =  ['T= ', mat2str(T),'; ',char(10),'z=T*x;',char(10),'SS=',mat2str(SS),';',char(10),'ss = SS*z;',char(10),'Ueqn=',mat2str(Ueqn),';',char(10),'ueqn = Ueqn*z;',char(10),'Alpha = ',num2str(AlphaGain),';',char(10),ssurf,'QC=',char(QC),';',char(10),'u=ueqn+QC;',char(10)];


end
