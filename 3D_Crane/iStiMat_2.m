function [ID,iGloStiMatrix] = iStiMat_2(L,Ac,E,RowriMatrix)
LocMat=zeros(12,12);
if RowriMatrix(2)<2*L
    J=((Ac)^2)/(2*pi);
    I=((Ac)^2)/(4*pi);
    a=RowriMatrix(2)/2;
    G=79300;
    tor=G*J/(2*a);
    bend=E*I/a;
else 
     J=0;
    I=0;
    a=RowriMatrix(2)/2;
    G=79300;
    tor=G*J/(2*a);
    bend=E*I/a;
end

LocMat(1,[1 7])=((Ac*E)/RowriMatrix(2))*[1 -1];
LocMat(2,[2 6 8 12])=[(3/(2*a^2))*bend (3/(2*a))*bend -(3/(2*a^2))*bend (3/(2*a))*bend];
LocMat(3,[3 5 9 11])=[(3/(2*a^2))*bend -(3/(2*a))*bend -(3/(2*a^2))*bend -(3/(2*a))*bend];
LocMat(4,[4 10])=[tor -tor];
LocMat(5,[5 9 11])=[2*bend (3/(2*a))*bend bend];
LocMat(6,[6 8 12])=[2*bend -(3/(2*a))*bend bend];
LocMat(7,[7])=[(Ac*E)/RowriMatrix(2)];
LocMat(8,[8 12])=[(3/(2*a^2))*bend -(3/(2*a^2))*bend];
LocMat(9,[9 11])=[(3/(2*a^2))*bend (3/(2*a))*bend];
LocMat(10,[10])=[tor];
LocMat(11,[11])=[2*bend];
LocMat(12,[12])=[2*bend];

LMdiag=diag(LocMat);
LocMat=LocMat+LocMat';
LocMat(1:size(LocMat,1)+1:end)=LMdiag;

T3=[RowriMatrix(4:6); RowriMatrix(7:9); RowriMatrix(10:12)];
T=blkdiag(T3,T3,T3,T3);
iGloStiMatrix=T'*LocMat*T;

ID=RowriMatrix(1);
end