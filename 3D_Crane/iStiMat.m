function [ID,iGloStiMatrix] = iStiMat(Ac,E,RowriMatrix)
LocMat=((Ac*E)/RowriMatrix(2))*[1 -1;-1 1];
M=zeros(2,6);
M(1,1:3)=RowriMatrix(3:5);
M(2,4:6)=RowriMatrix(3:5);

ID=RowriMatrix(1);
iGloStiMatrix=M'*LocMat*M;
end