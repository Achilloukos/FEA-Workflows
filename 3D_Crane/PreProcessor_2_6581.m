%Fasoulas Achillefs AEM:6581

%A=6;
%B=5;
%C=8;
%D=1;
phi=60+-0.9; %deg, y axis rotation
theta=45+6.8;  %deg, z axis rotation
L=1.5*1.18*1000; %mm
n=5+mod(6,1); %floorN
A0=6*(0.5+18/100)*100; %mm^2
Ad=0.5*A0; %mm^2
Ar=1.5*A0; %mm^2
A1=1.2*1.65*1000;
B1=1.58*1000;

floorN= n;
f=1;

NodeMatrixstd = standartfloornodes(floorN,f,L);
NodeMatrixstd(1,:)=[0 0 -L/2 0];
NodeMatrixstd(2,:)=[0 0 L/2 0];
NodeMatrixstd(23,:)=[23 L/2 -L/2 6*L];
NodeMatrixstd(24,:)=[24 L/2 +L/2 6*L];
NodeMatrixstd(25,:)=[25 1.5*L 0 5.5*L];
NodeMatrixstd(:,1)=[1:size(NodeMatrixstd,1)];

NodeMatrixstd(26:28,1)=[26:28];
NodeMatrixstd(26:28,4)=B1;
NodeMatrixstd(26:28,2)=-A1;
NodeMatrixstd(26,3)=-0.5*L;
NodeMatrixstd(27,3)=0;
NodeMatrixstd(28,3)=0.5*L;

NodeMatrix=NodeMatrixstd;

EleMatrix=elmnt2(NodeMatrix,L)
EleMatrix2=elmnt_2(NodeMatrix,L);

MatRot1=RotMat(0,90-phi,90-theta);
NodeMatrixRot=NodeMatrix;
for i=1:25
   NodeMatrixRot(i,2:4)=MatRot1*[NodeMatrix(i,2); NodeMatrix(i,3); NodeMatrix(i,4)];
end
MatRot2=RotMat(0,0,90-theta);
NodeMatrixRot(26,2:4)= MatRot2*[NodeMatrix(26,2); NodeMatrix(26,3); NodeMatrix(26,4)];  
NodeMatrixRot(27,2:4)= MatRot2*[NodeMatrix(27,2); NodeMatrix(27,3); NodeMatrix(27,4)];
NodeMatrixRot(28,2:4)= MatRot2*[NodeMatrix(28,2); NodeMatrix(28,3); NodeMatrix(28,4)];

NodeMatrix=NodeMatrixRot;

BCMatrix2=zeros(size(NodeMatrix,1),size(NodeMatrix,2)+3);
BCMatrix2(:,1)=NodeMatrix(:,1);
BCMatrix2([1 26 27 28],2:4)=1;
BCMatrix2(2,4)=1;
BCMatrix2([26 27 28],5:7)=1;

LoadVector2=zeros(168,1);
LoadVector2(147)=-20000;

OriMatrix=orient(NodeMatrix,EleMatrix,A0,L);
OriMatrix2=orient_2(NodeMatrix,EleMatrix2,A0,L);
Abeam = (sum(OriMatrix(1:(end-4),2).*OriMatrix(1:(end-4),6)))/sum(OriMatrix2(1:(end-4),2));
OriMatrix2(1:(end-4),3)=Abeam;


%txt table
sizeNodeMatrix=size(NodeMatrix);
sizeEleMatrix2=size(EleMatrix2);
sizeOriMatrix2=size(OriMatrix2);
sizeBCMatrix2=size(BCMatrix2);
sizeLoadVector2=size(LoadVector2);

OriMatrix2(end+1:end+((sizeLoadVector2(1))-sizeOriMatrix2(1)),:)=0;
EleMatrix2(end+1:end+((sizeLoadVector2(1))-sizeEleMatrix2(1)),:)=0;
NodeMatrix(end+1:end+((sizeLoadVector2(1))-sizeNodeMatrix(1)),:)=0;
BCMatrix2(end+1:end+((sizeLoadVector2(1)-size(BCMatrix2,1))),:)=0;
LoadVector2(end+1:end+((sizeEleMatrix2(1))-sizeLoadVector2(1)),:)=0;

Sizes=[sizeNodeMatrix sizeEleMatrix2 sizeOriMatrix2 sizeBCMatrix2 sizeLoadVector2 zeros(1,158)]';
Table=table(NodeMatrix,EleMatrix2,OriMatrix2,BCMatrix2,LoadVector2,Sizes);

writetable(Table,'PreText_2.txt');

clear all