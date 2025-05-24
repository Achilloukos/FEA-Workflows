%Fasoulas Achillefs AEM:6581

%This is the preprocessor for the truss elements problem

%Initial Data

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

%Creation of the Matrix that contains the coordinates of the Crane's Nodes before rotating



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

%Creation of the matrix that contains the IDs of the Nodes for each element
EleMatrix=elmnt2(NodeMatrix,L);

% Coordinate Transformation (rotating the structure around y and z axis)
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

%Boundary Conditions (without the MPCs yet)
BCMatrix=zeros(size(NodeMatrix));
BCMatrix(:,1)=NodeMatrix(:,1);
BCMatrix([1 26 27 28],2:4)=1;
BCMatrix(2,4)=1;

%Load Vector
LoadVector=zeros(84,1);
LoadVector(75)=-20000;
%Creation of a matrix that will be used for the transformation from local stifness matrix to global, later

OriMatrix=orient(NodeMatrix,EleMatrix,A0,L);

%Creation of .txt file
sizeNodeMatrix=size(NodeMatrix);
sizeEleMatrix=size(EleMatrix);
sizeOriMatrix=size(OriMatrix);
sizeBCMatrix=size(BCMatrix);
sizeLoadVector=size(LoadVector);
NodeMatrix(end+1:end+((sizeEleMatrix(1))-sizeNodeMatrix(1)),:)=0;
BCMatrix(end+1:end+((size(EleMatrix,1)-size(BCMatrix,1))),:)=0;
LoadVector(end+1:end+((sizeEleMatrix(1))-sizeLoadVector(1)),:)=0;
Sizes=[sizeNodeMatrix sizeEleMatrix sizeOriMatrix sizeBCMatrix sizeLoadVector zeros(1,92)]';
Table=table(NodeMatrix,EleMatrix,OriMatrix,BCMatrix,LoadVector,Sizes);

writetable(Table,'PreText.txt');

clear all
