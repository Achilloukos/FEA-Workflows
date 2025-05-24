%Fasoulas Achillefs AEM:6581

theta=45+6.8;

%Importing data from preprocessor
PreCell=DataImport('PreText.txt');
NodeMatrix=PreCell{1};
EleMatrix=PreCell{2};
OriMatrix=PreCell{3};
BCMatrix=PreCell{4};
LoadVector=PreCell{5};

A0=max(OriMatrix(:,6))/1.5;
L=min(OriMatrix(:,2));
theta=45+6.8;

%LoadVector to Matrix
LoadMat(:,1)=NodeMatrix(:,1);
LoadMat(:,2:4)=[LoadVector(1:3:82) LoadVector(2:3:83) LoadVector(3:3:84)];

%Creation of local stiffness matrix
StiffCell=cell(size(OriMatrix,1),2);
for i=1:size(OriMatrix,1)
    [StiffCell{i,1}, StiffCell{i,2}]=iStiMat(OriMatrix(i,6),210000,OriMatrix(i,1:5));
end

%Creation of global stiffness matrix
GloStiMatrix=zeros(3*size(NodeMatrix,1));

for i=1:size(StiffCell,1)
    IDi=StiffCell(i,1);
    StiMi=StiffCell{i,2};
    Nodesi=[NodeMatrix(EleMatrix(i,2),1) NodeMatrix(EleMatrix(i,3),1)];
    DoFsi=[(3*Nodesi(1)-2):3*Nodesi(1) (3*Nodesi(2)-2):3*Nodesi(2)];
    for j=1:size(StiMi,1)
        for k=1:size(StiMi,2)
            GloStiMatrix(DoFsi(j),DoFsi(k))=StiMi(j,k)+GloStiMatrix(DoFsi(j),DoFsi(k));
        end
    end
end

GloStiMatrixForces=GloStiMatrix; %This will be used for the calculation of reaction forces

%BCVector to Matrix
BCVector=zeros(84,1);
BCVector(1:3:82)=BCMatrix(:,2);
BCVector(2:3:83)=BCMatrix(:,3);
BCVector(3:3:84)=BCMatrix(:,4);

%BCVector(4:5)=1;
DoF=1:3*size(NodeMatrix,1);

ElimiVector=find(BCVector);
GloStiMatrix(ElimiVector,:)=[];
GloStiMatrix(:,ElimiVector)=[];
LoadVector(ElimiVector)=[];
DoF(ElimiVector)=[];


DeformTemp=GloStiMatrix\LoadVector;
%Multi Point Constraints and solving the linear system
T=eye(size(GloStiMatrix,1));

T(2,1)=-1/tand(90-theta);
T(:,2)=[];

GloStiMatrix=T'*GloStiMatrix*T;
LoadVector=T'*LoadVector;
DeformTempNew=GloStiMatrix\LoadVector;

%Displacements calculation
DisplacementsVector=zeros(3*size(NodeMatrix,1),1);
DisplacementsVector(5)=-(1/tand(90-theta))*DeformTempNew(1);
BCVector(5)=1;
DisplacementsVector(find(BCVector==0))=DeformTempNew(:);
 
%Reaction Forces calculation
ReactionForcesDoF=[ElimiVector; 4; 5];
GloStiMatrixForces=GloStiMatrixForces(ReactionForcesDoF,:);
ReactionForces=GloStiMatrixForces*DisplacementsVector;
AllReactionForces=zeros(84,1);
AllReactionForces(ReactionForcesDoF)=ReactionForces;
ReactionForcesMat=zeros(28,4);
ReactionForcesMat(:,1)=NodeMatrix(:,1);
ReactionForcesMat(:,2:4)=[AllReactionForces(1:3:82) AllReactionForces(2:3:83) AllReactionForces(3:3:84)]; 

%Calculation of the deformation (deltaL/L0) of each element
DeforMatrix=zeros(size(NodeMatrix,1),4);
DeforMatrix(:,1)=NodeMatrix(:,1);
DeforMatrix(:,2)=DisplacementsVector(1:3:3*size(NodeMatrix,1)-2);
DeforMatrix(:,3)=DisplacementsVector(2:3:3*size(NodeMatrix,1)-1);
DeforMatrix(:,4)=DisplacementsVector(3:3:3*size(NodeMatrix,1));

NodeMatrixDef=NodeMatrix;
NodeMatrixDef(:,2:4)=NodeMatrixDef(:,2:4)+DeforMatrix(:,2:4);
OriMatrixDef=orient(NodeMatrixDef,EleMatrix,A0,L);
DeformElement=-(OriMatrix(:,1:2)-OriMatrixDef(:,1:2));

DeformElementPercentage=zeros(size(EleMatrix,1),2); %DeformElementPercentage is the DeltaL/L0 of each element
DeformElementPercentage(:,1)=EleMatrix(:,1);
DeformElementPercentage(:,2)=DeformElement(:,2)./OriMatrix(:,2);

%Calculation of the normal stress of each element
EleSigma=zeros(size(EleMatrix,1),2);
EleSigma(:,1)=EleMatrix(:,1);
EleSigma(:,2)=DeformElementPercentage(:,2)*210000;

%Creation of .txt file
sizeReactionForcesMat=size(ReactionForcesMat);
sizeEleSigma=size(EleSigma);
sizeDeforMatrix=size(DeforMatrix);
sizeDeformElementPercentage=size(DeformElementPercentage);

ReactionForcesMat(end+1:end+((sizeEleSigma(1))-sizeReactionForcesMat(1)),:)=0;
DeforMatrix(end+1:end+(sizeEleSigma(1)-sizeDeforMatrix(1)),:)=0;

Sizes=[sizeReactionForcesMat sizeEleSigma  sizeDeforMatrix sizeDeformElementPercentage zeros(1,94)]';

Table=table(ReactionForcesMat,EleSigma,DeforMatrix,DeformElementPercentage,Sizes);

writetable(Table,'SolText.txt');

clear all
