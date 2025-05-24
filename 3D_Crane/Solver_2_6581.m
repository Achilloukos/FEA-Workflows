%Fasoulas Achillefs AEM:6581

PreCell=DataImport('PreText_2.txt');
NodeMatrix=PreCell{1};
EleMatrix2=PreCell{2};
OriMatrix2=PreCell{3};
BCMatrix2=PreCell{4};
LoadVector2=PreCell{5};

A0=max(OriMatrix2(:,6))/1.5;
L=min(OriMatrix2(:,2));
theta=45+6.8;

StiffCell=cell(size(OriMatrix2,1),2);
for i=1:size(OriMatrix2,1)
    [StiffCell{i,1}, StiffCell{i,2}]=iStiMat_2(L,OriMatrix2(i,3),210000,OriMatrix2(i,1:12));
end

%Global Stifness Matrix calculation

GloStiMatrix2=zeros(6*size(NodeMatrix,1));

for i=1:size(StiffCell,1)
    IDi=StiffCell(i,1);
    StiMi=StiffCell{i,2};
    Nodesi=[NodeMatrix(EleMatrix2(i,2),1) NodeMatrix(EleMatrix2(i,3),1)];
    DoFsi=[(6*Nodesi(1)-5):6*Nodesi(1) (6*Nodesi(2)-5):6*Nodesi(2)];
    for j=1:size(StiMi,1)
        for k=1:size(StiMi,2)
            GloStiMatrix2(DoFsi(j),DoFsi(k))=StiMi(j,k)+GloStiMatrix2(DoFsi(j),DoFsi(k));
        end
    end
end

GloStiMatrixForces2=GloStiMatrix2; %This will be used for the calculation of reaction forces

BCVector2=zeros(168,1);
BCVector2(1:6:163)=BCMatrix2(:,2);
BCVector2(2:6:164)=BCMatrix2(:,3);
BCVector2(3:6:165)=BCMatrix2(:,4);
BCVector2(4:6:166)=BCMatrix2(:,5);
BCVector2(5:6:167)=BCMatrix2(:,6);
BCVector2(6:6:168)=BCMatrix2(:,7);

DoF=1:6*size(NodeMatrix,1);

ElimiVector=find(BCVector2);
GloStiMatrix2(ElimiVector,:)=[];
GloStiMatrix2(:,ElimiVector)=[];

LoadVector2(ElimiVector)=[];
DoF(ElimiVector)=[];

%MPC using master-slave method

DeformTemp=GloStiMatrix2\LoadVector2;
T2=eye(size(GloStiMatrix2,1));

T2(8-3,7-3)=-1/tand(90-theta);
T2(:,8-3)=[];

GloStiMatrix2=T2'*GloStiMatrix2*T2;
LoadVector2=T2'*LoadVector2;
DeformTempNew2=GloStiMatrix2\LoadVector2;

%Node Displacements Calculation
DisplacementsVector2=zeros(6*size(NodeMatrix,1),1);
DisplacementsVector2(8)=-(1/tand(90-theta))*DeformTempNew2(4);
BCVector2(8)=1;
DisplacementsVector2(find(BCVector2==0))=DeformTempNew2(:)
DispMat2=zeros(28,7);
DispMat2(:,1)=NodeMatrix(:,1);
DispMat2(:,2:4)=[DisplacementsVector2(1:6:163) DisplacementsVector2(2:6:164) DisplacementsVector2(3:6:165)] 
DispMat2(:,5:7)=[DisplacementsVector2(4:6:166) DisplacementsVector2(5:6:167) DisplacementsVector2(6:6:168)] 
%Reaction Forces calculation
ReactionForcesDoF2=[ElimiVector; 7; 8];
GloStiMatrixForces2=GloStiMatrixForces2(ReactionForcesDoF2,:);
ReactionForces2=GloStiMatrixForces2*DisplacementsVector2;
AllReactionForces2=zeros(84,1);
AllReactionForces2(ReactionForcesDoF2)=ReactionForces2;
ReactionForcesMat2=zeros(28,7);
ReactionForcesMat2(:,1)=NodeMatrix(:,1);
ReactionForcesMat2(:,2:4)=[AllReactionForces2(1:6:163) AllReactionForces2(2:6:164) AllReactionForces2(3:6:165)]; 
ReactionForcesMat2(:,5:7)=[AllReactionForces2(4:6:166) AllReactionForces2(5:6:167) AllReactionForces2(6:6:168)];

%txt table
sizeReactionForcesMat2=size(ReactionForcesMat2);
sizeDispMat2=size(DispMat2);

% ReactionForcesMat2(end+1:end+((sizeDispMat2(1))-sizeReactionForcesMat(1)),:)=0;
% DeforMatrix(end+1:end+(sizeDispMat2(1)-sizeDeforMatrix(1)),:)=0;

Sizes=[sizeReactionForcesMat2 sizeDispMat2 zeros(1,24)]';

Table=table(ReactionForcesMat2,DispMat2,Sizes);

writetable(Table,'SolText_2.txt');
