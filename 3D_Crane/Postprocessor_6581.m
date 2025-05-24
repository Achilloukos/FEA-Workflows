%Fasoulas Achillefs AEM:6581

%Importing data from preprocessor and solver
SolCell=DataImport('PreText.txt');
NodeMatrix=SolCell{1};
EleMatrix=SolCell{2};
OriMatrix=SolCell{3};
BCMatrix=SolCell{4};
LoadVector=SolCell{5};

SolCell=DataImport('SolText.txt');
ReactionForcesMat=SolCell{1};
EleSigma=SolCell{2};
DeforMatrix=SolCell{3};
DeformElementPercentage=SolCell{4};

LoadMat(:,1)=NodeMatrix(:,1);
LoadMat(:,2:4)=[LoadVector(1:3:82) LoadVector(2:3:83) LoadVector(3:3:84)];

%Deformed structure using a scale factor of 75 for visual reasons
DefNodeMatrix=zeros(size(NodeMatrix));
DefNodeMatrix(:,1)=NodeMatrix(:,1);
DefNodeMatrix(:,2:4)=NodeMatrix(:,2:4)+75*DeforMatrix(:,2:4);


Sigma=EleSigma(:,2);
SigmaRound=round(Sigma);
Colormap=jet(2*max(SigmaRound)+1);
NodeMatrix=NodeMatrix(1:28,:);

%-Time for some 3D-visualization finaly-%
figure (1)
colormap jet
cb=colorbar;
title(cb,'Stress [N/mm^{2}]')
cb.TickLabels=num2cell(round(linspace(-max(SigmaRound),max(SigmaRound),11)));
plot3(NodeMatrix(:,2),NodeMatrix(:,3),NodeMatrix(:,4),'.')
text(NodeMatrix(:,2),NodeMatrix(:,3),NodeMatrix(:,4),num2str(NodeMatrix(:,1)))
axis equal
axis on
hold on
for i = 1:size(EleMatrix,1)
    
    %Assigning each stress to an rgb vector depending on its value
    color=Colormap(SigmaRound(i)+max(SigmaRound)+1,:);
    
    r1 = EleMatrix(i,2);
    r2 = EleMatrix(i,3);
    lx = [DefNodeMatrix(r1,2) DefNodeMatrix(r2,2)];
    ly = [DefNodeMatrix(r1,3) DefNodeMatrix(r2,3)];
    lz = [DefNodeMatrix(r1,4) DefNodeMatrix(r2,4)];
    
    lineobj=line(lx,ly,lz,'Color',color);
    if OriMatrix(i,6)==max(OriMatrix(:,6))
        lineobj.LineWidth=2.5;
    else 
        lineobj.LineWidth=1;
    end
    r1 = EleMatrix(i,2);
    r2 = EleMatrix(i,3);
    lx = [NodeMatrix(r1,2) NodeMatrix(r2,2)];
    ly = [NodeMatrix(r1,3) NodeMatrix(r2,3)];
    lz = [NodeMatrix(r1,4) NodeMatrix(r2,4)];
    
    lineobj=line(lx,ly,lz,'Color','k','LineStyle','--'); %structure before deformation
    hold on
end

%Forces Vector
quiver3(NodeMatrix(:,2),NodeMatrix(:,3),NodeMatrix(:,4),ReactionForcesMat(:,2),ReactionForcesMat(:,3),ReactionForcesMat(:,4),'k')
quiver3(DefNodeMatrix(:,2),DefNodeMatrix(:,3),DefNodeMatrix(:,4),LoadMat(:,2),LoadMat(:,3),LoadMat(:,4),'k')
%labels etc
xlabel('x-coordinate [mm]') 
ylabel('y-coordinate [mm]') 
zlabel('z-coordinate [mm]')
title('Stresses of each Truss Element on the deformed structure')
cb=colorbar;
title(cb,'Stress [N/mm^{2}]')
cb.TickLabels=num2cell(round(linspace(-max(SigmaRound),max(SigmaRound),11)));



%Collective Results Tables
figure(3)

sbplt=subplot(1,2,1);
pos1=get(sbplt,'position');
delete(sbplt);
sbplt=subplot(1,2,2);
pos2=get(sbplt,'position');
delete(sbplt);


ResultsTable_1=splitvars(table(EleMatrix,DeformElementPercentage(:,2:end),EleSigma(:,2:end)));
ResultsTable_1.Properties.VariableNames={'Element_ID','Node_1','Node_2','Strain','Normal_Stress_MPa'}

ResultsTable_2=splitvars(table(DeforMatrix,LoadMat(:,2:end),ReactionForcesMat(:,2:end)));
ResultsTable_2.Properties.VariableNames={'Node_ID','Displacement_x','Displacement_y','Displacement_z','ExternalLoad_x','ExternalLoad_y','ExternalLoad_z','ReactionForces_x','ReactionForces_y','ReactionForces_z'}
    
uitable('Data',ResultsTable_1{:,:},'ColumnName',ResultsTable_1.Properties.VariableNames,'Units', 'Normalized', 'Position',pos1);
uitable('Data',ResultsTable_2{:,:},'ColumnName',ResultsTable_2.Properties.VariableNames,'Units', 'Normalized', 'Position',pos2);

hold off