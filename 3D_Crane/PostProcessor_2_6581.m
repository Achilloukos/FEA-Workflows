%Fasoulas Achillefs AEM:6581

%   Data Import
SolCell=DataImport('PreText_2.txt');
NodeMatrix=SolCell{1};
EleMatrix2=SolCell{2};
OriMatrix2=SolCell{3};
BCMatrix2=SolCell{4};
LoadVector2=SolCell{5};

SolCell=DataImport('SolText_2.txt');
ReactionForcesMat2=SolCell{1};
DispMat2=SolCell{2};

LoadMat2(:,1)=NodeMatrix(:,1);
LoadMat2(:,2:4)=[LoadVector2(1:6:163) LoadVector2(2:6:164) LoadVector2(3:6:165)];
LoadMat2(:,5:7)=[LoadVector2(4:6:166) LoadVector2(5:6:167) LoadVector2(6:6:168)];

DefNodeMatrix2=zeros(size(NodeMatrix));
DefNodeMatrix2(:,1)=NodeMatrix(:,1);
DefNodeMatrix2(:,2:4)=NodeMatrix(:,2:4)+10*DispMat2(:,2:4);

%3D-Ploting for Displacements Visualization

figure (1)
plot3(NodeMatrix(:,2),NodeMatrix(:,3),NodeMatrix(:,4),'.')
text(NodeMatrix(:,2),NodeMatrix(:,3),NodeMatrix(:,4),num2str(NodeMatrix(:,1)))


axis equal
axis on
hold on

NodeDisp=cellfun(@norm,num2cell(DispMat2(:,2:4),2));
maxNodeDisp=max(NodeDisp);
NodeDispColoromap=jet(1000+1);


for i = 1:size(EleMatrix2,1)
    
    r1 = EleMatrix2(i,2);
    r2 = EleMatrix2(i,3);
    lx = [DefNodeMatrix2(r1,2) DefNodeMatrix2(r2,2)];
    ly = [DefNodeMatrix2(r1,3) DefNodeMatrix2(r2,3)];
    lz = [DefNodeMatrix2(r1,4) DefNodeMatrix2(r2,4)];
    
    xxx=linspace(lx(1),lx(2),11);
    yyy=linspace(ly(1),ly(2),11);
    zzz=linspace(lz(1),lz(2),11);
    ddd=linspace(NodeDisp(r1),NodeDisp(r2),10);
    
    for j=1:length(xxx)-1
        Percent=1000*ddd(j)/maxNodeDisp + 1;
        LineObj=line(xxx(j:j+1),yyy(j:j+1),zzz(j:j+1),'Color',NodeDispColoromap(round(Percent),:));
        LineObj.LineWidth=2;
    
    
    if OriMatrix2(i,3)==min(OriMatrix2(:,3))
        LineObj.LineWidth=1;
    end
    end 
                
    r1 = EleMatrix2(i,2);
    r2 = EleMatrix2(i,3);
    lx = [NodeMatrix(r1,2) NodeMatrix(r2,2)];
    ly = [NodeMatrix(r1,3) NodeMatrix(r2,3)];
    lz = [NodeMatrix(r1,4) NodeMatrix(r2,4)];
    
    lineobj=line(lx,ly,lz,'Color','k','LineStyle','--');
    hold on
end
quiver3(NodeMatrix(:,2),NodeMatrix(:,3),NodeMatrix(:,4),ReactionForcesMat2(:,2),ReactionForcesMat2(:,3),ReactionForcesMat2(:,4),'k')
quiver3(DefNodeMatrix2(:,2),DefNodeMatrix2(:,3),DefNodeMatrix2(:,4),LoadMat2(:,2),LoadMat2(:,3),LoadMat2(:,4),'k')
xlabel('x-coordinate [mm]') 
ylabel('y-coordinate [mm]') 
zlabel('z-coordinate [mm]')
title('Total Displacement on each part of the beams [mm]')
colormap jet
cb=colorbar;
title(cb,'Displacement [mm]')
cb.TickLabels=num2cell(round(linspace(min(NodeDisp),max(NodeDisp),11)));

hold off

%3D-Ploting for Angular Displacement Visualization

figure(2)
TitleCell={'Elastic Line theta-x','Elastic Line theta-y','Elastic Line theta-z'}

for k=1:3

subplot(2,2,k)
plot3(NodeMatrix(:,2),NodeMatrix(:,3),NodeMatrix(:,4),'.')
text(NodeMatrix(:,2),NodeMatrix(:,3),NodeMatrix(:,4),num2str(NodeMatrix(:,1)))

axis equal
axis on
hold on
xlabel('x-coordinate [mm]') 
ylabel('y-coordinate [mm]') 
zlabel('z-coordinate [mm]')
title(TitleCell{k})
colormap jet
cb=colorbar;
title(cb,'Angular Displacement [deg]')


Arctan=cellfun(@atan,num2cell(DispMat2(:,4+k)));
maxArctan=max(Arctan);
minArctan=min(Arctan);
rangeArctan=max(Arctan)-min(Arctan);
ArctanColoromap=jet(1000+1);

cb.TickLabels=num2cell(round(linspace(min(Arctan),max(Arctan),6),5));
for i = 1:size(EleMatrix2,1)
    
    r1 = EleMatrix2(i,2);
    r2 = EleMatrix2(i,3);
    lx = [DefNodeMatrix2(r1,2) DefNodeMatrix2(r2,2)];
    ly = [DefNodeMatrix2(r1,3) DefNodeMatrix2(r2,3)];
    lz = [DefNodeMatrix2(r1,4) DefNodeMatrix2(r2,4)];
    
    xxx=linspace(lx(1),lx(2),11);
    yyy=linspace(ly(1),ly(2),11);
    zzz=linspace(lz(1),lz(2),11);
    
    
    if EleMatrix2(i)==51 || EleMatrix2(i)==52 || EleMatrix2(i)==53 || EleMatrix2(i)==54
        ddd=zeros(1,10);
    else
        ddd=linspace(Arctan(r1),Arctan(r2),10);
    end
    
    for j=1:length(xxx)-1
        Percent=1000*(ddd(j)-minArctan)/rangeArctan + 1;
        LineObj=line(xxx(j:j+1),yyy(j:j+1),zzz(j:j+1),'Color',ArctanColoromap(round(Percent),:));
        LineObj.LineWidth=2;
    
    
    if OriMatrix2(i,3)==min(OriMatrix2(:,3))
        LineObj.LineWidth=1;
    end
    end 
end
end

%Collective Results Table
sbplt=subplot(2,2,4);
pos1=get(sbplt,'position');
delete(sbplt);

ResultsTable_2=splitvars(table(DispMat2,(ReactionForcesMat2(:,2:4)+LoadMat2(:,2:4))));
ResultsTable_2.Properties.VariableNames={'Node_ID','Displacement_x','Displacement_y','Displacement_z','Angle_x','Angle_y','Angle_z','ReactionForces_x','ReactionForces_y','ReactionForces_z'}
uitable('Data',ResultsTable_2{:,:},'ColumnName',ResultsTable_2.Properties.VariableNames,'Units', 'Normalized', 'Position',pos1);