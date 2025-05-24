% Achillefs Fasoulas 6581

%---------Importing Data--------
Cell = DataImport('PreProcessor_6581.txt');
NodeMatrix = Cell{1};
EleMatrix = Cell{2};
BcMat = Cell{3};
LoadMat =Cell{4};
LoadMator = LoadMat;


Cell2 = DataImport('Solver_6581.txt');
SSCell=Cell2{1};
DisplaceMat=Cell2{2};
Ktheor=round(Cell2{3},3);
Kexp=round(Cell2{4},3);
Deviation=100*round((Kexp-Ktheor)/Ktheor,3);


idx = listdlg('ListString' , {'<HTML>Von Mises &sigma;<sub>x</sub> &sigma;<sub>y</sub> &tau;<sub>xy</sub>'...
               ,'<HTML>&sigma;<sub>1</sub> &sigma;<sub>2</sub>','<HTML>&epsilon;<sub>x</sub> &epsilon;<sub>y</sub> &gamma;<sub>xy</sub>','Diplacements'},...
               'PromptString' ,{'Select which diagrams to plot' },'SelectionMode','single');

epsilon=SSCell(:,2:4);

maxEpsilonx=max(epsilon(:,1));
minEpsilonx=min(epsilon(:,1));

maxEpsilony=max(epsilon(:,2));
minEpsilony=min(epsilon(:,2));

maxGammaxy=max(epsilon(:,3));
minGammaxy=min(epsilon(:,3));

sigma=SSCell(:,5:7);

maxSigmax=max(sigma(:,1));
minSigmax=min(sigma(:,1));

maxSigmay=max(sigma(:,2));
minSigmay=min(sigma(:,2));

maxTafxy=max(sigma(:,3));
minTafxy=min(sigma(:,3));

principal=SSCell(:,8:10);

maxSigma1=max(principal(:,1));
minSigma1=min(principal(:,1));

maxSigma2=max(principal(:,2));
minSigma2=min(principal(:,2));


vonMises=SSCell(:,11);

maxvM=max(vonMises);
minvM=min(vonMises);

%Plot charts for normal, shear and vonMises stresses
if ismember(1,idx)
figure(1)

annotation('textbox' , [0.02 .1 .1 .2], 'String',...
        {['K_t_h_e_o_r = ',  num2str(Ktheor)] , ['K_e_x_p = ' num2str(Kexp)] ,...
        ['Deviation = ' num2str(Deviation) '%']},'EdgeColor' , 'none')

a1=subplot(2,2,1);
axis on
hold on
xlabel('x-coordinate [mm]') 
ylabel('y-coordinate [mm]') 
title('von Mises Stress')
colormap jet
cb=colorbar;
title(cb,'Stress [N/mm^2]')
cb.TickLabels=num2cell(round(linspace(minvM,maxvM,6),5));

a2=subplot(2,2,2);
axis on
hold on
xlabel('x-coordinate [mm]') 
ylabel('y-coordinate [mm]') 
title('\sigma_x')
colormap jet
cb=colorbar;
title(cb,'Normal Stress [N/mm^2]')
cb.TickLabels=num2cell(round(linspace(minSigmax,maxSigmax,6),5));

a3=subplot(2,2,3);
axis on
hold on
xlabel('x-coordinate [mm]') 
ylabel('y-coordinate [mm]') 
title('\sigma_y')
colormap jet
cb=colorbar;
title(cb,'Normal Stress [N/mm^2]')
cb.TickLabels=num2cell(round(linspace(minSigmay,maxSigmay,6),5));


a4=subplot(2,2,4);
axis on
hold on
xlabel('x-coordinate [mm]') 
ylabel('y-coordinate [mm]') 
title('\tau_x_y')
colormap jet
cb=colorbar;
title(cb,'Shear Stress [N/mm^2]')
cb.TickLabels=num2cell(round(linspace(minTafxy,maxTafxy,6),5));

end

%Plot charts for Principal stresses
if ismember(2,idx)
figure(2)

annotation('textbox' , [0.02 .1 .1 .2], 'String',...
        {['K_t_h_e_o_r = ',  num2str(Ktheor)] , ['K_e_x_p = ' num2str(Kexp)] ,...
        ['Deviation = ' num2str(Deviation) '%']},'EdgeColor' , 'none')
a5=axes;
axis on
hold on
xlabel('x-coordinate [mm]') 
ylabel('y-coordinate [mm]') 
title('\sigma_1')
colormap jet
cb=colorbar;
title(cb,'Principal Stress [N/mm^2]')
cb.TickLabels=num2cell(round(linspace(minSigma1,maxSigma1,11),5));

figure(3)
a6=axes;
axis on
hold on
xlabel('x-coordinate [mm]') 
ylabel('y-coordinate [mm]') 
title('\sigma_2')
colormap jet
cb=colorbar;
title(cb,'Principal Stress [N/mm^2]')
cb.TickLabels=num2cell(round(linspace(minSigma2,maxSigma2,11),5));

end

%Plot charts for Strains/Deformations
if ismember(3,idx)
figure(3)

annotation('textbox' , [0.02 .1 .1 .2], 'String',...
        {['K_t_h_e_o_r = ',  num2str(Ktheor)] , ['K_e_x_p = ' num2str(Kexp)] ,...
        ['Deviation = ' num2str(Deviation) '%']},'EdgeColor' , 'none')

a7=subplot(2,2,1);
axis on
hold on
xlabel('x-coordinate [mm]') 
ylabel('y-coordinate [mm]') 
title('\epsilon_x')
colormap jet
cb=colorbar;
title(cb,'Strain [-]')
cb.TickLabels=num2cell(round(linspace(minEpsilonx,maxEpsilonx,6),5));

a8=subplot(2,2,2);
axis on
hold on
xlabel('x-coordinate [mm]') 
ylabel('y-coordinate [mm]') 
title('\epsilon_y')
colormap jet
cb=colorbar;
title(cb,'Strain [-]')
cb.TickLabels=num2cell(round(linspace(minEpsilony,maxEpsilony,6),5));

a9=subplot(2,2,3);
axis on
hold on
xlabel('x-coordinate [mm]') 
ylabel('y-coordinate [mm]') 
title('\gamma_x_y')
colormap jet
cb=colorbar;
title(cb,'Angular Strain [-]')
cb.TickLabels=num2cell(round(linspace(minGammaxy,maxGammaxy,6),5));
end

%Visualizations for everything besides Displacements
Colormap=jet(1001);
for i = 1: size(EleMatrix,1)
        N1 = NodeMatrix(EleMatrix(i,2),2:3);
        N2 = NodeMatrix(EleMatrix(i,3),2:3);
        N3 = NodeMatrix(EleMatrix(i,4),2:3);
        
            
%       text(axes,(N1(1)+N2(1)+N3(1))/3 , (N1(2)+N2(2)+N3(2))/3 , num2str(EleMatrix(i,1)),...
%           'Color' , 'r' , 'FontSize',15 )
        x=[N1(1) N2(1) N3(1)]  ;    
        y=[N1(2) N2(2) N3(2)]  ;
        
        if ismember(1,idx)
        %von Mises
        vM=vonMises(i);
        vMperc=(vM-minvM)/(maxvM-minvM);
        
        color=Colormap(round(1000*vMperc)+1,:);
        patch(a1,x,y,color);
        
        %Sigma_x
        sigma_x=sigma(i,1);
        sigma_x_perc=(sigma_x-minSigmax)/(maxSigmax-minSigmax);
        
        color=Colormap(round(1000*sigma_x_perc)+1,:);
        patch(a2,x,y,color);
        
        %Sigma_y
        sigma_y=sigma(i,2);
        sigma_y_perc=(sigma_y-minSigmay)/(maxSigmay-minSigmay);
        
        color=Colormap(round(1000*sigma_y_perc)+1,:);
        patch(a3,x,y,color);
        
        %Taf_xy
        taf_xy=sigma(i,3);
        taf_xy_perc=(taf_xy-minTafxy)/(maxTafxy-minTafxy);
        
        color=Colormap(round(1000*taf_xy_perc)+1,:);
        patch(a4,x,y,color);
        end
        
        if ismember(2,idx)
        %Principal_1
        sigma_1=principal(i,1);
        sigma_1_perc=(sigma_1-minSigma1)/(maxSigma1-minSigma1);
        
        color=Colormap(round(1000*sigma_1_perc)+1,:);
        patch(a5,x,y,color);
        
%         quiver(a5,mean(x),mean(y),principal(i,1)*cosd(principal(i,3)),principal(i,1)*sind(principal(i,3)),'Color','k')
        
        %Principal_2
        sigma_2=principal(i,2);
        sigma_2_perc=(sigma_2-minSigma2)/(maxSigma2-minSigma2);
        
        color=Colormap(round(1000*sigma_2_perc)+1,:);
        patch(a6,x,y,color);
        
%         quiver(a6,mean(x),mean(y),principal(i,2)*cosd(principal(i,3)-90),principal(i,2)*sind(principal(i,3)-90),'Color','k')
        end
        
        if ismember(3,idx)
    %Strains
        
        %Epsilon_x
        epsilon_x=epsilon(i,1);
        epsilon_x_perc=(epsilon_x-minEpsilonx)/(maxEpsilonx-minEpsilonx);
        
        color=Colormap(round(1000*epsilon_x_perc)+1,:);
        patch(a7,x,y,color);
        
        %Epsilon_y
        epsilon_y=epsilon(i,2);
        epsilon_y_perc=(epsilon_y-minEpsilony)/(maxEpsilony-minEpsilony);
        
        color=Colormap(round(1000*epsilon_y_perc)+1,:);
        patch(a8,x,y,color);
        
        %Gamma_xy
        gamma_xy=epsilon(i,3);
        gamma_xy_perc=(gamma_xy-minGammaxy)/(maxGammaxy-minGammaxy);
        
        color=Colormap(round(1000*gamma_xy_perc)+1,:);
        patch(a9,x,y,color);
        end
       

end
% hold on

if ismember(2,idx)
for i = 1: size(EleMatrix,1)
    N1 = NodeMatrix(EleMatrix(i,2),2:3);
    N2 = NodeMatrix(EleMatrix(i,3),2:3);
    N3 = NodeMatrix(EleMatrix(i,4),2:3);
    %Polygon Coordinates
    x=[N1(1) N2(1) N3(1)];
    y=[N1(2) N2(2) N3(2)];
    
    quiver(a5,mean(x),mean(y),principal(i,1)*cos(principal(i,3)),principal(i,1)*sin(principal(i,3)),...
             3 ,'Color' , 'k');
    
        
    quiver(a6,mean(x),mean(y),principal(i,2)*cos(principal(i,3)-pi/2),principal(i,2)*sin(principal(i,3)-pi/2),...
             3 ,'Color' , 'k');
    
end
hold on
end

if ismember(4,idx)
TotalDispl=sqrt(DisplaceMat(:,2).^2+DisplaceMat(:,3).^2);
PlotFieldonMesh(NodeMatrix(:,2:3)+250*DisplaceMat(:,2:3),EleMatrix(:,2:4),TotalDispl);

annotation('textbox' , [0.02 .1 .1 .2], 'String',...
        {['K_t_h_e_o_r = ',  num2str(Ktheor)] , ['K_e_x_p = ' num2str(Kexp)] ,...
        ['Deviation = ' num2str(Deviation) '%']},'EdgeColor' , 'none')

axis equal;
colormap jet
cb=colorbar;
title(cb,'[mm]')
end


