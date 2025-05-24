%Fasoulas Achillefs AEM:6581
%Protypo provlima %

%-------PreProcessor-------%
testNodeMat = [1 0 0 ;
               2 40 0 ;
               3 40 30 ;
               4 0 30 ];

testelem = [ 1 1 2;
            2 3 2 ;
            3 1 3;
            4 4 3];
testorient = [1 40 1 0 1;
              2 30 0 -1 1;
              3 50 0.8 0.6 1;
              4 40 1 0 1];
          
%-------Solver-------%
w = 29.5*10^6/600;
k1 = w*15*[1 0 -1 0 ;
      0 0 0 0 ;
      -1 0 1 0;
      0 0 0 0 ];
k2 = w*[0 0 0 0 ;
      0 20 0 -20;
      0 0 0 0;
      0 -20 0 20];
k3 = w*[7.68 5.76 -7.68 -5.76;
        5.76 4.32 -5.76 -4.32 ;
        -7.68 -5.76 7.68 5.76 ;
        -5.76 -4.32 5.76 4.32];    
k4 = w*[15 0 -15 0;
        0 0 0 0 ;
        -15 0 15 0;
        0 0 0 0 ];   
testcell = {k1, k2 ,k3 ,k4};
StiffnessMatrix = zeros(8); 
for i = 1:4
    IDi = testelem(i,1);
    SMi = testcell{i};
    Nodes = testelem(i,2:3);
    DoFi = [2*Nodes(1)-1:2*Nodes(1) 2*Nodes(2)-1:2*Nodes(2)];
    for j = 1:4
        for k = 1:4
           StiffnessMatrix(DoFi(j),DoFi(k)) =  StiffnessMatrix(DoFi(j),DoFi(k))+ SMi(j,k);
        end
    end
end
testLoad = zeros(2*size(testNodeMat,1),1);
testLoad(6) = - 25000;
testLoad(3) = 20000;
LoadMat = [0 0 ; 20000 0 ; 0 -25000; 0 0];
testBcs = zeros(2*size(testNodeMat,1),1);
testBcs([1 2 4 7 8])= 1;
Eliminate =find(testBcs);
testLoad(Eliminate) = [];
StiffnessMatrixR = StiffnessMatrix;
StiffnessMatrixR(Eliminate,:)=[];
StiffnessMatrixR(:,Eliminate)=[];
displacement = myinv(StiffnessMatrixR)*testLoad;
DispF = [ 0 0; displacement(1) 0; displacement(2:3)';0 0];
for i =1:4
    Nodes = testelem(i,2:3);
    DoF = [2*Nodes(1)-1:2*Nodes(1) 2*Nodes(2)-1:2*Nodes(2)];
    cosines = testorient(i,3:4);
    M = [cosines 0 0 ; 
        0 0 cosines];
    sigma(i) = (29.5*10^6/testorient(i,2)) * ([-1 1]*M*DispF(DoF)');
end

%Visualization-PostProcessor
hold on
plot(testNodeMat(:,2),testNodeMat(:,3),'o')
text(testNodeMat(:,2)+1,testNodeMat(:,3)-1,num2str(testNodeMat(:,1)))
sigmap = round(sigma);
colormy = jet(1001);
for i=1:4
    percent = 1000*(sigmap(i)-min(sigmap))/(max(sigmap)-min(sigmap))+1;
    r1 = testelem(i,2);
    r2 = testelem(i,3);
    lx = [testNodeMat(r1,2) testNodeMat(r2,2)];
    ly = [testNodeMat(r1,3) testNodeMat(r2,3)];
    l=line(lx,ly,'Color',colormy(round(percent),:));
    l.LineWidth = 2;
end
colormap jet
cb=colorbar;
cb.TickLabels = num2cell(linspace(min(sigmap),max(sigmap),11));
quiver(testNodeMat(:,2),testNodeMat(:,3),LoadMat(:,1),LoadMat(:,2),'LineWidth',1.5,'Color','k')
hold off
disp("Element Stress: ")
disp([[1 2 3 4]' sigma'])
disp("Node Displacement: ")
disp([[1 2 3 4]' DispF])




    