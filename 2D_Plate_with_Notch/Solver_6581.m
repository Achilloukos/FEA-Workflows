% Achillefs Fasoulas 6581

[s, H, W , r, t] = PreprocessorData();
Cell = DataImport('PreProcessor_6581.txt');
NodeMatrix = Cell{1};
EleMatrix = Cell{2};
BCMat = Cell{3};
LoadMat =Cell{4};
LoadMator = LoadMat;

F = 1000; %[N] 
M = 10000; %N*mm

LoadMatInit=LoadMat;
E = 21000; %N/mm^2
nu = 0.28; %Poisson's ratio


%Degrees of freedom indexing
keys = 1:size(NodeMatrix,1);
values = cell(1,size(NodeMatrix,1));
for i =1:size(NodeMatrix,1)
    values{i} = 2*i-1:2*i;
end
DoFmap = containers.Map(keys,values);
clear keys values

%Creating Stiffness Matrix for each Element
thic = 4; %mm
D = E/(1-nu^2) * [1 nu 0 ; nu 1 0 ; 0 0 (1-nu)*0.5];


StiffnessCell = cell(size(EleMatrix,1) ,2);
for i = 1:size(EleMatrix,1)
    N1 = NodeMatrix(EleMatrix(i,2),2:3);
    N2 = NodeMatrix(EleMatrix(i,3),2:3);
    N3 = NodeMatrix(EleMatrix(i,4),2:3);
    A = Tarea(N1,N2,N3);
    [x ,y] = difs(N1,N2,N3);
    
    dJ = x(1,3)*y(2,3)-y(1,3)*x(2,3);
    
    B = (1/dJ)*[ y(2,3) 0 y(3,1) 0 y(1,2) 0;
               0 x(3,2) 0 x(1,3) 0 x(2,1);
               x(3,2) y(2,3) x(1,3) y(3,1) x(2,1) y(1,2)];
      
    StiffnessCell(i,1:2) = {i , thic*A*B'*D*B};   
    
end
clear N1 N2 N3 A B dJ x y 


%Summing Global Stiffness Matrix
StiffnessMatrix = zeros(2*size(NodeMatrix,1));

for i = 1:size(StiffnessCell,1)
    IDi = StiffnessCell{i,1};
    SMi = StiffnessCell{i,2};
    Nodes = EleMatrix(IDi,2:4);
    DoFi = [DoFmap(Nodes(1)) DoFmap(Nodes(2)) DoFmap(Nodes(3))];
    for j = 1:size(SMi,1)
        for k = 1:size(SMi,2)
           StiffnessMatrix(DoFi(j),DoFi(k)) =  StiffnessMatrix(DoFi(j),DoFi(k))+ SMi(j,k);
        end
    end
end
clear i j k IDi SMi Nodes DoFi

%Boundary condition Reshape to vector
BCMat = BCMat(:,2:3)';
BCMat = int32(BCMat(:));

%Load Matrix Reshape to vector
LoadMat = LoadMat(:,2:3)';
LoadMat = LoadMat(:);

%find DoF to eliminate
DoF =(1:2*size(NodeMatrix,1))';
Eliminate = find(BCMat);

%Eliminate all appropriate DoF
DoF(Eliminate) = [];
LoadMat(Eliminate) = [];
StiffnessMatrixReduced = StiffnessMatrix;
StiffnessMatrixReduced(Eliminate,:) = [];
StiffnessMatrixReduced(:,Eliminate) = [];

%Solving Linear System
DisplaceVec=StiffnessMatrixReduced\LoadMat;

%Organizing Displacements
DisplaceMatFull = zeros(2*size(NodeMatrix,1),1);
DisplaceMatFull(DoF) = DisplaceVec;
DisplaceMat = zeros(size(NodeMatrix));
DisplaceMat(:,1) = NodeMatrix(:,1);
DisplaceMat(:,2:3) = reshape(DisplaceMatFull,size(DisplaceMat(:,2:3),2),size(DisplaceMat(:,2:3),1))';

%Stresses && Strains Calculation
% SSCell=cell(size(EleMatrix,1),5);
SSCell=zeros(size(EleMatrix,1),10);


for i = 1:size(EleMatrix,1)
   q = [DisplaceMatFull(DoFmap(EleMatrix(i,2)))' DisplaceMatFull(DoFmap(EleMatrix(i,3)))' DisplaceMatFull(DoFmap(EleMatrix(i,4)))'];
   
   N1 = NodeMatrix(EleMatrix(i,2),2:3);
    N2 = NodeMatrix(EleMatrix(i,3),2:3);
    N3 = NodeMatrix(EleMatrix(i,4),2:3);
   [x ,y] = difs(N1,N2,N3);
    
    dJ = x(1,3)*y(2,3)-y(1,3)*x(2,3);
    
    B = (1/dJ)*[ y(2,3) 0 y(3,1) 0 y(1,2) 0;
               0 x(3,2) 0 x(1,3) 0 x(2,1);
               x(3,2) y(2,3) x(1,3) y(3,1) x(2,1) y(1,2)]; 
   
   
   epsilon = B*q';
   sigma = D*epsilon;
   
   R = sqrt(0.5*(sigma(1)-sigma(2))^2 + sigma(3)^2);
   avg = (sigma(1)+sigma(2))/2;
   theta=0.5*atan(2*sigma(3)/(sigma(1)-sigma(2)));
   principal = [avg+R avg-R theta]';

   
   vonMises = sqrt(sigma(2)^2 - sigma(2)*sigma(1) + sigma(1)^2 + 3*sigma(3)^2);
   
   SSCell(i,1:11) = [i epsilon' sigma' principal' vonMises];
   
end



if min(LoadMatInit(:,2:3))==0
     %Shape Coefficient for Tension
    Sigma_Onomastiki=F/((H-2*t)*thic);
    Sigmaximum=max(SSCell(:,5));
    Theoretical_Shape_Coefficient=(0.78+2.243*sqrt(t/r))*(0.993+0.18*(2*t/H)-1.06*(2*t/H)^2+1.71*(2*t/H)^3)*(1-2*t/H);
    Calculated_Shape_Coefficient=Sigmaximum/Sigma_Onomastiki;

    
else
   
    %Shape Coefficient for Bending
    Sigma_Onomastiki=6*M/(thic*(H-2*t)^2);
    Sigmaximum=max(SSCell(:,5));
    Theoretical_Shape_Coefficient=1.113+1.957*sqrt(t/r)+(-2.579-4.017*sqrt(t/r)-0.013*(t/r))*(2*t/H)+(4.1+3.922*sqrt(t/r)+0.083*(t/r))*(2*t/H)^2+(-1.528-1.893*sqrt(t/r)-0.066*(t/r))*(2*t/H)^3;
    Calculated_Shape_Coefficient=Sigmaximum/Sigma_Onomastiki;
end

clear i N1 N2 N3 x y dJ B epsilon sigma R avg principal theta


DataExport('Solver_6581.txt' ,SSCell , DisplaceMat,Theoretical_Shape_Coefficient,Calculated_Shape_Coefficient)
% save('Solver.mat','SSCell','DisplaceMat')

function [xmat, ymat] = difs (N1,N2,N3)

mat = [N1' N2' N3'];
xmat = repelem(mat(1,:),3,1)' - repelem(mat(1,:),3,1);
ymat = repelem(mat(2,:),3,1)' - repelem(mat(2,:),3,1);

    
end
    



function [A] = Tarea(N1,N2,N3)
    %This function calculates triangle area based on coordinates
    
    A = 0.5*abs(N1(1)*N2(2) + N2(1)*N3(2) + N3(1)*N1(2)...
           -N1(2)*N2(1) - N2(2)*N3(1) - N3(2)*N1(1));
end