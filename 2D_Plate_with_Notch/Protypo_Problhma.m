%% PREPROCCESSOR
% Achillefs Fasoulas 6581

[idx ,tf] = listdlg('ListString', {'Tension' ,'Bending'}, 'SelectionMode' , 'single','ListSize',[150,100],...
               'PromptString', 'Select type of load');

if ~tf
    error('No Load Selected')
end


v1 = 0:2:80;
v2 = 0:2:40;
[x,y] = meshgrid(v1,v2);
x = x';
y = y';
NodeMatrix=zeros(length(x(:)),3);
NodeMatrix(:,1) = 1:size(NodeMatrix,1);
NodeMatrix(:,2) = x(:);
NodeMatrix(:,3) = y(:);

%Creating Element Matrix
Elematrix = zeros(2*(length(v1)-1)*(length(v2)-1) , 4);
Elematrix(:,1) = 1:size(Elematrix,1);
c = 1;



for i = 1:2:size(Elematrix,1)
    if mod(c,length(v1)) == 0
        c=c+1;
    end
    Elematrix(i,2:4) = [c ,c+1 ,c+length(v1)];
    Elematrix(i+1,2:4) = [c+1 , c+length(v1)+1, c+length(v1)];
    c = c+1;
end


BcMat = zeros(size(NodeMatrix));
BcMat(:,1) = 1:size(BcMat,1);
Px = find(NodeMatrix(:,2) == 0);
BcMat(Px(1),3) = 1;
BcMat(Px,2) = 1;
clear Px

LoadMat = zeros(size(NodeMatrix,1),3);
RightID = find(NodeMatrix(:,2) == v1(end));
[~ ,RightIDy] = sort(NodeMatrix(RightID,3));
RightID =RightID(RightIDy);
Deltay = abs(NodeMatrix(RightID(2:end),3) - NodeMatrix(RightID(1:end-1),3));

if idx==1
F = 1000; %[N] 
Fdist = F*(Deltay/v2(end));
for i = 1:size(RightID,1)-1
    LoadMat(RightID(i),2) = LoadMat(RightID(i),2) + Fdist(i)/2;
    LoadMat(RightID(i+1),2) = LoadMat(RightID(i+1),2) + Fdist(i)/2;    
end
end

if idx==2
 M = 10000; %N*mm
   for i = 1:size(RightID,1)-1
       
       f1 = (12*M)/(v2(end)^3)* (NodeMatrix(RightID(i),3)-v2(end)/2);
       f2 = (12*M)/(v2(end)^3)* (NodeMatrix(RightID(i+1),3)-v2(end)/2);
       F = (f1+f2)*Deltay(i)*0.5;
       LoadMat(RightID(i),2) = LoadMat(RightID(i),2) + F/2;
       LoadMat(RightID(i+1),2) = LoadMat(RightID(i+1),2) + F/2;
        
   end
end

plot(NodeMatrix(:,2),NodeMatrix(:,3),'.')

quiver(NodeMatrix(:,2),NodeMatrix(:,3),LoadMat(:,2),LoadMat(:,3),'r')
axis equal
% text(NodeMatrix(:,2),NodeMatrix(:,3),num2str(NodeMatrix(:,1)))
for i = 1: size(Elematrix,1)
        N1 = NodeMatrix(Elematrix(i,2),2:3);
        N2 = NodeMatrix(Elematrix(i,3),2:3);
        N3 = NodeMatrix(Elematrix(i,4),2:3);
        
       line([N1(1) N2(1);N2(1) N3(1);N3(1) N1(1)],...
                [N1(2) N2(2);N2(2) N3(2);N3(2) N1(2)],...
                'Color','b')
             
%       text(axes,(N1(1)+N2(1)+N3(1))/3 , (N1(2)+N2(2)+N3(2))/3 , num2str(Elematrix(i,1)),...
%           'Color' , 'r' , 'FontSize',15 )
        
end


%% SOLVER
E = 210000;
nu = 0.28;



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

StiffnessCell = cell(size(Elematrix,1) ,2);
for i = 1:size(Elematrix,1)
    N1 = NodeMatrix(Elematrix(i,2),2:3);
    N2 = NodeMatrix(Elematrix(i,3),2:3);
    N3 = NodeMatrix(Elematrix(i,4),2:3);
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
    Nodes = Elematrix(IDi,2:4);
    DoFi = [DoFmap(Nodes(1)) DoFmap(Nodes(2)) DoFmap(Nodes(3))];
    for j = 1:size(SMi,1)
        for k = 1:size(SMi,2)
           StiffnessMatrix(DoFi(j),DoFi(k)) =  StiffnessMatrix(DoFi(j),DoFi(k))+ SMi(j,k);
        end
    end
end
clear i j k IDi SMi Nodes DoFi


%Boundary condition Reshape to vector
BcMat = BcMat(:,2:3)';
BcMat = int32(BcMat(:));


%Load Matrix Reshape to vector
LoadMat = LoadMat(:,2:3)';
LoadMat = LoadMat(:);

%find DoF to eliminate
DoF =(1:2*size(NodeMatrix,1))';
Eliminate = find(BcMat);

%Eliminate appropriate DoF
DoF(Eliminate) = [];
LoadMat(Eliminate) = [];
StiffnessMatrixReduced = StiffnessMatrix;
StiffnessMatrixReduced(Eliminate,:) = [];
StiffnessMatrixReduced(:,Eliminate) = [];

%Solve Linear system
displacement = StiffnessMatrixReduced\LoadMat;

%Organise displacements
displacementsfull = zeros(2*size(NodeMatrix,1),1);
displacementsfull(DoF) = displacement;
displacementMat = zeros(size(NodeMatrix,1),2);
displacementMat = reshape(displacementsfull,size(displacementMat(:,1:2),2),size(displacementMat(:,1:2),1))';

% Computing stress and strain
ResultCell = zeros(size(Elematrix,1) ,13);

for i = 1:size(Elematrix,1)
    N1 = NodeMatrix(Elematrix(i,2),2:3);
    N2 = NodeMatrix(Elematrix(i,3),2:3);
    N3 = NodeMatrix(Elematrix(i,4),2:3); 
    
   q = [displacementsfull(DoFmap(Elematrix(i,2)))' displacementsfull(DoFmap(Elematrix(i,3)))' displacementsfull(DoFmap(Elematrix(i,4)))']';
   [x ,y] = difs(N1,N2,N3);
    
    dJ = x(1,3)*y(2,3)-y(1,3)*x(2,3);
    
    B = (1/dJ)*[ y(2,3) 0 y(3,1) 0 y(1,2) 0;
               0 x(3,2) 0 x(1,3) 0 x(2,1);
               x(3,2) y(2,3) x(1,3) y(3,1) x(2,1) y(1,2)]; 
   
   
   epsilon = B*q;
   sigma = D*epsilon;
   R = sqrt(0.5*(sigma(1)-sigma(2))^2 + sigma(3)^2);
   avg = (sigma(1)+sigma(2))/2;
   theta = 0.5*atan((2*sigma(3))/(sigma(1)-sigma(2)));
   principal = [avg+R avg-R theta]';
   
   vonMises = sqrt(sigma(2)^2 - sigma(2)*sigma(1) + sigma(1)^2 + 3*sigma(3)^2);
   
   
   ResultCell(i,1:11) = [i, epsilon', sigma', principal', vonMises];
   
end




%% POST PROCCESSOR
epsilon = ResultCell(:,2:4);
sigma=ResultCell(:,5:7);
principal=ResultCell(:,8:10);
vonMises=ResultCell(:,11);
Colormap =jet(1001);

idx = listdlg('ListString' , {'Von Mises Stress', '<HTML>&sigma;<sub>x</sub>' , '<HTML>&sigma;<sub>y</sub>' , '<HTML>&tau;<sub>xy</sub>'...
               ,'<HTML>&sigma;<sub>1</sub>','<HTML>&sigma;<sub>2</sub>','<HTML>&epsilon;<sub>x</sub>','<HTML>&epsilon;<sub>y</sub>','<HTML>&gamma;<sub>xy</sub>'},...
               'PromptString' ,{'Select which diagrams to plot' ,'(you can select multiple) '});

          
%Von Mises Subplot
if ismember(1,idx)          
figure(1);
a1 = axes;
maxvM=max(vonMises);
minvM=min(vonMises);
axis on
hold on
xlabel('x-coordinate [mm]') 
ylabel('y-coordinate [mm]') 
title('von Mises Stress')
colormap jet
cb=colorbar;
title(cb,'[N/mm^2]')
cb.TickLabels=num2cell(round(linspace(minvM,maxvM,11),5));
end

% Sigma x subplot
if ismember(2,idx)
figure(2);
a2 = axes;
maxSigmax=max(sigma(:,1));
minSigmax=min(sigma(:,1));
axis on
hold on
xlabel('x [mm]') 
ylabel('y [mm]')
title('\sigma_x')
colormap jet
cb=colorbar;
title(cb,'[N/mm^2]')
cb.TickLabels=num2cell(round(linspace(minSigmax,maxSigmax,11),5));
end

% Sigma_y
if ismember(3,idx)
figure(3);
a3 = axes;
maxSigmay=max(sigma(:,2));
minSigmay=min(sigma(:,2));
axis on
hold on
xlabel('x [mm]') 
ylabel('y [mm]') 
title('\sigma_y')
colormap jet
cb=colorbar;
title(cb,'[N/mm^2]')
cb.TickLabels=num2cell(round(linspace(minSigmay,maxSigmay,11),5));
end

%Taf xy subplot
if ismember(4,idx)
figure(4);
a4 = axes;
maxTafxy=max(sigma(:,3));
minTafxy=min(sigma(:,3));
axis on
hold on
xlabel('x [mm]') 
ylabel('y [mm]') 
title('\tau_x_y')
colormap jet
cb=colorbar;
title(cb,'[N/mm^2]')
cb.TickLabels=num2cell(round(linspace(minTafxy,maxTafxy,11),5));
end

% Sigma_1
if ismember(5,idx)
figure(5);
a5 = axes;
maxSigma1=max(principal(:,1));
minSigma1=min(principal(:,1));
axis on
hold on
xlabel('x [mm]') 
ylabel('y [mm]') 
title('\sigma_1')
colormap jet
cb=colorbar;
title(cb,'[N/mm^2]')
cb.TickLabels=num2cell(round(linspace(minSigma1,maxSigma1,11),5));
end 

% Sigma2
if ismember(6,idx)
figure(6);
a6 =  axes;
maxSigma2 = max(principal(:,2));
minSigma2 = min(principal(:,2));
axis on
hold on
xlabel('x [mm]') 
ylabel('y [mm]') 
title('\sigma_2')
colormap jet
cb=colorbar;
title(cb,'[N/mm^2]')
cb.TickLabels=num2cell(round(linspace(minSigma2,maxSigma2,11),5));
end

% epsilonx
if ismember(7,idx)
figure(7);
a7 = axes;
maxEpsilonx=max(epsilon(:,1));
minEpsilonx=min(epsilon(:,1));
axis on
hold on
xlabel('x [mm]') 
ylabel('y [mm]') 
title('\epsilon_x')
colormap jet
cb=colorbar;
% title(cb,'[N/mm^2]')
cb.TickLabels=num2cell(round(linspace(minEpsilonx,maxEpsilonx,11),5));
end

%epsilony
if ismember(8,idx)
figure(8)
a8 =axes;
maxEpsilony=max(epsilon(:,2));
minEpsilony=min(epsilon(:,2));
axis on
hold on
xlabel('x [mm]') 
ylabel('y [mm]') 
title('\epsilon_y')
colormap jet
cb=colorbar;
% title(cb,'[N/mm^2]')
cb.TickLabels=num2cell(round(linspace(minEpsilony,maxEpsilony,11),5));
end

% gamma xy
if ismember(9,idx)
figure(9);
a9 = axes;
maxGamaxy=max(epsilon(:,3));
minGamaxy=min(epsilon(:,3));
axis on
hold on
xlabel('x [mm]') 
ylabel('y [mm]') 
title('\gamma_x_y')
colormap jet
cb=colorbar;
% title(cb,'[N/mm^2]')
cb.TickLabels=num2cell(round(linspace(minGamaxy,maxGamaxy,11),5));
end


for i = 1: size(Elematrix,1)
    N1 = NodeMatrix(Elematrix(i,2),2:3);
    N2 = NodeMatrix(Elematrix(i,3),2:3);
    N3 = NodeMatrix(Elematrix(i,4),2:3);


%    text(axes,(N1(1)+N2(1)+N3(1))/3 , (N1(2)+N2(2)+N3(2))/3 , num2str(Elematrix(i,1)),...
%               'Color' , 'r' , 'FontSize',15 )
    
    %Polygon Coordinates
    x=[N1(1) N2(1) N3(1)];
    y=[N1(2) N2(2) N3(2)];
    
    if ismember(1,idx)
    %Von Mises Stress
    vM=vonMises(i);
    vMperc=(vM-minvM)/(maxvM-minvM);
    color=Colormap(round(1000*vMperc)+1,:);
    patch(a1,x,y,color);
    end
    
    if ismember(2,idx)
    %Sigma_x
    sigma_x=sigma(i,1);
    sigma_x_perc=(sigma_x-minSigmax)/(maxSigmax-minSigmax);
    color=Colormap(round(1000*sigma_x_perc)+1,:);
    patch(a2,x,y,color);
    end

    if ismember(3,idx)
    %Sigma_y
    sigma_y=sigma(i,2);
    sigma_y_perc=(sigma_y-minSigmay)/(maxSigmay-minSigmay);
    color=Colormap(round(1000*sigma_y_perc)+1,:);
    patch(a3,x,y,color);
    end

    if ismember(4,idx)
    %Taf xy
    tafxy = sigma(i,3);
    tafxy_perc=(tafxy-minTafxy)/(maxTafxy-minTafxy);
    color=Colormap(round(1000*tafxy_perc)+1,:);
    patch(a4,x,y,color);
    end
    
    if ismember(5,idx)
    %Sigma_1
    sigma_1 = principal(i,1);
    sigma_1_perc = (sigma_1-minSigma1)/(maxSigma1-minSigma1);
    color = Colormap(round(1000*sigma_1_perc)+1,:);
    patch(a5,x,y,color);    
    %Arrows
    quiver(a5,mean(x),mean(y),principal(i,1)*cos(principal(i,3)),principal(i,1)*sin(principal(i,3)),...
             3 ,'Color' , 'k');
    end
    
     if ismember(6,idx)   
    %Sigma_2
    sigma_2 = principal(i,2);
    sigma_2_perc = (sigma_2-minSigma2)/(maxSigma2-minSigma2);
    color = Colormap(round(1000*sigma_2_perc)+1,:);
    patch(a6,x,y,color);
    %Arrows
    quiver(a6,mean(x),mean(y),principal(i,2)*cos(principal(i,3)-pi/2),principal(i,2)*sin(principal(i,3)-pi/2),...;
             3 ,'Color' , 'k');
     end
     
    if ismember(7,idx)
    %epsilon x
    epsilon_x= epsilon(i,1);
    epsilon_x_perc = (epsilon_x-minEpsilonx)/(maxEpsilonx-minEpsilonx);
    color = Colormap(round(1000*epsilon_x_perc)+1,:);
    patch(a7,x,y,color);
    end
    
    
    if ismember(8,idx)
    %epsilon y
    epsilon_y= epsilon(i,2);
    epsilon_y_perc = (epsilon_y-minEpsilony)/(maxEpsilony-minEpsilony);
    color = Colormap(round(1000*epsilon_y_perc)+1,:);
    patch(a8,x,y,color);
    end
    
    if ismember(9,idx)
    % gamma xy 
    gamma_xy= epsilon(i,3);
    gamma_xy_perc = (gamma_xy-minGamaxy)/(maxGamaxy-minGamaxy);
    color = Colormap(round(1000*gamma_xy_perc)+1,:);
    patch(a9,x,y,color);
    end
end

for i = 1: size(Elematrix,1)
    N1 = NodeMatrix(Elematrix(i,2),2:3);
    N2 = NodeMatrix(Elematrix(i,3),2:3);
    N3 = NodeMatrix(Elematrix(i,4),2:3);


%    text(axes,(N1(1)+N2(1)+N3(1))/3 , (N1(2)+N2(2)+N3(2))/3 , num2str(Elematrix(i,1)),...
%               'Color' , 'r' , 'FontSize',15 )
    
    %Polygon Coordinates
    x=[N1(1) N2(1) N3(1)];
    y=[N1(2) N2(2) N3(2)];
     if ismember(5,idx)
  
    %Arrows
    quiver(a5,mean(x),mean(y),principal(i,1)*cos(principal(i,3)),principal(i,1)*sin(principal(i,3)),...
             3 ,'Color' , 'k');
    end
    
     if ismember(6,idx)   

    %Arrows
    quiver(a6,mean(x),mean(y),principal(i,2)*cos(principal(i,3)-pi/2),principal(i,2)*sin(principal(i,3)-pi/2),...;
             3 ,'Color' , 'k');
     end
end     
     




% Xvec = NodeMatrix(1:400,2);
% Yvec = NodeMatrix(1:400,3);
% Xmat = zeros(20);
% c=1;
% for i =1:19
%     Xmat(end-i,:) = Xvec(c:c+19);
%     c=c+20;
% end
% Ymat = zeros(20);
% c =1;
% for i =1:19
%     Ymat(end-i,:) = Yvec(c:c+19);
%     c=c+20;
% end
% 
% 
% 
% cvec = sqrt(displacementMat(1:400,1).^2+displacementMat(1:400,2).^2);
% cmat =zeros(20);
% c = 1;
% for i =1:19
%     cmat(end-i,:) = cvec(c:c+19);
%     c=c+20;
% end
% z = zeros(size(Xmat));
% 
% colormap jet
% surf(Xmat,Ymat,z,cmat, 'FaceColor', 'interp')

function [xmat, ymat] = difs (N1,N2,N3)
%This function calculates the differences xij = xi - xj
%and yij = yi - yj and organises them in a matrix
mat = [N1' N2' N3'];
xmat = repelem(mat(1,:),3,1)' - repelem(mat(1,:),3,1);
ymat = repelem(mat(2,:),3,1)' - repelem(mat(2,:),3,1);
   
end
    



function [A] = Tarea(N1,N2,N3)
    %This function calculates triangle aerea based on coordinates
    
    A = 0.5*abs(N1(1)*N2(2) + N2(1)*N3(2) + N3(1)*N1(2)...
           -N1(2)*N2(1) - N2(2)*N3(1) - N3(2)*N1(1));
end


