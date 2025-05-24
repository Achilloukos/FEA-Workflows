% Achillefs Fasoulas 6581

[idx ,tf] = listdlg('ListString', {'Tension' ,'Bending'}, 'SelectionMode' , 'single','ListSize',[150,100],...
               'PromptString', 'Select type of load');

if ~tf
    error('No Load Selected')
end

answer = inputdlg({'Enter number of nodes in x axis:','Enter number of nodes in y axis:'},'Preproccessor',[1 35]  ,{'11','11'});
if isempty(answer)
    error('No mesh size selected')
elseif ~min(cellfun(@str2num,answer)>0)
    error ('Number of Nodes must be Positive')
end 
           
           
[s, H, W , r, t] = PreprocessorData();

fig = figure(1);


xspace = round(str2num(answer{1}));
yspace = round(str2num(answer{2}));
[NL,EL] = MeshgridGenerator_3(xspace,yspace);
[NL,EL] = Mirror_Y(NL,EL);
[NL,EL] = Mirror_X(NL,EL);
plot(NL(:,2),NL(:,3),'.')
%text(NL(:,2),NL(:,3),num2str(NL(:,1)))
NodeMatrix = NL;
EleMatrix = EL;
 


%Drawing Elements

for i = 1: size(EleMatrix,1)
        N1 = NodeMatrix(EleMatrix(i,2),2:3);
        N2 = NodeMatrix(EleMatrix(i,3),2:3);
        N3 = NodeMatrix(EleMatrix(i,4),2:3);
        
       line([N1(1) N2(1);N2(1) N3(1);N3(1) N1(1)],...
                 [N1(2) N2(2);N2(2) N3(2);N3(2) N1(2)],...
                 'Color','b')
%              
%       text(axes,(N1(1)+N2(1)+N3(1))/3 , (N1(2)+N2(2)+N3(2))/3 , num2str(EleMatrix(i,1)),...
%           'Color' , 'r' , 'FontSize',15 )
        
end

% BCMat UI
userin = questdlg('Do you want to manually select points for boundary conditions? ', 'Boundary Conditions' , 'Yes', 'No ','No ');

BCMat=zeros(size(NodeMatrix));
BCMat(:,1)=1:size(BCMat,1);
if userin == 'Yes' %#ok<BDSCA>
text(NodeMatrix(:,2),NodeMatrix(:,3),num2str(NodeMatrix(:,1)))
uiwait(msgbox('Select Nodes to lock in X direction and press ENTER',''))
[xi, yi] = getpts(fig);
SelectedId = zeros(size(xi));
for i = 1:length(xi)    
   [~ , SelectedId(i)] = min((NodeMatrix(:,2)-xi(i)).^2 +(NodeMatrix(:,3)-yi(i)).^2);

end
BCMat(SelectedId,2) =1;

uiwait(msgbox('Select Nodes to lock in Y direction and press ENTER',''))
[xi, yi] = getpts(fig);
SelectedId = zeros(size(xi));
for i = 1:length(xi)    
   [~ , SelectedId(i)] = min((NodeMatrix(:,2)-xi(i)).^2 +(NodeMatrix(:,3)-yi(i)).^2);

end
BCMat(SelectedId,3) =1;
end

%Automatic Mode 


if userin=='No '
    Px=find(NodeMatrix(:,2)==-W/2);
    BCMat(Px,2)=1;
    BCMat(1,3)=1;
end

%External Load Matrix

LoadMat = zeros(size(NodeMatrix,1),3);
RightID = find(NodeMatrix(:,2) == W/2);
[~ ,RightIDy] = sort(NodeMatrix(RightID,3));
RightID =RightID(RightIDy);
Deltay = abs(NodeMatrix(RightID(2:end),3) - NodeMatrix(RightID(1:end-1),3));



if idx == 1
    F = 1000; %[N] 
    Fdist = F*(Deltay/H);
    for i = 1:size(RightID,1)-1
        LoadMat(RightID(i),2) = LoadMat(RightID(i),2) + Fdist(i)/2;
        LoadMat(RightID(i+1),2) = LoadMat(RightID(i+1),2) + Fdist(i)/2;    
    end
    
end

if idx == 2
   M = 10000; %N*mm
   for i = 1:size(RightID,1)-1
       
       f1 = (12*M)/(H^3)* NodeMatrix(RightID(i),3);
       f2 = (12*M)/(H^3)* NodeMatrix(RightID(i+1),3);
       Fi = (f1+f2)*Deltay(i)*0.5;
       LoadMat(RightID(i),2) = LoadMat(RightID(i),2) + Fi/2; 
       LoadMat(RightID(i+1),2) = LoadMat(RightID(i+1),2) + Fi/2;
        
   end   
    
end
LoadMat(:,1) = 1:size(LoadMat,1);

clear RightID RightIDy Deltay Fdist i
DataExport('PreProcessor_6581.txt',NodeMatrix , EleMatrix , BCMat , LoadMat)

% save('PreProcessor.mat','NodeMatrix','EleMatrix','LoadMat','BCMat')

uiwait(msgbox('Preproccessor ran successfully continuing to Solver','Progress'))

close all

