% Achillefs Fasoulas 6581

function [NodeMatrix , Elematrix] = MeshgridGenerator_3(xspace, yspace)
%This function creates a xspace by yspace grid of points on the quarter
%geometry and connect the in triangular elements

%Initial Data
[~, H, W , r, t] = PreprocessorData();
y1 = round(7/10*yspace);
y2 = yspace  - y1;
%Coordinates of points on Perimeter

lin1=linspace(0,r,round(xspace*0.4));
lin2=linspace(r,W/2,xspace - length(lin1)+1);

down_side_x = [lin1 lin2(2:end)];

down_side_y=linspace(0,0,xspace);

up_side_x = down_side_x ;

up_side_y=zeros(length(up_side_x),1);

for i = 1:length(up_side_y)
    if up_side_x(i) <=r
        up_side_y(i) = -sqrt(r^2-(up_side_x(i))^2)+H/2-t+r;
    elseif up_side_x(i) > r
        up_side_y(i) = H/2-t+r;
        
    end   
end


%Creating NodeMatrix
xspace2 = length(down_side_x(down_side_x>=r));
NodeMatrix = zeros(xspace*y1+(y2-1)*xspace2,3);
NodeMatrix(:,1) = 1:size(NodeMatrix,1);
c = 1;
for i = 1: xspace %:xspace^2
    NodeMatrix(c:xspace:xspace*y1-xspace+c,3) = linspace(down_side_y(c),up_side_y(c),y1);
    NodeMatrix(c:xspace:xspace*y1-xspace+c,2) = up_side_x(c);
    c = c + 1;
end
c = 1; 
for i = 1:xspace2
    y = linspace(H/2-t+r,H/2,y2);
    NodeMatrix((xspace)*(y1)+c:xspace2:xspace*y1+(xspace2)*y2-xspace2+c-1,3) = y(2:end);
    NodeMatrix((xspace)*(y1)+c:xspace2:xspace*y1+(xspace2)*y2-xspace2+c-1,2) = down_side_x(xspace-xspace2+c);  
    c = c+1;
end


%Creating Element Matrix
Elematrix = zeros(2*(xspace-1)*(y1-1)+2*(xspace2-1)*(y2-1) , 4);
Elematrix(:,1) = 1:size(Elematrix,1);


c = 1;
%Drawing Nodes
plot(NodeMatrix(:,2),NodeMatrix(:,3),'.')
text(NodeMatrix(:,2),NodeMatrix(:,3),num2str(NodeMatrix(:,1)))


for i = 1:2:2*(xspace-1)*(y1-1)
    if mod(c,xspace) == 0
        c=c+1;
    end
    Elematrix(i,2:4) = [c ,c+1 ,c+xspace];
    Elematrix(i+1,2:4) = [c+1 , c+xspace+1, c+xspace];
    c = c+1;
end

c = c+xspace-xspace2+1;
for i = 2*(xspace-1)*(y1-1)+1:2:2*(xspace-1)*(y1-1)+2*(xspace2-1)*(y2-1)
    if mod(c-xspace*y1,xspace2) == 0 
        c=c+1;
    end
    Elematrix(i,2:4) = [c ,c+1 ,c+xspace2];
    Elematrix(i+1,2:4) = [c+1 , c+xspace2+1, c+xspace2];
    c = c+1;
end


%Drawing Elements

for i = 1: size(Elematrix,1)
        N1 = NodeMatrix(Elematrix(i,2),2:3);
        N2 = NodeMatrix(Elematrix(i,3),2:3);
        N3 = NodeMatrix(Elematrix(i,4),2:3);
        
       line([N1(1) N2(1);N2(1) N3(1);N3(1) N1(1)],...
                [N1(2) N2(2);N2(2) N3(2);N3(2) N1(2)],...
                'Color','b')
             
%       text((N1(1)+N2(1)+N3(1))/3 , (N1(2)+N2(2)+N3(2))/3 , num2str(Elematrix(i,1)),...
%           'Color' , 'r' , 'FontSize',15 )
        
end
clear i N1 N2 N3

end