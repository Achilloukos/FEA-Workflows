function [NodeMatrixstd] = standartfloornodes(floorN,initfloor,L)
%%This function creates the node matrix for standart cubical elements.
%Inputs: Number of floors, Initial floor , characteristic cube dimention
%Outputs: Node matrix   ID|x|y|z|
NodeMatrixstd=zeros(4*floorN+5,4);
initfloor=1;
for i=3:4:4*floorN+2
 
   %Insert x values
   NodeMatrixstd(i,2)=L/2;
   NodeMatrixstd(i+1,2)=L/2;
   NodeMatrixstd(i+2,2)=-L/2;
   NodeMatrixstd(i+3,2)=-L/2;
   %inssert y valyues
   NodeMatrixstd(i,3)=-L/2;
   NodeMatrixstd(i+1,3)=L/2;
   NodeMatrixstd(i+2,3)=L/2;
   NodeMatrixstd(i+3,3)=-L/2;
   %insert z values
   NodeMatrixstd([i:i+3],4)=initfloor*L;
   initfloor=initfloor+1;
end
NodeMatrixstd(:,1)=[1:size(NodeMatrixstd,1)];
end

