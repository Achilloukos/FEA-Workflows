function OriMatrix = orient(NodeMatrix,EleMatrix,A0,L)
OriMatrix=zeros(size(EleMatrix,1),6); %Columns: ID|Length|cos1|cos2|cos3|Ac
OriMatrix(:,1)=EleMatrix(:,1);
for i=1:size(EleMatrix,1)
    Vector=NodeMatrix(EleMatrix(i,3),2:4)-NodeMatrix(EleMatrix(i,2),2:4);
    Length=sqrt(sum(Vector.^2));
    OriMatrix(i,2)=Length;
    OriMatrix(i,3)=Vector(1)/Length;
    OriMatrix(i,4)=Vector(2)/Length;
    OriMatrix(i,5)=Vector(3)/Length;
    if Length-L<0.1
        OriMatrix(i,6)=1.5*A0;
    else
        OriMatrix(i,6)=0.5*A0;
    end
        
    
end