function OriMatrix = orient_2(NodeMatrix,EleMatrix,A0,L)
OriMatrix=zeros(size(EleMatrix,1),12); %Columns: ID|Length|Ac|cos(x,X)|cos(x,Y)|cos(x,Z)|
OriMatrix(:,1)=EleMatrix(:,1);
% NodeMatrix=zeros(size(NodeMatrix));
% NodeMatrix(
for i=1:size(EleMatrix,1)
    Vector21=NodeMatrix(EleMatrix(i,3),2:4)-NodeMatrix(EleMatrix(i,2),2:4);
    v0=[1000*rand() 1000*rand() 1000*rand()];
    v1=cross(v0,Vector21);
    v2=cross(v1,Vector21);
    P3=NodeMatrix(EleMatrix(i,2),2:4)+Vector21*0.5+v1;
    Vector31=P3-NodeMatrix(EleMatrix(i,2),2:4);
    Length=sqrt(sum(Vector21.^2));
    OriMatrix(i,2)=Length;
    OriMatrix(i,4)=Vector21(1)/Length;
    OriMatrix(i,5)=Vector21(2)/Length;
    OriMatrix(i,6)=Vector21(3)/Length;
    lx=OriMatrix(i,4)
    mx=OriMatrix(i,5)
    nx=OriMatrix(i,6)
    
    A123=sqrt((Vector21(2)*Vector31(3)-Vector31(2)*Vector21(3))^2+(Vector21(3)*Vector31(1)-Vector31(3)*Vector21(1))^2+(Vector21(1)*Vector31(2)-Vector31(1)*Vector21(2))^2);
    OriMatrix(i,7)=(Vector21(2)*Vector31(3)-Vector31(2)*Vector21(3))/A123;
    OriMatrix(i,8)=(Vector21(3)*Vector31(1)-Vector31(3)*Vector21(1))/A123;
    OriMatrix(i,9)=(Vector21(1)*Vector31(2)-Vector31(1)*Vector21(2))/A123;
    lz=OriMatrix(i,7);
    mz=OriMatrix(i,8);
    nz=OriMatrix(i,9);
    
    OriMatrix(i,10)=mz*nx-nz*mx;
    OriMatrix(i,11)=nz*lx-lz*nx;
    OriMatrix(i,12)=lz*mx-mz*lx;
    ly=OriMatrix(i,10);
    my=OriMatrix(i,11);
    ny=OriMatrix(i,12);
    
    if Length>2*L
    
        OriMatrix(i,3)=0.5*A0;
    end
    
    
end