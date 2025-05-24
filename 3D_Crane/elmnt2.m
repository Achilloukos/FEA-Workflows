function EleMatrixstd= elmnt2(A,L)
b=1;
for i=1:size(A,1)
    for j=i+1:size(A,1)
        
        Vector=A(i,2:4)-A(j,2:4);
        
        Length=sqrt(Vector(1)^2 +Vector(2)^2 + Vector(3)^2);
          if (Length==L && i~=28 && j~=28) || Length==sqrt(2)*L || Length==sqrt(1.5)*L || Length==sqrt(2.25)*L || Length==sqrt(1.25)*L %|| (i==14 && j==26) || (j==27 && i==21) || (j==27 && i==22) || (j==28 && i==13)
            
            EleMatrixstd(b,1:3)=[b i j];
            b=b+1;
        end
    end
end
EleMatrixstd(99:102,1)=[99:102];
EleMatrixstd(99:102,2)=[26 27 27 28];
EleMatrixstd(99:102,3)=[14 21 22 13];
end