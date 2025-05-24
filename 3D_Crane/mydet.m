function [DET] = mydet(inmat)
%mydet calcultes the determinant of a square matrix 
if size(inmat,1) == size(inmat,2)
    if size(inmat,1)==2
        DET = inmat(1)*inmat(4)-inmat(2)*inmat(3);
    else
        n = size(inmat,1);
        DET = 0;
        for i = 1:n
            temp = inmat;
            temp(1,:) = [];
            temp(:,i)=[];
            DET = DET + (-1)^(i+1)*mydet(temp)*inmat(1,i);
            
        end
       
      
    end
else
    disp("ERROR: matrix is not square")
end

end

