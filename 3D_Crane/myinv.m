function [invmat] = myinv(inmat)
    tic
    if mydet(inmat)~= 0
        n = size(inmat,1);
        ADJ = zeros(n);
        for i =1:n
            for j = 1:n
                temp = inmat;
                temp(i,:) = [];
                temp(:,j)=[];
                ADJ(j,i) = (-1)^(i+j)*mydet(temp);
            end
        end
        invmat = (1/mydet(inmat))*ADJ;
        t2 = toc;
        disp(toc)
    
    
    else
        disp("ERROR: Matrix can't be inverted")
        t2 = toc;
        disp(toc)
    end
end

