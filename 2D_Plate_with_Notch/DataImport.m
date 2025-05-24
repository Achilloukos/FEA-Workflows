% Achillefs Fasoulas 6581

function [S] = DataImport(filename)
Data = readtable(filename);
Dim = Data{:,end};
Dim = nonzeros(Dim);
L = length(Dim);
for i = 1:L/2
    j = 2*i;
    S{i} = Data{1:Dim(j-1),1:Dim(j)};
    Data(:,1:Dim(j)) = [];
end
   
end