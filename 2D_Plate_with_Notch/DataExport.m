% Achillefs Fasoulas 6581

function [] = DataExport(filename ,varargin)
% This function takes as input a filename and any number of variables to be
% written in a text file named 'filename'
%Use DataImport.m to read the data

vec = zeros(2*size(varargin,2),1);
c = 1;
for i = 1:size(varargin,2)
    vec(c:c+1) = size(varargin{i});
    c = c + 2;
end
 rows = vec(1:2:end-1);
maxrows = max(rows);
Matrix = zeros(maxrows, sum(vec(2:2:end))+1);


coli = 1;
colf = size(varargin{1},2);

for i = 1:size(varargin,2)
    Matrix(1:size(varargin{i},1),coli:colf) = varargin{i};
    coli = colf + 1;
    if i ~= size(varargin,2)
     colf = colf +size(varargin{i+1},2);
    end

end
Matrix(1:size(vec),end) = vec;

writetable(table(Matrix),filename)



end
