% Achillefs Fasoulas 6581

function [s, H, W , r, t] = PreprocessorData()
% q=input('Insert AEM: ');
q=6581;
q=num2str(q);
A=str2num(q(1));
B=str2num(q(2));
C=str2num(q(3));
D=str2num(q(4));

s = 4;
H = 100;
W = 170;
r = 5+ 10*(10*D + C)/100;
t = 15 +15*(10*B+A)/100;

end