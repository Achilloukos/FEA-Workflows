function [Mat] = RotMat(thetax, thetay,thetaz)
%This function calculates the rotation matrix in general 3d cartesian
%coodrinates 
%INPUTS IN DEGREEEEEEEEES!!!!!!
Rx = [ 1 0 0; 0 cosd(thetax) -sind(thetax) ;0 sind(thetax) cosd(thetax) ];
Ry = [cosd(thetay) 0 sind(thetay) ; 0 1 0 ; -sind(thetay) 0 cosd(thetay)];
Rz = [cosd(thetaz) -sind(thetaz) 0 ; sind(thetaz) cosd(thetaz) 0 ; 0 0 1];

Mat = Rz*Ry*Rx;


end

