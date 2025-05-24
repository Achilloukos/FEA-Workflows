% Achillefs Fasoulas 6581

function [NodeMatrixFinal,ElematrixFinal] = Mirror_X(NodeMatrix,Elematrix)
% This fumction mirrors Nodes & elements around the x axis


NodeMatrixMirrorX=NodeMatrix;
NodeMatrixMirrorX(:,3)=NodeMatrixMirrorX(:,3).*(-1);
NodeMatrixMirrorX(NodeMatrixMirrorX(:,3) == 0 , :) = [];
NodeMatrixMirrorX(:,1)=size(NodeMatrix,1)+1:size(NodeMatrix,1)+size(NodeMatrixMirrorX,1);

NodeMatrixFinal=zeros(size(NodeMatrix,1)+size(NodeMatrixMirrorX,1),3);
NodeMatrixFinal(1:size(NodeMatrix,1),:)=NodeMatrix;
NodeMatrixFinal(size(NodeMatrix,1)+1:end,:)=NodeMatrixMirrorX;
NodeMatrixFinal(:,1) = 1:size(NodeMatrixFinal,1);

% Transforming Elematrix to coordinate form
ElematrixC = zeros(size(Elematrix,1),7);
ElematrixC(:,2:7) = [NodeMatrix(Elematrix(:,2),2:3), NodeMatrix(Elematrix(:,3),2:3) NodeMatrix(Elematrix(:,4),2:3)];

% Mirroring Elematrix
ElematrixCMirror = ElematrixC;
ElematrixCMirror(:,[3 5 7]) = -(1)*ElematrixCMirror(:,[3 5 7]);


%Returning to node form
ElematrixMirror = zeros(size(ElematrixCMirror,1),4);
for i = 1: size(ElematrixCMirror)
    for j = 1:3
        Ni = ElematrixCMirror(i,2*j:2*j+1);
        P1x = find(NodeMatrixFinal(:,2) == Ni(1));
        P1y = find(NodeMatrixFinal(:,3) == Ni(2));
        NIDj = intersect(P1x,P1y);
        ElematrixMirror(i,j+1) = NIDj;
    end
end

ElematrixFinal = [Elematrix ; ElematrixMirror(:,[1 4 2 3])];
ElematrixFinal(:,1) = 1:size(ElematrixFinal,1);

% % % EleMatrixMirrorX=zeros(size(EleMatrix));
% % % EleMatrixMirrorX(:,1)=EleMatrix(:,1)+size(EleMatrix,1);
% % % 
% % % c=NodeMatrixMirrorX(1,1);
% % % 
% % % for i=1:2:size(EleMatrix,1)
% % %      if mod(c,xspace)==0
% % %         c=c+1;
% % %      end
% % %     EleMatrixMirrorX(i,2:4)=[c c+xspace c+xspace+1];
% % %     EleMatrixMirrorX(i+1,2:4)=[c+1 c c+xspace+1];
% % %     c=c+1;
% % % end
% % % Ele=zeros(2*size(EleMatrix,1),4);
% % % Ele(1:size(EleMatrix,1),:)=EleMatrix;
% % % Ele(size(EleMatrix,1)+1:2*size(EleMatrix,1),:)=EleMatrixMirrorX;
% % % 
% % % EleMatrixFinal=Ele;
% % % for i=1:xspace
% % %     EleMatrixFinal(EleMatrixFinal==NodeMatrix(end,1)+i)=i;
% % % end
% % %  EleMatrixFinal(:,1)=1:2*size(EleMatrix,1);
% % %  
% % %  for i=size(NodeMatrix,1)+xspace+1:size(NodeMatrixFinal,1)+xspace
% % %     EleMatrixFinal(EleMatrixFinal==i)=EleMatrixFinal(EleMatrixFinal==i)-xspace;
% % %  end
% % %   EleMatrixFinal(:,1)=1:2*size(EleMatrix,1);
% % % end
end
