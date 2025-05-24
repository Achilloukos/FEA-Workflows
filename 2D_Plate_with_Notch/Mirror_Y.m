% Achillefs Fasoulas 6581

function[NodeMatrixFinal,ElematrixFinal] = Mirror_Y(NodeMatrix,Elematrix)
% This fumction mirrors Nodes & elements around the y axis


NodeMatrixMirrorY=NodeMatrix;
NodeMatrixMirrorY(:,2)=NodeMatrixMirrorY(:,2).*(-1);
NodeMatrixMirrorY(NodeMatrixMirrorY(:,2) == 0 , :) = [];
NodeMatrixMirrorY(:,1)=size(NodeMatrix,1)+1:size(NodeMatrix,1)+size(NodeMatrixMirrorY,1);


NodeMatrixFinal=zeros(size(NodeMatrix,1)+size(NodeMatrixMirrorY,1),3);
NodeMatrixFinal(1:size(NodeMatrix,1),:)=NodeMatrix;
NodeMatrixFinal(size(NodeMatrix,1)+1:end,:)=NodeMatrixMirrorY;
NodeMatrixFinal(:,1) = 1:size(NodeMatrixFinal,1);

% Transforming Elematrix to coordinate form
ElematrixC = zeros(size(Elematrix,1),7);
ElematrixC(:,2:7) = [NodeMatrix(Elematrix(:,2),2:3), NodeMatrix(Elematrix(:,3),2:3) NodeMatrix(Elematrix(:,4),2:3)];

% Mirroring Elematrix
ElematrixCMirror = ElematrixC;
ElematrixCMirror(:,[2 4 6]) = -(1)*ElematrixCMirror(:,[2 4 6]);


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

% % % EleMatrixMirrorY=zeros(size(EleMatrix));
% % % EleMatrixMirrorY(:,1)=EleMatrix(:,1)+size(EleMatrix,1);
% % % 
% % % c=NodeMatrixMirrorY(1,1);
% % % 
% % % for i=1:2:size(EleMatrix,1)
% % %      if mod(c,xspace)==0
% % %         c=c+1;
% % %      end
% % %     EleMatrixMirrorY(i,2:4)=[c c+xspace c+xspace+1];
% % %     EleMatrixMirrorY(i+1,2:4)=[c+1 c c+xspace+1];
% % %     c=c+1;
% % % end
% % % 
% % % 
% % % Ele = zeros(2*size(EleMatrix,1),4);
% % % Ele(1:size(EleMatrix,1),:) = EleMatrix;
% % % Ele(size(EleMatrix,1)+1:2*size(EleMatrix,1),:) = EleMatrixMirrorY;
% % % 
% % % 
% % % ElematrixFinal = Ele;
% % % c = 1;
% % % for i = 1:xspace:xspace*yspace+1
% % %     for j = i:1:xspace+i-1
% % %         if i == j
% % %             ElematrixFinal(ElematrixFinal == size(NodeMatrix,1)+i) = i;
% % %         else
% % %             ElematrixFinal(ElematrixFinal == size(NodeMatrix,1)+j) = j +size(NodeMatrix,1)-c;
% % %         end
% % %     end
% % %     c = c+1;
% % % end
% % % 
% % % ElematrixFinal(ElematrixFinal > size(NodeMatrix,1)+xspace*yspace - yspace) = ...
% % %     ElematrixFinal(ElematrixFinal > size(NodeMatrix,1)+xspace*yspace - yspace) - yspace;
% % % 
% % % ElematrixFinal(:,1) = 1:size(ElematrixFinal,1);





end