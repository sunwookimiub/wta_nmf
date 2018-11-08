function [maxIY,maxIX]=CWTA(X,P)
th=0.1;
[~,T]=size(X);
[L,M]=size(P);
iP=sub2ind(size(X), repmat(P(:), [T,1]), reshape(repmat((1:T), [numel(P),1]), [T*numel(P),1]));
Y=X(iP);
Y=reshape(Y,[size(P), T]);
[maxY, maxIY]=max(Y,[],2);
[sortedY, sortedIY]=sort(Y,2,'descend');
nClrWinnerIdx=sortedY(:,1,:)-sortedY(:,2,:)<th;
% maxIY(nClrWinnerIdx)=-rand(numel(find(nClrWinnerIdx)),1);
maxIY=reshape(maxIY,[L,T]); 
for ii=1:L
    maxIX(ii,:)=P(ii,maxIY(ii,:));
end
maxIX(nClrWinnerIdx)=-rand(numel(find(nClrWinnerIdx)),1);
end