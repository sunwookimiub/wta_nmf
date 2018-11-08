function [maxIY,maxIX]=WTA(X,P)
[~,T]=size(X);
[L,M]=size(P);
iP=sub2ind(size(X), repmat(P(:), [T,1]), reshape(repmat((1:T), [numel(P),1]), [T*numel(P),1]));
Y=X(iP);
Y=reshape(Y,[size(P), T]);
[maxY, maxIY]=max(Y,[],2);
maxIY=reshape(maxIY,[L,T]); 
for ii=1:L
    maxIX(ii,:)=P(ii,maxIY(ii,:));
end
end