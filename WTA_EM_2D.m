function [W,H,errt]=WTA_EM_2D(X,K,P)
% K: number of topics
% P: L X M hash table, where L is the num of permuations and M is the random samples per permutation 

[F,T]=size(X);
[L,M]=size(P);

W=rand(F,K);
W=W./sum(W,1);
H=rand(K,T);
H=H./sum(H,1);
PP=rand([size(X),K]);
PP=PP./sum(PP,3);

maxIter=200;
[~,iX]=CWTA(X(:),P);
% [~,iX2]=WTA(X',P2);
 R=X.*PP;
for ii=1:maxIter
%     R=X.*PP;
    R2=reshape(R,[F*T,K]);
    [~,iR]=CWTA(R2,P);
    match=iR==iX;
%     iR3=reshape(iR, [L,1,K]);

    cnt=zeros(F,T,K);
    for ll=1:L
        [ri, ci]=ind2sub([F,T],iX(ll));
        for kk=1:K          
            if iX(ll)<0
                continue;
            else
                cnt(ri, ci, kk)=cnt(ri, ci, kk)+match(ll,kk);            
            end
        end
    end
        
    idx=find(sum(cnt,3)==0);
    PP2=zeros(F,T);
    PP2(idx)=1/K;
    PP=cnt./(sum(cnt,3)+eps);    
    for kk=1:K        
        PP(:,:,kk)=PP(:,:,kk)+PP2;
    end
    R=X.*PP;
    
    W=reshape(sum(R,2),[F,K]);
    W=W./sum(W,1);    
    H=reshape(sum(R,1),[T,K]);
    H=H';
    H=H./sum(H,1);
    
    
    
    
    
    
%     W=W.*((X./(W*H+eps))*H');
%     W=W./sum(W,1);
%     H=H.*(W'*(X./(W*H+eps)));
%     H=H./sum(H,1);
    errt(ii)=sum(sum(X.*log(X./(W*H+eps)+eps)));
end

end

