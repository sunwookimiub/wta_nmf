function [W,H,errt]=WTA_EM(X,K,P1,P2, W, H)
% K: number of topics
% P: L X M hash table, where L is the num of permuations and M is the random samples per permutation 

[F,T]=size(X);
[L1,M1]=size(P1);
[L2,M2]=size(P2);

% W=rand(F,K);
% W=W./sum(W,1);
% H=rand(K,T);
% H=H./sum(H,1);
% PP=zeros([size(X),K]);

maxIter=200;
[~,iX1]=WTA(X,P1);
[~,iX2]=WTA(X',P2);
alpha=0.5;
for ii=1:maxIter
    
    [~,iW]=WTA(W,P1);
    iW3=reshape(iW, [L1,1,K]);
    match3=iW3==iX1;
    
    cnt13=zeros(F,T,K);
    for ll=1:L1
        for tt=1:T
            if iX1(ll,tt)<0
                continue;
            else
                for kk=1:K
                    if iW3(ll,1,kk)<0
                        continue;
                    else
                        cnt13(iX1(ll,tt),tt,kk)=cnt13(iX1(ll,tt),tt,kk)+match3(ll,tt,kk);
                    end
                end
            end
        end
    end
    
    [~,iH]=WTA(H',P2);
    iH3=reshape(iH, [L2,1,K]);
    match3=iH3==iX2;
    
    cnt23=zeros(F,T,K);
    
    for ll=1:L2
        for ff=1:F
            if iX2(ll,ff)<0
                continue;
            else
                for kk=1:K
                    if iH3(ll,1,kk)<0
                        continue;
                    else
                        cnt23(ff,iX2(ll,ff),kk)=cnt23(ff,iX2(ll,ff),kk)+match3(ll,ff,kk);
                    end
                end
            end
        end
    end
%     
%     idx=find(sum(cnt13,3)==0);
%     PP2=zeros(F,T);
%     PP2(idx)=1/K;
%     maxCnt=max(max(sum(cnt13,3)));
%     sumCnt=sum(cnt13,3);
%     cnt13=cnt13+(maxCnt-sumCnt)/K;
    PPW=(cnt13)./(sum(cnt13,3));
    PPW(isnan(PPW))=1/K;
%     for kk=1:K        
%         PPW(:,:,kk)=PPW(:,:,kk)+PP2;
%     end
    
%     idx=find(sum(cnt23,3)==0);
%     PP2=zeros(F,T);
%     PP2(idx)=1/K;        
%     maxCnt=max(max(sum(cnt23,3)));
%     sumCnt=sum(cnt23,3);
%     cnt23=cnt23+(maxCnt-sumCnt)/K;
    PPH=(cnt23)./(sum(cnt23,3));    
    PPH(isnan(PPH))=1/K;
%     for kk=1:K        
%         PPH(:,:,kk)=PPH(:,:,kk)+PP2;
%     end
%     R=X.*(PPW);
%     R=X.*(PPH);
    R=X.*(alpha*PPW+(1-alpha)*PPH);
    
%     figure;
%     for ii=1:K
%         subplot(2,4,ii);
%         colormap('gray');
%         imagesc(reshape(R(:,1,ii),[4,4]));    
%         axis equal;
%         axis off;
%     end
    
    W=reshape(sum(R,2),[F,K]);
    W=W./sum(W,1);    
    H=reshape(sum(R,1),[T,K]);
    H=H';
    H=H./sum(H,1);
    
    
%     idxXR=sub2ind([F,T], iX1(:)', kron(1:T, ones(1,L1)));
%     idxCnt=sub2ind([L1,T], kron(ones(1,T), 1:L1), kron(1:T, ones(1,L1)));
%     
%     cnt3=zeros(F,T,K);
%     for kk=1:K
%         cnt2=zeros(F,T);
%         match2=match3(:,:,kk);
%         cnt2(idxXR)=cnt2(idxXR)+match2(idxCnt);
%         cnt3(:,:,kk)=cnt2;
%     end
%     
%     [~,iH]=WTA(H',P2);
%     iH3=reshape(iH, [L2,1,K]);
%     match3=iH3==iX2;
%     idxXR=sub2ind([T,F], iX2(:)', kron(1:F, ones(1,L2)));
%     idxCnt=sub2ind([L2,F], kron(ones(1,F), 1:L2), kron(1:F, ones(1,L2)));
%     
%     for kk=1:K
%         cnt2=zeros(T,F);
%         match2=match3(:,:,kk);
%         cnt2(idxXR)=cnt2(idxXR)+match2(idxCnt);
%         cnt3(:,:,kk)=cnt3(:,:,kk)+cnt2';
%     end
%         
%     idx=find(sum(cnt3,3)==0);
%     PP2=zeros(F,T);
%     PP2(idx)=1/K;
%     PP=cnt3./(sum(cnt3,3)+eps);    
%     for kk=1:K        
%         PP(:,:,kk)=PP(:,:,kk)+PP2;
%     end
%     R=X.*PP;
%     
%     W=reshape(sum(R,2),[F,K]);
%     W=W./sum(W,1);    
%     H=reshape(sum(R,1),[T,K]);
%     H=H';
%     H=H./sum(H,1);
%     
    
    
    
    
    
%     W=W.*((X./(W*H+eps))*H');
%     W=W./sum(W,1);
%     H=H.*(W'*(X./(W*H+eps)));
%     H=H./sum(H,1);
    errt(ii)=sum(sum(X.*log(X./(W*H+eps)+eps)));
end

end

