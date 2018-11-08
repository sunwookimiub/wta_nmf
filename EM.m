function [W,H,errt]=EM(X,K)

[F,T]=size(X);

W=rand(F,K);
W=W./sum(W,1);
H=rand(K,T);
H=H./sum(H,1);

maxIter=100;

for ii=1:maxIter    
    W=W.*((X./(W*H+eps))*H');
    W=W./sum(W,1);
    H=H.*(W'*(X./(W*H+eps)));
    H=H./sum(H,1);
    errt(ii)=sum(sum(X.*log(X./(W*H+eps)+eps)));
end

end