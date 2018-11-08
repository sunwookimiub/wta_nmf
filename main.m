%% create GT basis
N=4;
W3=zeros(N,N,N*2);
% vertical bars
for ii=1:N
    W3(:,ii,ii)=1;
end
% horizontal bars
for ii=1:N
    W3(ii,:,ii+N)=1;
end
figure;
title('Basis in 2D arrays')
for ii=1:2*N
    subplot(2,4,ii);
    colormap('gray');
    imagesc(W3(:,:,ii));    
    axis equal;
    axis off;
end
% vectorize basis images
W=reshape(W3, [N*N, N*2]);
W=W+(normrnd(0,0.01, N*N,N*2)).^2; 
figure;
colormap('gray');
imagesc(W);
axis equal;
axis off;
title('Basis in 1D arrays')

%% Create mixture
T=200;
nActTop=2; % the maximum active topics per documents
H=rand(N*2, T);
for t=1:T
    idx=randperm(N*2,N*2-nActTop);
    H(idx,t)=0;
end
H=(H+normrnd(0,0.01, N*2,T).^2)./sum(H,1);
X=W*H+eps;
X3=reshape(X, [N, N, T]);
figure;
for ii=1:T
    subplot(2,4,ii);
    colormap('gray');
    imagesc(X3(:,:,ii));    
    axis equal;
    axis off;
end
%% EM
[Wem, Hem, errt]=EM(X,N*2);
Wem3=reshape(Wem, [N, N, N*2]);
figure; plot(errt);
figure;
for ii=1:2*N
    subplot(2,4,ii);
    colormap('gray');
    imagesc(Wem3(:,:,ii));    
    axis equal;
    axis off;
end
%% 
% create hash table
L1=500; %2*N*100;
M1=6;
P1=zeros(L1,M1);
for ii=1:L1
    P1(ii, :)=randperm(N*N,M1);
end
L2=500; %T*100;
M2=10;
P2=zeros(L2,M2);
for ii=1:L2
    P2(ii, :)=randperm(T,M2);
end
%%
% [maxIY, maxIX]=CWTA(X,P1);
Wi=rand(size(X,1),N*2);
% Wi(randperm(N*N*N*2,N*2))=1;
Wi=Wi./sum(Wi,1);
% Wi=W;
Hi=rand(N*2,size(X,2));
% Hi(randperm(T*N*2,N*2))=1;
Hi=Hi./sum(Hi,1);
[Wwta, Hwta, errWta]=WTA_EM(X,2*N,P1,P2,Wi,Hi);
Wwta3=reshape(Wwta, [N, N, N*2]);
figure; plot(errWta);
figure;
for ii=1:2*N
    subplot(2,4,ii);
    colormap('gray');
    imagesc(Wwta3(:,:,ii));    
    axis equal;
    axis off;
end

% Xwta=Wwta*Hwta;
% Xwta3=reshape(Xwta, [N, N, T]);
% figure;
% for ii=1:T
%     subplot(2,4,ii);
%     colormap('gray');
%     imagesc(Xwta3(:,:,ii));    
%     axis equal;
%     axis off;
% end

%% WTA EM 2D
% create 2D hash table
L=5000; M=4;
P=zeros(L,M);
for ii=1:L
    P(ii, :)=randperm(N*N*T,M);
end
% [maxIY, maxIX]=WTA(X,P);
[Wwta, Hwta, errWta]=WTA_EM_2D(X,2*N,P);
Wwta3=reshape(Wwta, [N, N, N*2]);
figure; plot(errWta);
figure;
for ii=1:2*N
    subplot(2,4,ii);
    colormap('gray');
    imagesc(Wwta3(:,:,ii));    
    axis equal;
    axis off;
end