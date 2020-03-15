

i_set=3; 

Net=Networks(i_set);
pcm=     Net.pcm;  
[N, ~, T]=       size(pcm);
st= Partitions(i_set).contigStates;

alleigs=zeros(10,T); 
Us1=zeros(N,T); Us2=zeros(N,T); Us3=zeros(N,T);
dt=1; 

for t=1:T

A= pcm(:,:,t);
A= A./(1+svds(A,1)); 
A= A-eye(size(A,1));        % Stabilize continuous case
[U, tt] = schur(A,'real');   
eigvals = exp(diag(tt)*dt); 
alleigsB(:,t)=1-eigvals(end-10+1:end).^2;
alleigsA(:,t)=1-eigvals(1:10).^2;
U1=U(:,1:10); U2=U(:,11:end-10); U3=U(:,end-10+1:end); 
Us1(:,t)= mean(U1.^2,2); 
Us2(:,t)= mean(U2.^2,2); 
Us3(:,t)= mean(U3.^2,2); 

end

%%
arrayfun(@(x) mean(mean(Us1(:,st==x))), [1:3])
arrayfun(@(x) mean(mean(Us2(:,st==x))), [1:3])
arrayfun(@(x) mean(mean(Us3(:,st==x))), [1:3])

arrayfun(@(x) mean(alleigsA(:,st==x), 'all'), [1:3]) %transient
arrayfun(@(x) mean(alleigsB(:,st==x), 'all'), [1:3]) %persistent


%%

figure(8); imagesc(Metric_matrices(3).modalCtrl)
figure(2); imagesc(Metric_matrices(3).aveCtrl)
figure(3); imagesc(Metric_matrices(3).tModalCtrl)

figure(4); imagesc(alleigs)
figure(5); 
subplot(1,3,1); imagesc(Us1(:,st==1))
subplot(1,3,2); imagesc(Us1(:,st==2))
subplot(1,3,3); imagesc(Us1(:,st==3))



