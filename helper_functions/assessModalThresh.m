% Assess best modal threshold
load('Data/Networks')
nSets=size(Networks,2);

thresh=(0:.05:.5);
%targ=.75;
all_mx=zeros(nSets,4);

for i_set=1:nSets
    i_set
    mx_thresh=zeros(115,2);
    pcm=Networks(i_set).pcm;

    for i_t= 1:size(pcm,3)
        A=pcm(:,:,i_t);

        [mn, stds]=deal(zeros(length(thresh),2));

        for t=1:length(thresh)
            [mp, mt]=modal_persist_trans(A, thresh(t), 1);
            mn(t,:)= [mean(mp), mean(mt)];
            stds(t,:)= [std(mp), std(mt)];
        end

        [~,is]=max(stds);
        mx_thresh(i_t,:)=thresh(is);

%         figure(1)
%         clf;
%         subplot(121)
%         hold on
%         fill([thresh',flipud(thresh)'],[mn(:,1)-stds(:,1), ...
%             flipud(mn(:,1)+stds(:,1))],'r', 'linestyle','none');
%         errorbar(thresh, mn(:,1), stds(:,1))
%         
%         subplot(122)
%         errorbar(thresh, mn(:,2), stds(:,2))

    end
    
    % Compute centile
%     f=histogram(mx_thresh(:,1),'BinEdges', thresh);
%     nless = sum(f.Values < f.Values(16));
%     nequal = sum(f.Values == f.Values(16));
%     centile1 = 100 * (nless + 0.5*nequal) / length(f.Values);
%     
%     f=histogram(mx_thresh(:,2),'BinEdges', thresh);
%     nless = sum(f.Values < f.Values(16));
%     nequal = sum(f.Values == f.Values(16));
%     centile2 = 100 * (nless + 0.5*nequal) / length(f.Values);
    
    all_mx(i_set,:)= [mean(mx_thresh),std(mx_thresh)]; %, centile1, centile2];
    
end

disp('done')