% Display single trial

i_preict=38;
i_ict=26;

mn_amount= zeros(12, 2); 

for ipt=(1:12)
    
i_preict=ipt+12+24;
i_ict=ipt+24;

ps=zeros(1000,2);
b=zeros(1,2);
percentRank = @(YourArray, TheProbes) reshape( mean( bsxfun(@le,...
    YourArray(:), TheProbes(:).') ) * 100, size(TheProbes));
l_preict=length(Partitions(i_preict).states);
l_ict=length(Partitions(i_ict).states);


for i=1:10000
    [vict, nict]=RunLength(Partitions(i_ict).states(randperm(l_ict)));
    [vpre, npreict]=RunLength(Partitions(i_preict).states(randperm(l_preict)));
    
    % maximum measures
%      mx_ict= arrayfun(@(x)max(nict(vict==x)), unique(vict));
%      mx_preict= arrayfun(@(x)max(npreict(vpre==x)), unique(vpre));
%      ps(i,:)=[mean(mx_ict),mean(mx_preict)];
    
    % max overall
     ps(i,:)=[max(nict),max(npreict)];

end
disp('done!')

[vbict,nbict]=RunLength(Partitions(i_ict).states);
[vbpre,nbpreict]=RunLength(Partitions(i_preict).states);

% bmx_ict= arrayfun(@(x)max(nbict(vbict==x))/sum(nbict(vbict==x)), unique(vbict));
% bmx_preict= arrayfun(@(x)max(nbpreict(vbpre==x))/sum(nbpreict(vbpre==x)), unique(vbpre));
% b=[mean(bmx_ict),mean(bmx_preict)];

% maximum measures
b=[max(nbict), max(nbpreict)];

rnk_ict= percentRank(ps(:,1),b(1,1))
rnk_preict= percentRank(ps(:,2),b(1,2))

%
figure(5)
clf; hold on

for typ=[1,2]
    subplot(1,2,typ)
    hold on;
    bns=histogram(ps(:,typ)', 'Normalization','probability');
    stem(b(1,typ), .1); 
    title(sprintf('%s, percentile: %0.2f', Partitions(i_ict).ID, percentRank(ps(:,typ),max(b(1,typ)))))

    yyaxis right
    pd=fitdist(ps(:,typ), 'poisson');
    y=pdf(pd, ps(:,typ));
    plot(ps(:,typ), y, 'o')
    
    mn_amount(ipt, typ)=b (1,typ)/pd.lambda;
    
end

end
%%

figure(6); clf;
imagesc(Partitions(i_ict).states)
title(sprintf('max run length: %d', sum(b1(1:3))))

