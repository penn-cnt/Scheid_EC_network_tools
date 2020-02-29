function [output]= postHocSignedRank(data, type)
% output: [col1, col2, mean diff, U-val, zval, p-val]

n= size(data,2); 
ctr=1; 
output=zeros(n*(n-1)/2, 6); 

for i=1:n
    for j= i+1:n
        mndiff=mean(data(:,i))-mean(data(:,j));
        if strcmp(type, 'signedrank')
            [p, ~, sts] =signrank(data(:,i), data(:,j)); 
            output(ctr,:)=[i, j, mndiff, sts.signedrank, sts.zval, p]; 
        elseif strcmp(type, 'ttest')
            [~,p, ~,sts] =ttest(data(:,i), data(:,j)); 
            output(ctr,:)=[i, j, mndiff, sts.tstat, sts.df, p]; 
        end
        ctr=ctr+1;
    end

end