function [control_B, sigma]=getSpreadControl(allCoords, i_soz, sigma)
% allCoords- An Nx2 vector of the x,y electrode coords
% i_soz- either a logical vector or array of SOZ indices
% sigma- controls the gaussian dropoff
    
    % Get distance of nodes to closest soz node. 
    dists=min(pdist2(allCoords, allCoords(i_soz,:)), [], 2);
    
    if nargin<3
        % Create gaussian dropoff, choose sigma s.t. full-width of 5mm at half max
        % Since spacing b/w electrodes is 1.3 mm, FWHM= 5mm/1.3= 2.355*sigma
        sigma=(5/1.3)/2.355; 
    end
    
    x = [0:.1:max(dists)+.5];
    y = normpdf(x,0,sigma); y=y./max(y);
    
    xind=arrayfun(@(z)find(z<=x',1, 'first'), dists);
    control_B=diag(y(xind)); 
    
%     figure(3)
%     clf; hold on;
%     scatter(allCoords(:,1), allCoords(:,2), 100, y(xind), 'filled')
%     %scatter(allCoords(i_soz,1), allCoords(i_soz,2), 150*y(xind), 'r', 'filled')
    
    

