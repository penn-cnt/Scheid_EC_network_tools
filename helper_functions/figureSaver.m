% Save Figures how you want them
figs = findobj('Type', 'figure');

singCol= 3.125; %PNAS single column width
doubCol= 7.5;   %PNAS double column width

for ff= figs'
           
    if strcmp(ff.Name, 'methodsMetric')
        formatFigure(ff,singCol,3.63)
    
    % Odd figures are column figs
    elseif mod(ff.Number,2)== 1
        formatFigure(ff,singCol,3)
    else
        % Even figures are row figs
        formatFigure(ff,doubCol,3)
    end
    
    %saveas(ff, sprintf('Figs/%s.png', ff.Name))
end


%%
function formatFigure(f, w, h, units)
% input f: figure handle\
% w, h: width and heigh of figure in inches

fsize = 12;
%get axes of children
if nargin < 4
    % Scales input
    w=w/14;
    h=h/8.8;
    units='normalized';
end

a = get(f,'Children');
for i = 1:numel(a)
    c=class(a(i));
    
    % Mess with Axes
    if contains(c, 'Axes')
        set(a(i), ...
          'Box'         , 'off'     , ...
          'XMinorTick'  , 'off'      , ...
          'YMinorTick'  , 'off'      , ...
          'YGrid'       , 'off'      , ...
          'FontSize'    , fsize);

        set( a(i)                      , ...
            'FontName'   , 'Helvetica' );
        hXLabel = get(a(i),'XLabel');
        hYLabel = get(a(i),'YLabel');
        set([hXLabel, hYLabel]  , ...
            'FontSize'   , fsize         );
       
    end

    % Mess with legend
    if contains(c, 'Legend')
        set(a(i),'FontSize',fsize, 'Position', 'Best')
        % lt = get(l,'Title');
        % set(lt,'String','Channels')
        % set(lt,'FontSize',fsize-8);
    end
end

set(f,'PaperPositionMode','auto')  
ppos = get(f,'PaperPosition');
su = get(f,'Units');
pu = get(f,'PaperUnits');  
set(f,'Units',pu);
spos = get(f,'Position');
set(f,'Position',[spos(1) spos(2) ppos(3) ppos(4)])
set(gcf,'Units',su)  

set(f, 'Units', units, 'OuterPosition', [0, 0, w, h])

t = get(gca,'Title');
set(t,'FontSize',fsize);

%set edge color to flat
try
    a = get(gca,'Children');
    set(a,'EdgeColor','flat')
catch
end
end


     