function plotCI(x,y1,y2,C)
if nargin == 3
    C = [.85 .85 .85];
end    
f = fill([x',fliplr(x')], [y1',fliplr(y2')],C);
set(f,'EdgeColor','none'),
end