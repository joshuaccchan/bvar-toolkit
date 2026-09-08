function plotCI_othercolor(x,y1,y2,option) 

f = fill([x',fliplr(x')], [y1',fliplr(y2')], option, 'FaceAlpha', 0.15, 'EdgeColor', 'none');

end