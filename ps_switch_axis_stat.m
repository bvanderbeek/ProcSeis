function ps_switch_axis_stat(Ax,stat)

cla(Ax);
if stat == 0
    fill(Ax,[0,1,1,0,0],[0,0,1,1,0],'g');
elseif stat == 1
    fill(Ax,[0,1,1,0,0],[0,0,1,1,0],'r');
elseif stat == 2
    fill(Ax,[0,1,1,0,0],[0,0,1,1,0],'y');
end
Ax.XLim    = [0,1];
Ax.YLim    = [0,1];
Ax.Visible = 'off';
drawnow;