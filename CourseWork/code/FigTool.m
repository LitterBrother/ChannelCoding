function FigTool(p)                   %w/h
width = 800;            %figure width
d_left_down = 800;      %distance to the left down cornner of the screen
d_bottom = 50;         %distance to the bottom of the screen
height = width/p;       %figure height
set(gcf,'position',[d_left_down,d_bottom,width,height]);
set(gca,'Position',[.12 .12 0.8 0.8]);
set(gca,'FontName','Times New Roman','FontSize',22);
set(gca,'linewidth',1.2);
end