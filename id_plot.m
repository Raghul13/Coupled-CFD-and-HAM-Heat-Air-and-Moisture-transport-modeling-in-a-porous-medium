function id_plot(x,y,id_type)

dx = x(2,2)-x(1,1);
dy = y(2,2)-y(1,1);

imagesc(x(1,:),y(:,1)',id_type);
hold on;
for i = 1:length(id_type(:))
    text(x(i),y(i),num2str(id_type(i)),'VerticalAlignment','middle','HorizontalAlignment','center')
end
for i = 1:size(x,1)
    if i == 1
        plot([x(1,1)-dx/2 x(1,end)+dx/2],[y(i,1)-dy/2 y(i,1)-dy/2],'k')
    end
    plot([x(1,1)-dx/2 x(1,end)+dx/2],[y(i,1)+dy/2 y(i,1)+dy/2],'k')
end
for i = 1:size(x,2)
    if i == 1
        plot([x(1,i)-dx/2 x(1,i)-dx/2],[y(1,1)-dy/2 y(end,1)+dy/2],'k')
    end
    plot([x(1,i)+dx/2 x(1,i)+dx/2],[y(1,1)-dy/2 y(end,1)+dy/2],'k')
end
set(gca,'YDir','normal')
axis([x(1)-dx/2 x(end)+dx/2 y(1)-dy/2 y(end)+dy/2])
colormap jet
pause;
close all

end