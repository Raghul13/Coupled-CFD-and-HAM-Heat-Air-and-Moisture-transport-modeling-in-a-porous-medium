function fv_plotting(x,y,phi)

dx = x(2,2)-x(1,1);
dy = y(2,2)-y(1,1);

x_mod = x-dx/2;
y_mod = y-dy/2;
x_mod(:,end+1) = x_mod(:,end)+dx;
y_mod(:,end+1) = y_mod(:,end);
x_mod(end+1,:) = x_mod(end,:);
y_mod(end+1,:) = y_mod(end,:)+dy;

phi(:,end+1) = nan(size(phi,1),1);
phi(end+1,:) = nan(1,size(phi,2));

pc = pcolor(x_mod,y_mod,phi); set(pc,'LineStyle','none');