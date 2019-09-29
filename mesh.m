function [xMesh,yMesh,Nxt,Nyt,Nxu,Nyu,Nxp,Nyp,Nxv,Nyv,xu_grid,xv_grid,yu_grid,yv_grid,DX,DY,DXv,DYv,u_type,v_type,p_type,xu,xv,xCenter,yCenter,xp,yp,xt,yt,yu,yv,T_type,t_type,dxt,dyt,Tc_type,Tu_type,Tb_type,Tr_type,Tl_type]=mesh

deltaX = 0.1;
Nx =32;

% Defining the domain - x=[0-14] y=[0-9]

x1_start = [4 0];
x1_end = [0 0];
[xVertex1, alph1] = geoDist(x1_start,x1_end,deltaX,Nx);%20

x2_start = [4 0];
x2_end = [4.5 0];
[xVertex2,alpha2] = geoDist(x2_start,x2_end,deltaX,Nx/8);%4

x3_start = [5 0];
x3_end = [4.5 0];
[xVertex3,alpha3] = geoDist(x3_start,x3_end,deltaX,Nx/8);%4

x4_start = [5 0];
x4_end = [6.5 0]; %7
[xVertex4,alpha4] = geoDist(x4_start,x4_end,deltaX,Nx/2);

x5_start = [8 0];
x5_end = [6.5 0];
[xVertex5,alpha6] = geoDist(x5_start,x5_end,deltaX,Nx/2);%4

x6_start = [8 0];
x6_end = [9.5 0];
[xVertex6,alpha7] = geoDist(x6_start,x6_end,deltaX,Nx/2);%4

x7_start = [10 0];
x7_end = [9.5 0];
[xVertex7,alpha8] = geoDist(x7_start,x7_end,deltaX,Nx/6);

x8_start = [10 0];
x8_end = [14 0];
[xVertex8,alpha8] = geoDist(x8_start,x8_end,deltaX,Nx);


          

xMesh = [xVertex1(end:-1:1,1)' xVertex2(2:end,1)' xVertex3(end-1:-1:1,1)' xVertex4(2:end,1)' xVertex5(end-1:-1:1,1)' xVertex6(2:end,1)' xVertex7(end-1:-1:1,1)' xVertex8(2:end,1)'];


% dy=0.15;
% ym1=0:dy:0.5;
% %ym1(end)=0.5;
 deltaY = 0.1;
% %ym=0:deltaY:0.2;

Ny = 32;



y1_start = [0, 1];
y1_end = [0, 0];

[yVertex1, alph2] = geoDist(y1_start, y1_end, deltaY,Ny/4);



y2_start = [0, 1];
y2_end = [0, 2];
yVertex2 = geoDist(y2_start, y2_end, deltaY, Ny/4);


y3_start = [0, 3];
y3_end = [0, 2];
[yVertex3,alph3] = geoDist(y3_start, y3_end, deltaY, Ny/4);

y4_start = [0, 3];
y4_end = [0, 4];
yVertex4 = geoDist(y4_start, y4_end, deltaY, Ny/4);

y5_start = [0, 5];
y5_end = [0, 4];
[yVertex5,alph5] = geoDist(y5_start, y5_end, deltaY, Ny/4);

y6_start = [0, 5];
y6_end = [0, 7];
[yVertex6,alph6] = geoDist(y6_start, y6_end, deltaY, Ny/2);


y7_start = [0, 9];%9-8
 y7_end = [0, 7];
 [yVertex7,alph7] = geoDist(y7_start, y7_end, deltaY, Ny/2);

yMesh =  [yVertex1(end:-1:1,2)' yVertex2(2:end,2)' yVertex3(end-1:-1:1,2)' yVertex4(2:end,2)' yVertex5(end-1:-1:1,2)' yVertex6(2:end,2)' yVertex7(end-1:-1:1,2)'];

% Distance between the center of the cell
xCenter = (xMesh(1:end-1) + xMesh(2:end))/2;
yCenter = (yMesh(1:end-1) + yMesh(2:end))/2;
[xGrid, yGrid] = meshgrid(xMesh, yMesh);
%P_grid
xt=xCenter;
yt=yCenter;
dx_t=xt(end)-xt(end-1);
xt(end+1)=xt(end)+dx_t; %adding ghost cell to the right
xt=[xt(1)-dx_t xt];% adding ghost cell to the left
Xt=[xt(1)-dx_t xt];
Yt=[abs(yt(1)-yt(2)) yt];
Nxt=length(xt);
Nyt=length(yt);
[dxt,dyt] = meshgrid(abs(diff(Xt)),abs(diff(Yt)));
Nxp=length(xCenter);
Nyp=length(yCenter);
[dxp, dyp] = meshgrid(abs(diff(xCenter)), abs(diff(yCenter)));

[xp,yp]=meshgrid(xCenter,yCenter);
%U_grid=xMesh
xu=xMesh;
yu=yCenter;
Nxu=length(xu);
Nyu=length(yu);
[xu_grid,yu_grid]=meshgrid(xu,yu);
[Dx1,Dy1]=meshgrid(abs(diff(xMesh)),abs(diff(yCenter))); %distance between the center points

%V_grid=yMesh
yv=yMesh;
xv=xCenter;

dx_v=xv(end)-xv(end-1);
xv(end+1)=xv(end)+dx_v; %adding ghost cell to the right
xv=[xv(1)-dx_v xv];% adding ghost cell to the left
Nxv=length(xv);
Nyv=length(yv);
[Dx2,Dy2]=meshgrid(abs(diff(xv)),abs(diff(yv)));
[xv_grid,yv_grid]=meshgrid(xv,yv);

%Identifying the distance between the vertices of U cell
 uxvertex=[xCenter(2) xCenter xCenter(end-1)];
 [DX,DY]=meshgrid(abs(diff(uxvertex)),abs(diff(yMesh)));
%Identifying the distance between the vertices of V cell
 vyvertex=[yCenter(2) yCenter yCenter(end-1)];
 vxvertex=[xu(2) xu xu(end-1)];
 [DXv,DYv]=meshgrid(abs(diff(vxvertex)),abs(diff(vyvertex)));
 
%  dx_v=xv(end)-xv(end-1);
% xv(end+1)=xv(end)+dx_v; %adding ghost cell to the right
% xv=[xv(1)-dx_v xv];% adding ghost cell to the left

    figure(1);
  plot(xGrid,yGrid,'.-k','MarkerSize',20,'LineWidth',2); hold on;
  
  plot(xGrid',yGrid','.-k','MarkerSize',20,'LineWidth',2);
    fill([4,10,10,9,9,5,5,4,4],[1,1,3,3,5,5,3,3,1],'white')
%      figure(2);
    % plot(xu_grid,yu_grid,'.-b','MarkerSize',20,'LineWidth',2); 
    %plot(xu_grid',yu_grid','.-b','MarkerSize',20,'LineWidth',2)
%     fill([4,10,10,9,9,5,5,4,4],[1,1,3,3,5,5,3,3,1],'white')
% %     
% % figure(3);
  %plot(xv_grid,yv_grid,'.-g','MarkerSize',20,'LineWidth',2);
%     
   % plot(xv_grid',yv_grid','.-g','MarkerSize',20,'LineWidth',2);
% fill([4,10,10,9,9,5,5,4,4],[1,1,3,3,5,5,3,3,1],'white')   


 u_type = zeros(Nyu,Nxu);
 
%Identifying the obstruction coordinates
 fx1 = find(xu==4);
 fx2= find(xu==10);
 fx3=find(xu==5); %5
 fx4=find(xu==8); %9
 fy1=find(yu>1,1);
 fy2=find(yu>3,1);
 fy3=find(yu>5,1);
 


% u_type(fy1-1,fx1:fx2)=7; %bottom slip
% u_type(fy1:fy2-1,fx1)=8;% left slip
% u_type(fy2,fx1:fx3-1)=9;% top slip
u_type(fy2-1,fx3:fx4)=7;%bottom
u_type(fy2:fy3-1,fx3)=8;%left slip
u_type(fy3,fx3:fx4)=9;%top slip
u_type(fy2:fy3-1,fx4)=11;%right slip
%u_type(fy2,fx4+1:fx2)=9;%top slip
%u_type(fy1:fy2-1,fx2)=11;%right
%u_type(fy1:fy2-1,fx1+1:fx2-1)=10;
u_type(fy2:fy3-1,fx3+1:fx4-1)=14; %remove

% 
%   u_guess(fy2:fy3-1,fx3+1)=0;%left slip
%   u_guess(fy2:fy3-1,fx3+1)=0;
gx1= find(xv>4,1);
gx2=find(xv>10,1);
gx3=find(xv>5,1)-1;
gx4=find(xv>8,1);

gy1=find(yv==1);
gy2=find(yv==3);
gy3=find(yv==5);

v_type=zeros(Nyv,Nxv);

%v_type(gy1,gx1:gx2-1)=7; %bottom %gx
%v_type(gy1:gy2,gx1-1)=8;%left v=0
%v_type(gy2,gx1:gx3)=9;% top
v_type(gy2,gx3+1:gx4-1)=7;%bottom
v_type(gy2:gy3,gx3)=8;%left
v_type(gy3,gx3+1:gx4-1)=9;%top
v_type(gy2:gy3,gx4)=11;%right
% v_type(gy2,gx4:gx2-1)=9;
% v_type(gy1:gy2,gx2)=11;%right
% v_type(gy1+1:gy2-1,gx1:gx2-1)=10; %interior
v_type(gy2+1:gy3-1,gx3+1:gx4-1)=14;


p_type=zeros(Nyp,Nxp);

hx1=find(xCenter>4,1);
hx2=find(xCenter>10,1);
hx3=find(xCenter>5,1);
hx4=find(xCenter>8,1);

hy1=find(yCenter>1,1);
hy2=find(yCenter>3,1);
hy3=find(yCenter>5,1);


% p_type(hy1-1,hx1:hx2-1)=7; %bottom
% p_type(hy1:hy2-1,hx1-1)=8; %left
% p_type(hy2,hx1:hx3-2)=9;%top
p_type(hy2-1,hx3:hx4-1)=7;%bottom
p_type(hy2:hy3-1,hx3-1)=8; %left
% p_type(hy2,hx3-1)=12; %top(s)-left(e)
p_type(hy3,hx3:hx4-1)=9;%top
p_type(hy2:hy3-1,hx4)=11;%right
% p_type(hy2+1:hy3-1,hx4)=11;%right
% p_type(hy2,hx4)=13; %top(s)-right(w)
% p_type(hy2,hx4+1:hx2-1)=9;%top
% p_type(hy1:hy2-1,hx2)=11;%right

%p_type(hy1:hy2-1,hx1:hx2-1)=10;%interior
p_type(hy2:hy3-1,hx3:hx4-1)=14;



T_type=zeros(Nyt,Nxt);
t_type=zeros(Nyt,Nxt);
tx1=find(xt>4,1);
tx2=find(xt>10,1);
tx3=find(xt>5,1);
tx4=find(xt>8,1);

ty1=find(yt>1,1);
ty2=find(yt>3,1);
ty3=find(yt>5,1);

% T_type(ty1-1,tx1:tx2-1)=7; %bottom
% T_type(ty1:ty2-1,tx1-1)=8; %left
% T_type(ty2,tx1:tx3-2)=9;%top
T_type(ty2-1,tx3:tx4-1)=7;
T_type(ty2:ty3-1,tx3-1)=8; %left
% T_type(ty2,tx3-1)=12; %top(s)-left(e)
T_type(ty3,tx3:tx4-1)=9;%top
T_type(ty2:ty3-1,tx4)=11;%right
% T_type(ty2+1:ty3-1,tx4)=11;%right
% T_type(ty2,tx4)=13; %top(s)-right(w)
% T_type(ty2,tx4+1:tx2-1)=9;%top
% T_type(ty1:ty2-1,tx2)=11;%right
t_type(ty2:ty3-1,tx3)=15;
t_type(ty2:ty3-1,tx4-1)=16;
Tc_type=zeros(Nyt,Nxt);
jx=find(xt>6.5,1); jy=find(yt>4,1);
Tc_type(jy,tx3:tx4-1)=17;
Tc_type(ty2:ty3-1,jx)=18;

% T_type(ty1:ty2-1,tx1:tx2-1)=10;%interior
T_type(ty2:ty3-1,tx3:tx4-1)=14;
Tu_type=zeros(Nyt,Nxt);
Tu_type(ty3-1,tx3:tx4-1)=19; %Top border cells

Tb_type=zeros(Nyt,Nxt);
Tb_type(ty2,tx3:tx4-1)=20; %Bottom edged cells

Tl_type=zeros(Nyt,Nxt);
Tl_type(ty2:ty3-1,tx3)=21;

Tr_type=zeros(Nyt,Nxt);
Tr_type(ty2:ty3-1,tx4-1)=22;



 
