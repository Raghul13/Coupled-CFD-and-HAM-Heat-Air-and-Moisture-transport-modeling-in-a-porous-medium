clear; close all; clc;
 
% Extracting the grid details from the mesh.m file
 [xMesh,yMesh,Nxt,Nyt,Nxu,Nyu,Nxp,Nyp,Nxv,Nyv,xu_grid,xv_grid,yu_grid,yv_grid,DX,DY,DXv,DYv,u_type,v_type,p_type,xu,xv,xCenter,yCenter,xp,yp,xt,yt,yu,yv,T_type,t_type,dxt,dyt]=mesh;

% Source Code

%% Solving the X,Y - Momentum and Continuity Equation using the SIMPLE Algorithm
Re = 1;
y_u=yMesh;
x_u=xMesh;
Ny_u=Nyu;
Nx_u=Nxu;
Ny_v=Nyv;
Nx_v=Nxv;
Ny_p=Nyp;
Nx_p=Nxp;

%% Set the relaxation parameters

alpha_u = .1;
alpha_v = .1;
alpha_p = .0005; %.0005
alpha_T=0.002;
%% Set the boundary speeds

u_bot =  0; v_bot = 0;
u_top =  0; v_top = 0;
u_lef = 1 ; v_left = 0; 
u_rig =  0; v_right = 0; 

%% Create the guessed values

u_guess =ones(Ny_u,Nx_u); 
v_guess = zeros(Ny_v,Nx_v); 
p_guess = zeros(Ny_p,Nx_p);
rho=1;

u_guess=u_lef*u_guess;
u_guess(u_type==8)=0;
u_guess(u_type==11)=0;
% u_guess(u_type==10)=0;
u_guess(u_type==14)=0;
 

%% Set the stop conditions

II_max =150;
u_tol = 10^-5;
II = 1;
u_change = 1;

while II <= II_max && abs(u_change(end)) > u_tol

    u_guess_old = u_guess;
    v_guess_old = v_guess;
    p_guess_old = p_guess;    
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculation of u-star
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  



dxp=DX;
dxe=[DX(:,2:end) zeros(Nyu,1)];
dxw=[zeros(Nyu,1) DX(:,1:end-1)];

dyp=DY;
dyn = [DY(2:end,:); zeros(1,Nxu)]; 
dys = [zeros(1,Nxu); DY(1:end-1,:)];


 %% Calculation of the velocity at the midpoints (nonlinear terms)

   % ue = [zeros(Ny_u,1) (u_guess(:,2:end-1)+u_guess(:,3:end))/2 zeros(Ny_u,1)];
   ue = [zeros(Ny_u,1) ((u_guess(:,2:end-1).*DX(:,3:end))+(u_guess(:,3:end).*DX(:,2:end-1))) zeros(Ny_u,1)];
    uw = [zeros(Ny_u,1) ((u_guess(:,1:end-2).*DX(:,2:end-1))+(u_guess(:,2:end-1).*DX(:,1:end-2))) zeros(Ny_u,1)];
    vn = [zeros(Ny_u,1) (v_guess(2:end,2:end-2)+v_guess(2:end,3:end-1))/2 zeros(Ny_u,1)];
    vs = [zeros(Ny_u,1) (v_guess(1:end-1,2:end-2)+v_guess(1:end-1,3:end-1))/2 zeros(Ny_u,1)];
      
    %% Calculate the pressures at the east and west midpoint

    p_e = [zeros(Ny_u,1) p_guess(:,2:end) zeros(Ny_u,1)];
    p_w = [zeros(Ny_u,1) p_guess(:,1:end-1) zeros(Ny_u,1)];
    
    %% Coefficients for the A matrix
    
    Fe_u =rho*( ue./(dxp.*(dxp+dxe).^2));
    Fw_u = rho*(uw./(dxp.*(dxp+dxw).^2)); 
    Fn_u =2*rho*vn./(dyp.*(dyn+dyp));
    Fs_u = 2*rho*vs./(dyp.*(dys+dyp));
  DXe=dxp.*(dxp+dxe);
  DXw=dxp.*(dxp+dxw);
  DYn=dyp.*(dyp+dyn);
  DYs=dyp.*(dys+dyp);
    Deu = 2./(Re*DXe); Dwu = 2./(Re*DXw); Dnu = 2./(Re*DYn); Dsu = 2./(Re*DYs);        

    %% Create the A matrix and the right hand side
    
    A_u = zeros(length(u_guess(:)),length(u_guess(:)));
    b_u = zeros(length(u_guess(:)),1);
    
    
    
    for i = 1:Ny_u*Nx_u
        
        %% If the cell is ON the left
        if i <Ny_u+1 
            A_u(i,i) = 1;   
            b_u(i) = u_lef;
            
        %% If the cell is ON the right boundary
        elseif i > Ny_u*(Nx_u-1)            
            A_u(i,i) = 1; A_u(i,i-Ny_u) = -1;  
            b_u(i) = 0;
            
        %% If the cell is IN the domain
        else
            A_u(i,i)      =  dxe(i)*Fe_u(i)-dxw(i)*Fw_u(i)+dyn(i)*Fn_u(i)-dys(i)*Fs_u(i)+Deu(i)+Dwu(i)+Dnu(i)+Dsu(i);        
            A_u(i,i-Ny_u) = -(dxp(i)*Fw_u(i))-Dwu(i);   
            A_u(i,i+Ny_u) =  dxp(i)*Fe_u(i)-Deu(i);    
            b_u(i) = -(p_e(i)-p_w(i))/dxp(i);

            %% If the cell BOARDERS the top of the domain
            if mod(i,Ny_u) == 0
                A_u(i,i)      = A_u(i,i) - Dnu(i);
                b_u(i) = b_u(i) ;                
            else % If the cell does not boarder the top then 
                 % there is a North cell
                A_u(i,i+1)    =  dyn(i).*Fn_u(i)-Dnu(i); 
            end

            %% If the cell BOARDERS the bottom of the domain
            if mod(i,Ny_u) == 1
                A_u(i,i)      = A_u(i,i) - Dsu(i);    %+Dsu(i)        
                %b_u(i) = b_u(i) + 2*u_bot*Dsu(i);
                b_u(i)=b_u(i);
            else % If the cell does not boarder the bottom of 
                 % the domain then there is a South cell                
                A_u(i,i-1)    = -(dxp(i)*Fs_u(i))-Dsu(i); 
            end

        end
    end
    
   for i = 1:length(u_type(:))
       if u_type(i)== 7 %bottom
            A_u(i,:)=0;
            A_u(i,i) =  dxe(i)*Fe_u(i)-dxw(i)*Fw_u(i)+dyn(i)*Fn_u(i)-dys(i)*Fs_u(i)+Deu(i)+Dwu(i)+2*Dnu(i)+Dsu(i);      
            A_u(i,i-Ny_u) = -(dxp(i)*Fw_u(i))-Dwu(i);   
            A_u(i,i+Ny_u) =  (dxp(i)*Fe_u(i))-Deu(i); 
            A_u(i,i+1)  =  0;
            A_u(i,i-1)    =  -(dyp(i)*Fs_u(i))-Dsu(i);            
            b_u(i)        = -(p_e(i)-p_w(i))/dxp(i) + 2*u_bot*Dnu(i);
        elseif u_type(i)== 8 % left
             A_u(i,i)=1;
              A_u(i,i-Ny_u) = 0;   
               A_u(i,i+Ny_u) =  0;
               A_u(i,i+1)  =  0;
               A_u(i,i-1)    =  0;
             b_u(i)=0;
         elseif u_type(i)==11 % right
             A_u(i,i)=1;
             A_u(i,i-Ny_u) = 0;   
              A_u(i,i+Ny_u) =  0;
              A_u(i,i+1)  =  0;
              A_u(i,i-1)    =  0;
             b_u(i)=0;
         elseif u_type(i)== 9 %top
                A_u(i,i)      =  dxe(i)*Fe_u(i)-dxw(i)*Fw_u(i)+dyn(i)*Fn_u(i)-dys(i)*Fs_u(i)+Deu(i)+Dwu(i)+Dnu(i)+2*Dsu(i);      
                A_u(i,i-Ny_u) = -dxp(i)*Fw_u(i)-Dwu(i);   
                A_u(i,i+Ny_u)=  dxp(i)*Fe_u(i)-Deu(i);
                A_u(i,i-1)    =  0;
                A_u(i,i+1)  =  dyp(i)*Fn_u(i)-Dnu(i);            
                b_u(i)        = -(p_e(i)-p_w(i))/dxp(i) + 2*u_bot*Dsu(i);
       elseif u_type(i)==14 %interior             
             A_u(i,i)=1;
              A_u(i,i-Ny_u) = 0;   
              A_u(i,i+Ny_u) =  0;
              A_u(i,i+1)  =  0;
              A_u(i,i-1)    =  0;
             b_u(i) =0;
      end
           
   end
   
    %% Calculate the u star value    
    u_star = reshape(A_u\b_u,Ny_u,Nx_u);
    
    %% Update the u_star value with relaxation
    u_star = u_guess + alpha_u*(u_star-u_guess);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculation of v-star
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 dxpv=DXv;
dxwv=[zeros(Nyv,1) DXv(:,1:end-1)];
dxev=[DXv(:,2:end) zeros(Nyv,1)];

dypv=DYv;
dysv=[zeros(1,Nxv); DYv(1:end-1,:)];
dynv=[DYv(2:end,:); zeros(1,Nxv)];


    %% Calculation of the velocity at the midpoints (nonlinear terms)

    ue = [zeros(1,Nx_v); zeros(Ny_v-2,1) (u_guess(1:end-1,2:end)+u_guess(2:end,2:end))/2 zeros(Ny_v-2,1); zeros(1,Nx_v)];
    uw = [zeros(1,Nx_v);zeros(Ny_v-2,1) (u_guess(1:end-1,1:end-1)+u_guess(2:end,1:end-1))/2 zeros(Ny_v-2,1); zeros(1,Nx_v)];
    vn = [zeros(1,Nx_v); (v_guess(2:end-1,:).*DYv(3:end,:))+(v_guess(3:end,:).*DYv(2:end-1,:)); zeros(1,Nx_v)];
    vs = [zeros(1,Nx_v); (v_guess(1:end-2,:).*DYv(2:end-1,:)+v_guess(2:end-1,:).*DYv(1:end-2,:)); zeros(1,Nx_v)];

    %% Calculate of pressure at the north and south midpoints
 
    p_n = [zeros(1,Nx_v); zeros(Ny_v-2,1)  p_guess(2:end,:) zeros(Ny_v-2,1)  ; zeros(1,Nx_v)]; %(2:end)
    p_s = [zeros(1,Nx_v);zeros(Ny_v-2,1)  p_guess(1:end-1,:) zeros(Ny_v-2,1); zeros(1,Nx_v)]; %(1:end-1)

    %% Calculation of the coefficients for the matrix

    Fe_v = rho*ue./(dxpv.*(dxpv+dxev)); 
    Fw_v = rho*uw./(dxpv.*(dxpv+dxwv)); 
    Fn_v=rho*vn./(dypv.*(dypv+dynv).^2); 
    Fs_v = rho*vs./(dypv.*(dysv+dypv).^2);
    Dev = 2./(Re*(dxpv.*(dxpv+dxev))); Dwv = 2./(dxpv.*(dxpv+dxwv)*Re); Dnv = 2./(Re*dypv.*(dypv+dynv)); Dsv = 2./(dypv.*(dypv+dysv)) ;

    %% Create the A matrix and b vector for the vertical velocity calculation

    A_v = eye(length(v_guess(:)),length(v_guess(:)));
    b_v = zeros(length(v_guess(:)),1);

    %% Set the values of the matrix

    for i = 1:Ny_v*Nx_v

        %% If the cell is ON the top boundary    
        if mod(i,Ny_v) == 0 
            A_v(i,i) = 1;    
            b_v(i) = v_top;

        %% If the cell is ON the bottom boundary
        elseif mod(i,Ny_v) == 1          
            A_v(i,i) = 1; 
            b_v(i) = v_bot;
            
        %% If the cell is a Left Ghost cell
        elseif i < Ny_v
            A_v(i,i) = 1;
            A_v(i,i+Ny_v) = 1;
            b_v(i) = 2*v_left;
            
        %% If the cell is a Right Ghost cell
        elseif i > Ny_v*(Nx_v-1)
            A_v(i,i) = 1;
            A_v(i,i-Ny_v) = -1;
            b_v(i) = 0;

        %% If the cell is In the domain
        else
            A_v(i,i)      =  dxev(i)*Fe_v(i)-dxwv(i)*Fw_v(i)+dynv(i)*Fn_v(i)-dysv(i)*Fs_v(i)+Dev(i)+Dwv(i)+Dnv(i)+Dsv(i);
            A_v(i,i+1)    =  dypv(i)*Fn_v(i)-Dnv(i);  
            A_v(i,i-1)    = -(dypv(i)*Fs_v(i))-Dsv(i); 
            A_v(i,i-Ny_v) = -(dxpv(i)*Fw_v(i))-Dwv(i);             
            A_v(i,i+Ny_v) =  dxpv(i)*Fe_v(i)-Dev(i);
            b_v(i) = -(p_n(i)-p_s(i))/dypv(i);              
        end
    end
%     
       for i = 1:length(v_type(:))
          
          if v_type(i)== 7 %bottom
              A_v(i,:) = 0;
                A_v(i,i)=1;
                A_v(i,i-Ny_v) = 0;   
                A_v(i,i+Ny_v) =  0;
                A_v(i,i+1)  =  0;
                A_v(i,i-1)    =  0;
                b_v(i)=0;
         elseif v_type(i)== 8 %left
             A_v(i,i)      =  -dxwv(i)*Fw_v(i)+dynv(i)*Fn_v(i)-dysv(i)*Fs_v(i)+2*Dev(i)+Dwv(i)+Dnv(i)+Dsv(i);
             A_v(i,i+1)    =  dypv(i)*Fn_v(i)-Dnv(i);  
             A_v(i,i-1)    = -(dypv(i)*Fs_v(i))-Dsv(i);
                  A_v(i,i+Ny_v) = 0; 
             A_v(i,i-Ny_v) = -(dxpv(i)*Fw_v(i))-Dwv(i);
             b_v(i)= -(p_n(i)-p_s(i))/dypv(i);
         
         elseif v_type(i)==11 %right
             A_v(i,i)      =  dxev(i)*Fe_v(i)+dynv(i)*Fn_v(i)-dysv(i)*Fs_v(i)+Dev(i)+2*Dwv(i)+Dnv(i)+Dsv(i);
             A_v(i,i+1)    =  dypv(i)*Fn_v(i)-Dnv(i);  
             A_v(i,i-1)    = -(dypv(i)*Fs_v(i))-Dsv(i);
             A_v(i,i-Ny_v) = 0; 
             A_v(i,i+Ny_v) =  dxpv(i)*Fe_v(i)-Dev(i);
             b_v(i)= -(p_n(i)-p_s(i))/dypv(i);
  
         elseif v_type(i)== 9 %top
              A_v(i,i)= 1;
              A_v(i,i-Ny_v) = 0;   
             A_v(i,i+Ny_v) =  0;
             A_v(i,i+1)  =  0;
             A_v(i,i-1)    =  0;
              b_v(i)=0;
         elseif v_type(i)==14%interior
             
             A_v(i,i)=1;
             A_v(i,i-Ny_v) = 0;   
             A_v(i,i+Ny_v) =  0;
             A_v(i,i+1)  =  0;
             A_v(i,i-1)    =  0;
             b_v(i) =0;
          end
       end
         

A_v=sparse(A_v);
    %% Calculate the v-star value

    v_star = reshape(A_v\b_v,Ny_v,Nx_v);

    %% Use relaxation to set the v-star value

    v_star = v_guess + alpha_v*(v_star-v_guess);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculation of the pressure correction
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Get the velocity at the midpoints of the Pressure cell

    uw_star = u_star(:,1:end-1);
    ue_star = u_star(:,2:end);
    vn_star = v_star(2:end,2:end-1);
    vs_star = v_star(1:end-1,2:end-1);

%% Set the value of aP for the u and v velocities

    Ap_u = reshape(diag(A_u),Ny_u,Nx_u);
    Ap_v = reshape(diag(A_v),Ny_v,Nx_v);
   Ap_v = Ap_v(:,2:end-1);

%% Calculate the pressure correction matrix coefficients

[DXpress,DYpress]=meshgrid(abs(diff(xMesh)),abs(diff(yMesh)));

dx_p = DXpress;
dx_e = [DXpress(:,2:end) zeros(Ny_p,1)]; 
dx_w = [zeros(Ny_p,1) DXpress(:,1:end-1)]; 
dy_p = DYpress;
dy_n = [DYpress(2:end,:); zeros(1,Nx_p)]; 
dy_s = [zeros(1,Nx_p); DYpress(1:end-1,:)];

    Cw = 2./(Ap_u(:,1:end-1).*(dx_p.*(dx_p+dx_w)));
    Ce = 2./(Ap_u(:,2:end).*(dx_p.*(dx_p+dx_e)));
    Cn = 2./(Ap_v(2:end,:).*(dx_p.*(dy_p+dy_n)));
    Cs = 2./(Ap_v(1:end-1,:).*(dx_p.*(dy_p+dy_s)));

%% Create the pressure correction matrix

    A_p = zeros(length(p_guess(:)));
    b_p = zeros(length(p_guess(:)),1);

%% Set the pressure correction coefficients

    for i = 1:Ny_p*Nx_p       

        % Not on left boundary
        if i>Ny_p 
            A_p(i,i) = A_p(i,i) + Cw(i);
            A_p(i,i-Ny_p) = -Cw(i);
            b_p(i) = -(ue_star(i)-uw_star(i))/dx_p(i) -(vn_star(i)-vs_star(i))/dy_p(i);
       
            %b_p(i) = -(ue_star(i)-uw_star(i))/dx_p(i) -(vn_star(i)-vs_star(i))/dy_p(i);
        end
          

        % Not on right boundary
        if i<=Ny_p*(Nx_p-1)
            A_p(i,i) = A_p(i,i) + Ce(i);
            A_p(i,i+Ny_p) = -Ce(i);
             b_p(i) = -(ue_star(i)-uw_star(i))/dx_p(i) -(vn_star(i)-vs_star(i))/dy_p(i);
           
        end

        % Not at the top boundary
        if mod(i,Ny_p)~=0
            A_p(i,i) = A_p(i,i) + Cn(i);
            A_p(i,i+1) = -Cn(i);
             b_p(i) = -(ue_star(i)-uw_star(i))/dx_p(i) -(vn_star(i)-vs_star(i))/dy_p(i);
           
        end    

        % Not at the bottom boundary
        if mod(i,Ny_p)~=1
            A_p(i,i) = A_p(i,i) + Cs(i);
            A_p(i,i-1) = -Cs(i);
             b_p(i) = -(ue_star(i)-uw_star(i))/dx_p(i) -(vn_star(i)-vs_star(i))/dy_p(i);
           
        end

       % b_p(i) = -(ue_star(i)-uw_star(i))/dx_p(i) -(vn_star(i)-vs_star(i))/dy_p(i);
%         
        if p_type(i)== 7 %bottom
             A_p(i,i) = Cs(i)+Ce(i)+Cw(i) ;
             A_p(i,i+Ny_p) = -Ce(i);
             A_p(i,i+1) = 0;
            A_p(i,i-1) = -Cs(i);
            A_p(i,i-Ny_p) = -Cw(i);
            b_p(i) = -(ue_star(i)-uw_star(i))/dx_p(i) -(-vs_star(i))/dy_p(i);
          
       
        elseif p_type(i)== 8 %left
            A_p(i,i) = Cw(i)+Cn(i)+Cs(i);
         A_p(i,i+Ny_p) = 0;
        A_p(i,i+1) = -Cn(i);
        A_p(i,i-1) = -Cs(i);
        A_p(i,i-Ny_p) = -Cw(i); 
         b_p(i) = -(-uw_star(i))/dx_p(i) -(vn_star(i)-vs_star(i))/dy_p(i);
        
        elseif p_type(i)==11 %right
            A_p(i,i) = Ce(i)+Cn(i)+Cs(i);
        A_p(i,i+Ny_p) = -Ce(i);
        A_p(i,i+1) = -Cn(i);
        A_p(i,i-1) = -Cs(i);
        A_p(i,i-Ny_p) = 0; 
         b_p(i) = -(ue_star(i))/dx_p(i) -(vn_star(i)-vs_star(i))/dy_p(i);
         
        elseif p_type(i)==9  %top
            A_p(i,i) = Cn(i)+Ce(i)+Cw(i);
        A_p(i,i+Ny_p) = -Ce(i);
        A_p(i,i+1) = -Cn(i);
        A_p(i,i-1) = 0;
        A_p(i,i-Ny_p) = -Cw(i); 
         b_p(i) = -(ue_star(i)-uw_star(i))/dx_p(i) -(vn_star(i))/dy_p(i);
         
        elseif p_type(i)==14 %interior
             A_p(i,i) = 1;
              A_p(i,i+Ny_p) = 0;
        A_p(i,i+1) = 0;
        A_p(i,i-1) = 0;
        A_p(i,i-Ny_p) = 0; 
               b_p(i)=0;
        end
       
    end

%% Set the correction of the pressure at the middle of the domain to zero

%        ind = floor((Ny_p*Nx_p)-Nx_p/3- 3*Ny_p); 
%  ind = floor((Ny_p*Nx_p)- 3*Ny_p); 
   ind = floor((Ny_p*Nx_p)-Nx_p/2 - 3*Ny_p); 

     A_p(ind,:) = 0;
     A_p(ind,ind) = 1;
     b_p(ind) = 0;

%% Calculate the pressure correction and apply it with under-relaxation

    p_correction = reshape(A_p\b_p,Ny_p,Nx_p);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculation of the velocity correction
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculate the pressure correction on the east and west faces of u-cells

    pc_E = p_correction(:,2:end);
    pc_W = p_correction(:,1:end-1);

%% Calculate the velocity correction for all cells IN the domain

    u_correction = 0*u_guess;
    u_correction(:,2:end-1) = 2*(pc_W-pc_E)./((dx_p(:,1:end-1)+dx_p(:,2:end)).*Ap_u(:,2:end-1));
%     u_correction(u_type==10)=0;
     u_correction(u_type==14)=0;
     u_correction(u_type==11)=0;
     u_correction(u_type==8)=0;
     
%% Calculate the pressure correction on the north and south faces of v-cells

    pc_N = p_correction(2:end,:);
    pc_S = p_correction(1:end-1,:);

%% Calculate the velocity correction for all cells IN the domain

    v_correction = 0*v_guess;
    v_correction(2:end-1,2:end-1) = 2*(pc_S-pc_N)./((dy_p(1:end-1,:)+dy_p(2:end,:)).*Ap_v(2:end-1,:));
%     v_correction(v_type==10)=0;
    v_correction(v_type==14)=0;
%     v_correction(v_type==8)=0;
%     v_correction(v_type==11)=0;
    v_correction(v_type==9)=0;
    v_correction(v_type==7)=0;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Corect the velocity directly and pressure with under-relaxation
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    u_guess = u_star + u_correction;
    v_guess = v_star + v_correction;

    p_guess = p_guess + alpha_p*p_correction;

%% Impose flux conservation

    m_in = sum(u_guess(:,1).*dyp(:,1));
    m_out = sum(u_guess(:,end-1).*dyp(:,end-1));

    u_guess(:,end) = m_in/m_out*u_guess(:,end-1);
    
ue = u_guess(:,2:end);
    uw = u_guess(:,1:end-1);
    vn = v_guess(2:end,2:end-1);
  % vn = v_guess(2:end,1:end-1);
    vs = v_guess(1:end-1,2:end-1);  

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot the results
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    if mod(II,1) == 0
    figure(2);
    hold off; 
    u_temp = [u_bot+zeros(1,Nx_u); u_guess; u_top+zeros(1,Nx_u)];
    v_temp = [v_left+zeros(Ny_v,1)  v_guess  v_right+zeros(Ny_v,1)];

    x_vt = [0 xv 1]; [X_v, Y_v] = meshgrid(x_vt,yv);
    y_ut = [0 yu 1]; [X_u, Y_u] = meshgrid(xu,y_ut);
   x_min=xu(1);
   x_max=xu(end);
   y_min=yu(1);
   y_max=yu(end);
    [X,Y] = meshgrid(linspace(x_min,x_max,100),linspace(y_min,y_max,40));

    U = griddata(X_u(:),Y_u(:),u_temp(:),X,Y);
    V = griddata(X_v(:),Y_v(:),v_temp(:),X,Y);

    [x_pt,y_pt] = meshgrid(xu,yu);
    fv_plotting(x_pt,y_pt,reshape(u_guess,Ny_u,Nx_u)); hold on;
    quiver(X,Y,U,V,'k')
    colorbar
    
    end
    
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate the continuity residule and changes in u,v, and p
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ue = u_guess(:,2:end);
    uw = u_guess(:,1:end-1);
    vn = v_guess(2:end,2:end-1);
  % vn = v_guess(2:end,1:end-1);
    vs = v_guess(1:end-1,2:end-1);
    % vs = v_guess(1:end-1,1:end-1);

     
     continuity_residule(II) = sqrt(sum(((ue(:)-uw(:))./dx_p(:)+(vn(:)-vs(:))./dy_p(:)).^2));

    u_change(II) = sqrt(sum((u_guess(:)-u_guess_old(:)).^2))/Nx_u/Ny_u;
    v_change(II) = sqrt(sum((v_guess(:)-v_guess_old(:)).^2))/Nx_v/Ny_v;
    p_change(II) = sqrt(sum((p_guess(:)-p_guess_old(:)).^2))/Nx_p/Ny_p;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot the continuity residule
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    figure(1);
    hold off; 
    plot(u_change,'.-r');
    hold on;
    plot(v_change,'.-g');
    plot(p_change,'.-k');
    plot(continuity_residule,'.-b');
    set(gca,'YScale','log')
    drawnow;

    figure(5)
    subplot(3,2,1);
imagesc(u_type); set(gca,'YDir','normal'); colorbar
subplot(3,2,3);
imagesc(u_star); set(gca,'YDir','normal'); colorbar
subplot(3,2,5);
imagesc(u_correction); set(gca,'YDir','normal'); colorbar
subplot(3,2,2);
imagesc(v_type); set(gca,'YDir','normal'); colorbar
subplot(3,2,4);
imagesc(v_star); set(gca,'YDir','normal'); colorbar
subplot(3,2,6);
imagesc(v_correction); set(gca,'YDir','normal'); colorbar
drawnow

% T=At\bt;
% T=T+alpha_T.*T;
%     figure(6);
%        fv_plotting(Xt,Yt,reshape(T,[Nyt Nxt]));
%        
%        colormap jet;
%        drawnow;
    
    
    II = II + 1;
%         figure(6);
%         fv_plotting(Xt,Yt,reshape(T,[Nyt Nxt]));hold on;
%         quiver(X,Y,U,V,'k')
%         figure(7);
%         [c,h]=contour(reshape(T,[Nyt Nxt]),14,'k');
%        clabel(c,h);
% 
%     figure(8);
%    [C,H]=contour(u_guess,14,'k');
%    clabel(C,H);
end

% HEAT TRANSFER PART

uet= [zeros(Nyt,1) ue zeros(Nyt,1)];
 uwt= [zeros(Nyt,1) uw zeros(Nyt,1)];
 vnt=[zeros(Nyt,1) vn zeros(Nyt,1)];
 vst=[zeros(Nyt,1) vs zeros(Nyt,1)];
 k=0.02; % Thermal Conductivity of air
 
 [Xt,Yt]=meshgrid(xt,yt);
 %t_type=p_type;
 dt=0.5;
 t=0:dt:100;
nt=length(t);
%   dxt= 0.1*ones(Nyt,Nxt); 
%    dyt=0.1*ones(Nyt,Nxt);
dx_pt = dxt;
dy_pt = dyt;
dx_et = [dxt(:,2:end) zeros(Nyt,1)]; 
dx_wt = [zeros(Nyt,1) dxt(:,1:end-1)]; 
%dy_pt = dyt;
dy_nt = [dyt(2:end,:); zeros(1,Nxt)]; 
dy_st = [zeros(1,Nxt); dyt(1:end-1,:)];
 bt=zeros(Nyt*Nxt,1);
% Td=zeros(Nyt*Nxt,1);
Fe = uet./(dx_pt.*(dx_pt+dx_et)); 
Fw = uwt./(dx_pt.*(dx_pt+dx_wt)); 
Fn = vnt./(dy_pt.*(dy_pt+dy_nt)); 
Fs = vst./(dy_pt.*(dy_pt+dy_st));

 
De=2*k./(dx_pt.*(dx_pt+dx_et));
Dw=2*k./(dx_pt.*(dx_pt+dx_wt)); 
Dn=2*k./(dy_pt.*(dy_pt+dy_nt));
Ds=2*k./(dy_pt.*(dy_pt+dy_st));

 T_left=60;
 
 
 
 T=25*ones(Nyt,Nxt);
 
 
%T=5*T;
 %[a,b]=size(T(T_type==10));
% [a1,b1]=size(T(T_type==14));
%  
  %T(T_type==10)=5*ones(a,b);
  %T(T_type==14)=5*ones(a1,b1);
 T=reshape(T,[Nyt*Nxt 1]);

% Thermal properties of Eggplant
ks=0.5;
%  ks=0.2;
 
 rhos=870;
 cps=1800;
 rhol=750;
 cpl=750;
 kl=0.15;
% h=4.026;
 rs=1/rhos*cps;
De1=(2*ks*rs./(dx_pt.*(dx_pt+dx_et)));
Dw1=(2*ks*rs./(dx_pt.*(dx_pt+dx_et)));
Dn1=(2*ks*rs./(dy_pt.*(dy_pt+dy_nt)));
Ds1=(2*ks*rs./(dy_pt.*(dy_pt+dy_st)));

% ksf= ks*k/(k+ks);
ksf=(k+ks)/2;
% ksf=k;
Def=(2*ksf*rs./(dx_pt.*(dx_pt+dx_et)));
Dwf=(2*ksf*rs./(dx_pt.*(dx_pt+dx_et)));
Dnf=(2*ksf*rs./(dy_pt.*(dy_pt+dy_nt)));
Dsf=(2*ksf*rs./(dy_pt.*(dy_pt+dy_st)));
tc_type=zeros(Nyt,Nxt);
jx=find(xt>6.5,1); jy=find(yt>4,1);
tc_type(jy,jx)=1; %Center

Tbot=100;
 %bt=T;
 bt=zeros(Nyt*Nxt,1);
 Tt=zeros(nt,Nyt*Nxt);
 tic 
 %% Time Stepping Loop
   for n=2:nt
  At=zeros(Nxt*Nyt);
  bt=T/dt;
   %% Setting the Coeff A
   for i = 1:Nyt*Nxt
        
        %% If the cell is ON the left
        if i <Nyt+1 
            At(i,i) = 1;
            At(i,i+Nxt)=1;
            
            bt(i)=2*T_left;
        %% If the cell is ON the right boundary
        elseif i > Nyt*(Nxt-1)            
            At(i,i) = 1; 
            At(i,i-Nyt) = -1;  
            bt(i) = 0;
            
        %% If the cell is IN the domain
        else
            At(i,i)      = (dx_et(i)*Fe(i)-dx_wt(i)*Fw(i)+dy_nt(i)*Fn(i)-dy_st(i)*Fs(i))...
                           +(De(i)+Dw(i)+Dn(i)+Ds(i)+1/dt);        
            At(i,i-Nyt) = -(dx_pt(i)*Fw(i))-Dw(i);   
            At(i,i+Nyt) =  dx_pt(i)*Fe(i)-De(i);    
            bt(i)=bt(i);

            %% If the cell BOARDERS the top of the domain
            if mod(i,Nyt) == 0  % Free Slip
                At(i,i)      = At(i,i) - Dn(i);
                 bt(i) = bt(i);                
            else % If the cell does not boarder the top then 
                 % there is a North cell
                At(i,i+1)    =  dy_pt(i)*Fn(i)-Dn(i); 
            end

            %% If the cell BOARDERS the bottom of the domain
            if mod(i,Nyt) == 1 % Free Slip
%                 At(i,i)      = At(i,i) + Ds(i);      no slip
                  At(i,i)=At(i,i)-Ds(i);
               % bt(i) = bt(i)+2*Tbot*Ds(i);
               bt(i) = bt(i);
            else % If the cell does not boarder the bottom of 
                 % the domain then there is a South cell                
                At(i,i-1)    = -(dy_pt(i)*Fs(i))-Ds(i); 
            end

        end
   end
   
   %SOLID PART
      for i = 1:length(T_type(:))
            if T_type(i)== 7  %Bottom
              At(i,:)=0;
               At(i,i)      = (dx_et(i)*Fe(i)-dx_wt(i)*Fw(i)+dy_nt(i)*Fn(i)-dy_st(i)*Fs(i))...
                               +(Def(i)+Dwf(i)+Dn1(i)+Ds(i)+1/dt);     % Present Cell
               At(i,i-Nyt) = -(dx_pt(i)*Fw(i))-Dwf(i);                 % West Cell
            At(i,i+Nyt) =  dx_pt(i)*Fe(i)-Def(i);                      % East Cell
            At(i,i-1)=-(dy_pt(i)*Fs(i))-Ds(i);                         % South Cell
             At(i,i+1)=-Dn1(i); %Solid                                 % North Cell
            bt(i) = bt(i);
%               At(i,i)=1;
%               At(i,i+1)=1;
%               bt(i)= 60;
            %Td(i)=bt(i);
%                 At(i,i) = 1;      
%          At(i,i-Nyt) = 0;   
%          At(i,i+Nyt) =  0;    
%          At(i,i-1)    =  0; 
%          At(i,i+1)  =  -1;
 
        
        
         elseif T_type(i)== 8 % Left
              At(i,:)=0;
               At(i,i)      =(dx_et(i)*Fe(i)-dx_wt(i)*Fw(i)+dy_nt(i)*Fn(i)-dy_st(i)*Fs(i)+De1(i)+Dw(i)+Dnf(i)+Dsf(i)+1/dt);  
               At(i,i-Nyt) = -(dx_pt(i)*Fw(i))-Dw(i);   
            At(i,i+Nyt) = -De1(i); % solid   
            At(i,i+1)= dy_pt(i)*Fn(i)-Dnf(i);
            At(i,i-1)=-(dy_pt(i)*Fs(i))-Dsf(i);
            At(i,i+Nyt)=-De1(i); % Solid
            bt(i)=bt(i);
%             At(i,i)=1;
%             At(i,i+Nyt)=1;
%                 bt(i) = 60;
            
         elseif T_type(i)==11  % Right
              At(i,:)=0;
               At(i,i)      =(dx_et(i)*Fe(i)-dx_wt(i)*Fw(i)+dy_nt(i)*Fn(i)-dy_st(i)*Fs(i)+De(i)+Dw1(i)+Dnf(i)+Dsf(i)+1/dt);
              % At(i,i-Nyt) = -Dwf(i);   
            At(i,i+Nyt) =  dx_pt(i)*Fe(i)-De(i);    
            At(i,i+1)= dy_pt(i)*Fn(i)-Dnf(i);
            At(i,i-1)=-(dy_pt(i)*Fs(i))-Dsf(i);
            At(i,i-Nyt)=-Dw1(i);

%                At(i,i)=1;
%                At(i,i-Nyt)=1;
                bt(i) = bt(i);
                %Td(i)=bt(i);
%              At(i,i)=1;
%              At(i,i-Nyt) = -1;   
%                At(i,i+Nyt) =  0;
%                At(i,i+1)  =  0;
%                At(i,i-1)    =  0;
%               bt(i)=0;
         elseif T_type(i)== 9  % Top
              At(i,:)=0;
               At(i,i)      = (dx_et(i)*Fe(i)-dx_wt(i)*Fw(i)+dy_nt(i)*Fn(i)-dy_st(i)*Fs(i)+Def(i)+Dwf(i)+Dn(i)+Ds1(i)+1/dt);  
               At(i,i-Nyt) = -(dx_pt(i)*Fw(i))-Dwf(i); 
               At(i,i+1)= dy_pt(i)*Fn(i)-Dn(i);
               At(i,i-1)=-Ds1(i);
            At(i,i+Nyt) =  dx_pt(i)*Fe(i)-Def(i);

%                At(i,i)=1;
%                At(i,i-1)=1;
                bt(i) = bt(i);   
%                 
                
                
                
            elseif T_type(i)==14 %interior-solid
               At(i,:)=0;
%              At(i,i)=1;
%                 At(i,i-Nyt) = 0;   
%                 At(i,i+Nyt) =  0;
%                 At(i,i+1)  =  0;
%                 At(i,i-1)    =  0;
%               bt(i) =0;

            At(i,i)      = (De1(i)+Dw1(i)+Dn1(i)+Ds1(i)+1/dt);        
            At(i,i-Nyt) = -Dw1(i);   
            At(i,i+Nyt) = -De1(i); 
             At(i,i-1)=-Ds1(i);
              At(i,i+1)= -Dn1(i);
%               bt(i)=bt(i);
              
     %% PHASE CHANGE MODEL
           %Ts=24;
              %Tl=25.1;
              
%               if Tt(n-1,i)>=Ts && Tt(n-1,i)<=Tl
%               % Calculation of evaporation rate for source term
%                   m(n,i) = (Tt(n-1,i)-Ts)/(Tl-Ts);
% %                   mt(n,i) = (m(n,i)/dt)*0.09;
%                    mt(n,i) = (m(n,i)/dt);
%                    G= -1/(Tl-Ts);
%                    kp = m(n,i)*rhol*kl+(1-m(n,i))*rhos*ks;
%                    cpp =m(n,i)*rhol*cpl+(1-m(n,i))*rhos*cps;
%                    rsp= kp/cpp;
%                    rhop=m(n,i)*rhol+(1-m(n,i))*rhos;
%                    Dep=(2*rsp./(dx_pt.*(dx_pt+dx_et)));
%                    Dwp=(2*rsp./(dx_pt.*(dx_pt+dx_et)));
%                    Dnp=(2*rsp./(dy_pt.*(dy_pt+dy_nt)));
%                    Dsp=(2*rsp./(dy_pt.*(dy_pt+dy_st)));
%                      At(i,:)=0;
%                     At(i,i)      = (Dep(i)+Dwp(i)+Dnp(i)+Dsp(i)+1/dt);        
%             At(i,i-Nyt) = -Dwp(i);   
%             At(i,i+Nyt) = -Dep(i); 
%              At(i,i-1)=-Dsp(i);
%               At(i,i+1)= -Dnp(i);
% %                   dG(n,i)= G*(Tt(n-1,i)-Tt(n-2,i))/dt;
%                   T_S(n,i) = (mt(n,i)*90000)/rhop;
% %                     cp= 2400;
%                   bt(i)=bt(i)+(T_S(n,i)/cpp);
% %                      bc(i)=bc(i)+(dG(n,i)*0.05);
%               else

                   bt(i)=bt(i);
            
            end
      end
     toc 
      
     At=sparse(At);
      
      %bt_prev=T;
     tic 
      % Solving T by inverting the Coefficient Matrix At
       T=At\bt;
     toc
      %T=T+alpha_T*T; 
      Tt(n,:)=T';
%       theta_t(n,:)=abs((T-25)./(T_left-25));
      
      
%       del_x=reshape(dx_pt,[Nxt*Nyt,1]);
      %del_t=T(T_type==8)-T(t_type==15);
%       ht(:,n)= ks*((T_left-T(t_type==15))./del_x(T_type==8));
      T_Center(n,1)=T(tc_type==1);
      %SURFACE TEMPERATURE
%       T_surf(n,:)=T(T_type==9);
      
       figure(6);
       fv_plotting(Xt,Yt,reshape(T,[Nyt Nxt]));hold on;
       quiver(X,Y,U,V,'k')
       figure(7);
       [c,h]=contour(reshape(T,[Nyt Nxt]),14,'k');
       clabel(c,h);
        %quiver(X,Y,U,V,'k')
%       bt=T;
       colormap jet;
       drawnow;
   end
   figure(8);
   [C,H]=contour(u_guess,14,'k');
   clabel(C,H);
   
   
   
   
   
%    % SPECIES TRANSPORT
%    
%      D=1.322*10^3;
% %  
%      Tt=273+Tt;
%   dt=2;
%   t=0:dt:90;
%  nt=length(t);
% 
%   bc=zeros(Nyt*Nxt,1);
%  Fe = uet./(dx_pt.*(dx_pt+dx_et)); 
%  Fw = uwt./(dx_pt.*(dx_pt+dx_wt)); 
%  Fn = vnt./(dy_pt.*(dy_pt+dy_nt)); 
%  Fs = vst./(dy_pt.*(dy_pt+dy_st));
%  k=0.016;
%  T=reshape(T,[Nyt Nxt]);
%  De=2*k./(dx_pt.*(dx_pt+dx_et));
%  Dw=2*k./(dx_pt.*(dx_pt+dx_wt)); 
%  Dn=2*k./(dy_pt.*(dy_pt+dy_nt));
%  Ds=2*k./(dy_pt.*(dy_pt+dy_st));
% % 
%   M_left=0.021;
% %  x
% %  
% 
%   M=0.3*ones(Nyt,Nxt);
%   M = reshape(M,[Nyt*Nxt 1]);
% %  
% for i = 1:length(T_type(:))
%             if T_type(i)== 14
%                 M(i)=0;
%                 M(i)=0.92;
%             end
% end
%                 
%    M_in=M;
% for n=2:nt
% % % Thermal properties of Eggplant
% % %ks=0.0998;
%   ks=D*exp(-2830./Tt(n,:));
%  ks=reshape(ks,[Nyt Nxt]);
% % % rhos=390;
% % % cps=2.953;
% % % h=4.026;
% % % rs=1;
%  De1=(2*ks./(dx_pt.*(dx_pt+dx_et)));
%  Dw1=(2*ks./(dx_pt.*(dx_pt+dx_et)));
%  Dn1=(2*ks./(dy_pt.*(dy_pt+dy_nt)));
%  Ds1=(2*ks./(dy_pt.*(dy_pt+dy_st)));
% % 
% ksf= ks*k./(k+ ks);
% % %ksf=(k+ks)/2;
% % ksf=k;
% Def=(2*ksf./(dx_pt.*(dx_pt+dx_et)));
%  Dwf=(2*ksf./(dx_pt.*(dx_pt+dx_et)));
%  Dnf=(2*ksf./(dy_pt.*(dy_pt+dy_nt)));
%  Dsf=(2*ksf./(dy_pt.*(dy_pt+dy_st)));
% %  %bt=T;
% % % bt_prev=zeros(Nyt*Nxt,1);
%     
%    Ac=zeros(Nxt*Nyt);
%    bc=M/dt;
%     for i = 1:Nyt*Nxt
%          
%         %% If the cell is ON the left
%          if i <Nyt+1 
%              Ac(i,i) = 1;
%              Ac(i,i+Nxt)=1;
% %             
%              bc(i)=2*M_left;
%         %% If the cell is ON the right boundary
%          elseif i > Nyt*(Nxt-1)            
%              Ac(i,i) = 1; 
%              Ac(i,i-Nyt) = -1;  
%              bc(i) = 0;
% %             
%          %% If the cell is IN the domain
%          else
%              Ac(i,i)      = dx_et(i)*Fe(i)-dx_wt(i)*Fw(i)+dy_nt(i)*Fn(i)-dy_st(i)*Fs(i)+De(i)+Dw(i)+Dn(i)+Ds(i)+1/dt;        
%             Ac(i,i-Nyt) = -(dx_pt(i)*Fw(i))-Dw(i);   
%             Ac(i,i+Nyt) =  dx_pt(i)*Fe(i)-De(i);    
%             bc(i)=bc(i);
% % 
% %             %% If the cell BOARDERS the top of the domain
%             if mod(i,Nyt) == 0  % Free Slip
%                 Ac(i,i)      = Ac(i,i) - Dn(i);
%                  bc(i) = bc(i);                
%             else % If the cell does not boarder the top then 
% %                  % there is a North cell
%                 Ac(i,i+1)    =  dy_pt(i)*Fn(i)-Dn(i); 
%             end
% % 
% %             %% If the cell BOARDERS the bottom of the domain
%             if mod(i,Nyt) == 1 % Free Slip
% %                 Ac(i,i)      = Ac(i,i) + Ds(i);      no slip
%                   Ac(i,i)=Ac(i,i)-Ds(i);
% %                % bt(i) = bt(i)+2*Tbot*Ds(i);
%                bc(i) = bc(i);
%             else % If the cell does not boarder the bottom of 
% %                  % the domain then there is a South cell                
%                 Ac(i,i-1)    = -(dy_pt(i)*Fs(i))-Ds(i); 
%             end
% % 
%         end
%    end
% %    
% %    %SOLID PART
%       for i = 1:length(T_type(:))
%             if T_type(i)== 7  %Bottom
%               Ac(i,:)=0;
%                Ac(i,i)      = dx_et(i)*Fe(i)-dx_wt(i)*Fw(i)+dy_nt(i)*Fn(i)-dy_st(i)*Fs(i)+Def(i)+Dwf(i)+Dn1(i)+Ds(i)+1/dt;  
%                Ac(i,i-Nyt) = -(dx_pt(i)*Fw(i))-Dwf(i);   
%             Ac(i,i+Nyt) =  dx_pt(i)*Fe(i)-Def(i); 
%             Ac(i,i-1)=-(dy_pt(i)*Fs(i))-Ds(i);
%              Ac(i,i+1)=-Dn1(i); %Solid
%             bc(i) = bc(i);
% %             %Td(i)=bt(i);
% % %                 At(i,i) = 1;      
% % %          At(i,i-Nyt) = 0;   
% % %          At(i,i+Nyt) =  0;    
% % %          At(i,i-1)    =  0; 
% % %          At(i,i+1)  =  -1;
% %  
% %         
% %         
%          elseif T_type(i)== 8 % Left
%               Ac(i,:)=0;
%                Ac(i,i)      =(dx_et(i)*Fe(i)-dx_wt(i)*Fw(i)+dy_nt(i)*Fn(i)-dy_st(i)*Fs(i)+De1(i)+Dw(i)+Dnf(i)+Dsf(i)+1/dt);  
%                Ac(i,i-Nyt) = -(dx_pt(i)*Fw(i))-Dw(i);   
%             Ac(i,i+Nyt) = -De1(i); % solid   
%             Ac(i,i+1)= dy_pt(i)*Fn(i)-Dnf(i);
%             Ac(i,i-1)=-(dy_pt(i)*Fs(i))-Dsf(i);
% %             %At(i,i+Nyt)=-De1(i); % Solid
%                 bc(i) = bc(i);
% %             
%          elseif T_type(i)==11  % Right
%               Ac(i,:)=0;
%                Ac(i,i)      =(dx_et(i)*Fe(i)-dx_wt(i)*Fw(i)+dy_nt(i)*Fn(i)-dy_st(i)*Fs(i)+De(i)+Dw1(i)+Dnf(i)+Dsf(i)+1/dt);
%               % At(i,i-Nyt) = -Dwf(i);   
%             Ac(i,i+Nyt) =  dx_pt(i)*Fe(i)-De(i);    
%             Ac(i,i+1)= dy_pt(i)*Fn(i)-Dnf(i);
%             Ac(i,i-1)=-(dy_pt(i)*Fs(i))-Dsf(i);
%             Ac(i,i-Nyt)=-Dw1(i);
%                 bc(i) = bc(i);
% %                 %Td(i)=bt(i);
% % %              At(i,i)=1;
% % %              At(i,i-Nyt) = -1;   
% % %                At(i,i+Nyt) =  0;
% % %                At(i,i+1)  =  0;
% % %                At(i,i-1)    =  0;
% % %               bt(i)=0;
%          elseif T_type(i)== 9  % Top
%               Ac(i,:)=0;
%                Ac(i,i)      = (dx_et(i)*Fe(i)-dx_wt(i)*Fw(i)+dy_nt(i)*Fn(i)-dy_st(i)*Fs(i)+Def(i)+Dwf(i)+Dn(i)+Ds1(i)+1/dt);  
%                Ac(i,i-Nyt) = -(dx_pt(i)*Fw(i))-Dwf(i); 
%                Ac(i,i+1)= dy_pt(i)*Fn(i)-Dn(i);
%                Ac(i,i-1)=-Ds1(i);
%             Ac(i,i+Nyt) =  dx_pt(i)*Fe(i)-Def(i);    
%                 bc(i) = bc(i);
% %           
% %         elseif T_type(i)== 13 %top(S)-right(E)
% %              At(i,:)=0;
% %                  At(i,i) = (dx_et(i)*Fe(i)-dx_wt(i)*Fw(i)+dy_nt(i)*Fn(i)-dy_st(i)*Fs(i)+Def(i)+Dw1(i)+Dnf(i)+Ds1(i)+1/dt);       
% %          At(i,i-Nyt) =-Dw1(i);  
% %          At(i,i+Nyt)= dx_pt(i)*Fe(i)-Def(i);    
% %          At(i,i+1)  =  dy_pt(i)*Fn(i)-Dnf(i);
% %           At(i,i-1) =-Ds1(i);
% %         bt(i)=bt(i);
% %            %Td(i)=bt(i);
% %         elseif T_type(i)== 12 %top(S)-left(E)
% %                 At(i,:)=0;
% %                 At(i,i)      = (dx_et(i)*Fe(i)-dx_wt(i)*Fw(i)+dy_nt(i)*Fn(i)-dy_st(i)*Fs(i)+De1(i)+Dwf(i)+Dnf(i)+Ds1(i)+1/dt);      
% %          At(i,i-Nyt) =-(dx_pt(i)*Fw(i))-Dwf(i);   
% %          At(i,i+Nyt)= -De1(i);  % solid  
% %          At(i,i+1)  = dy_pt(i)*Fn(i)-Dnf(i) ;
% %           At(i,i-1)    = -Ds1(i); %solid
% %              bt(i)=bt(i);     
% % %                 
% %                 
% %                 
% %                 
%          elseif T_type(i)== 14 %interior-solid
%                Ac(i,:)=0;
% %              Ac(i,i)=1;
% %                 At(i,i-Nyt) = 0;   
% %                 At(i,i+Nyt) =  0;
% %                 At(i,i+1)  =  0;
% %                 At(i,i-1)    =  0;
% %               bt(i) =0;
% 
%             Ac(i,i)      = (De1(i)+Dw1(i)+Dn1(i)+Ds1(i)+1/dt);        
%             Ac(i,i-Nyt) = -Dw1(i);   
%             Ac(i,i+Nyt) = -De1(i); 
%              Ac(i,i-1)=-Ds1(i);
%               Ac(i,i+1)= -Dn1(i);
%               bc(i)=bc(i);
%             end
%       end
%       Ac=sparse(Ac);
% %       %bt_prev=T;
%        M=Ac\bc;
% %       %T=T+alpha_T*T; 
%        Mt(n,:)=M';
%        theta_m(n,:)=abs((M-M_left)./(M_left-M_in));
%        M_Center(n,1)=M(tc_type==1);
%        
%        % Non Dimensional moisture content @ center of the domain
% %        Theta_m(n)=abs(M_Center(n)-M_left)./(M_left-M_in);
%        figure(9);
%        
%        fv_plotting(Xt,Yt,reshape(M,[Nyt Nxt]));hold on;
%        quiver(X,Y,U,V,'k')
%        figure(10);
%    [c1,h1]=contour(reshape(M,[Nyt Nxt]),14,'k');
%     clabel(c1,h1);
%    %       quiver(X,Y,U,V,'k')
% % %       bt=T;
% %        colormap jet;
%         drawnow;
% end
% % Theta_m=abs(M_Center-M_left)./(M_left-M_in);
% Theta_m=abs((M_Center-M_left)/(M_left-M_in));
% Theta_t=abs((T_Center-T_left)/(T_left-25));
% Tt=Tt-273;
% for n = 2:nt
%  %Detemining the moisture content at surface level
%        for i = 1:length(T_type(:))
%         if T_type(i)== 9  % Top
%           M_surf(n,i)=(M(i)+M(i-1))/2; 
%            T_surf(n,i)=(Tt(n,i)+Tt(n,i-1))/2;
%         elseif T_type(i)==7
%             M_bot(n,i)=(M(i)+M(i+1))/2;
%             T_bot(n,i)=(Tt(n,i)+Tt(n,i+1))/2;
%         elseif T_type(i)==8
%             M_inlet(n,i)=(M(i)+M(i+Nyt))/2;
%             T_inlet(n,i)=(Tt(n,i)+Tt(n,i+Nyt))/2;       
%        
%         elseif T_type(i)==7
%             M_outlet(n,i)=(M(i)+M(i-Nyt))/2;
%             T_outlet(n,i)=(Tt(n,i)+Tt(n,i-Nyt))/2;
%         end
%        end
%        
% end
% T_surf=T_surf(:,T_type==9);
% for n=2:nt
% T_surf1(n,1)=sum(T_surf(n,:))/length(T_surf(n,:));
% end
% T_center=Tt(:,tc_type==1);
%        
% %    
% % figure(11);
% % for j=1:Nxt*Nyt
% %     plot(Tt(:,j),n,'-k');
% %     drawnow;
% % end
% % figure(12);
% % for j=1:Nxt*Nyt
% %     plot(Mt(:,j),n,'-k');
% %     drawnow;
% % end
% % 
%    
%   t_ha=[0;150;300;450;600;750;900;1050;1200;1350;1500];
%    
%    
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   % xt=xCenter;
% yt=yCenter;
%  k=15;
%  [Xt,Yt]=meshgrid(xt,yt);
%  t_type=p_type;
% dt=0.05;
% t=0:dt:4;
% nt=length(t);
% Nxt=Nx_p;
% Nyt=Ny_p;
%  
%   
%  bt=zeros(Nyt*Nxt,1);
%  
% Fe = ue./dx_p./(dx_p+dx_e); 
% Fw = uw./dx_p./(dx_p+dx_w); 
% Fn = vn./dy_p./(dy_p+dy_n); 
% Fs = vs./dy_p./(dy_p+dy_s);
% 
% %     De = 2/(k*dx^2); Dw = 2/(k*dx^2); Dn = 2/(k*dy.^2); Ds = 2/(k*dy.^2) ;
% 
% De=2*k./dx_p./(dx_p+dx_e);
% Dw=2*k./dx_p./(dx_p+dx_w); 
% Dn=2*k./dy_p./(dy_p+dy_n);
% Ds=2*k./dx_p./(dy_p+dy_s);
% 
%  T_left=4;
%  %T_left=20;
%  Tbot=1;
%  
%  T=ones(Nxt*Nyt,1);
%  %T=exp(-Xt.^2/4);
%  T=reshape(T,[Nyt*Nxt 1]);
% %T=5*T;
% T(p_type==10)=0;
% T(p_type==14)=0;
%   for n=2:nt
%   At=zeros(Nxt*Nyt);
%   bt=T;
%    for i = 1:Nyt*Nxt
%         
%         %% If the cell is ON the left
%         if i <Nyt+1 
%             At(i,i) = 1;
%             At(i,i+Nxt)=1;
%             
%             bt(i)=2*T_left;
%         %% If the cell is ON the right boundary
%         elseif i > Nyt*(Nxt-1)            
%             At(i,i) = 1; 
%             At(i,i-Nyt) = -1;  
%             bt(i) = 0;
%             
%         %% If the cell is IN the domain
%         else
%             At(i,i)      = dx_e(i)*Fe(i)-dx_w(i)*Fw(i)+dy_n(i)*Fn(i)-dy_s(i)*Fs(i)+De(i)+Dw(i)+Dn(i)+Ds(i);        
%             At(i,i-Nyt) = -(dx_p(i)*Fw(i))-Dw(i);   
%             At(i,i+Nyt) =  dx_p(i)*Fe(i)-De(i);    
%             bt(i)=bt(i);
% 
%             %% If the cell BOARDERS the top of the domain
%             if mod(i,Nyt) == 0 
%                 At(i,i)      = At(i,i) - Dn(i);
%                  bt(i) = bt(i) ;                
%             else % If the cell does not boarder the top then 
%                  % there is a North cell
%                 At(i,i+1)    =  dy_p(i)*Fn(i)-Dn(i); 
%             end
% 
%             %% If the cell BOARDERS the bottom of the domain
%             if mod(i,Nyt) == 1
%                 At(i,i)      = At(i,i) + Ds(i);            
%                 bt(i) = bt(i)+2*Tbot*Ds(i);
%             else % If the cell does not boarder the bottom of 
%                  % the domain then there is a South cell                
%                 At(i,i-1)    = -(dy_p(i)*Fs(i))-Ds(i); 
%             end
% 
%         end
%     end
% 
%       for i = 1:length(p_type(:))
%             if p_type(i)== 7
%               At(i,:)=0;
%                At(i,i)      = dx_e(i)*Fe(i)-dx_w(i)*Fw(i)+dy_n(i)*Fn(i)-dy_s(i)*Fs(i)+De(i)+Dw(i)+Dn(i)+Ds(i) +1/dt + Dn(i);  
%                At(i,i-Nyt) = -(dx_p(i)*Fw(i))-Dw(i);   
%             At(i,i+Nyt) =  dx_p(i)*Fe(i)-De(i); 
%             At(i,i-1)=-(dy_p(i)*Fs(i))-Ds(i);
%              At(i,i+1)=0;
%             bt(i) = bt(i)+2*Tbot*Dn(i);
% %                 At(i,i) = 1;      
% %          At(i,i-Nyt) = 0;   
% %          At(i,i+Nyt) =  0;    
% %          At(i,i-1)    =  0; 
% %          At(i,i+1)  =  -1;
% 
%         
%         
%          elseif p_type(i)== 8
%               At(i,:)=0;
%                At(i,i)      = dx_e(i)*Fe(i)-dx_w(i)*Fw(i)+dy_n(i)*Fn(i)-dy_s(i)*Fs(i)+De(i)+Dw(i)+Dn(i)+Ds(i) +1/dt + Dw(i);  
%                At(i,i-Nyt) = -(dx_p(i)*Fw(i))-Dw(i);   
%             %At(i,i+Nyt) =  dx_p(i)*Fe(i)-De(i);    
%             At(i,i+1)= dy_p(i)*Fn(i)-Dn(i);
%             At(i,i-1)=-(dy_p(i)*Fs(i))-Ds(i);
%             At(i,i+Nyt)=0;
%                 bt(i) = bt(i)+2*Tbot*Dw(i);
%             % At(i,i)=1;
% %              At(i,i-Nyt) = 0;   
% %                At(i,i+Nyt) = -1;
% %                At(i,i+1)  =  0;
% %                At(i,i-1)    =  0;
% %               bt(i)=0;
%          elseif p_type(i)==11
%               At(i,:)=0;
%                At(i,i)      = dx_e(i)*Fe(i)-dx_w(i)*Fw(i)+dy_n(i)*Fn(i)-dy_s(i)*Fs(i)+De(i)+Dw(i)+Dn(i)+Ds(i) +1/dt + De(i);
%               % At(i,i-Nyt) = -(dx_p(i)*Fw(i))-Dw(i);   
%             At(i,i+Nyt) =  dx_p(i)*Fe(i)-De(i);    
%             At(i,i+1)= dy_p(i)*Fn(i)-Dn(i);
%             At(i,i-1)=-(dy_p(i)*Fs(i))-Ds(i);
%             At(i,i-Nyt)=0;
%                 bt(i) = bt(i)+2*Tbot*De(i);
% %              At(i,i)=1;
% %              At(i,i-Nyt) = -1;   
% %                At(i,i+Nyt) =  0;
% %                At(i,i+1)  =  0;
% %                At(i,i-1)    =  0;
% %               bt(i)=0;
%          elseif p_type(i)== 9
%               At(i,:)=0;
%                At(i,i)      = dx_e(i)*Fe(i)-dx_w(i)*Fw(i)+dy_n(i)*Fn(i)-dy_s(i)*Fs(i)+De(i)+Dw(i)+Dn(i)+Ds(i) +1/dt + Ds(i);  
%                At(i,i-Nyt) = -(dx_p(i)*Fw(i))-Dw(i); 
%                At(i,i+1)= dy_p(i)*Fn(i)-Dn(i);
%                At(i,i-1)=0;
%             At(i,i+Nyt) =  dx_p(i)*Fe(i)-De(i);    
%                 bt(i) = bt(i)+2*Tbot*Ds(i);
% %               At(i,i)      =  1;      
% %          At(i,i-Nyt) =0;   
% %          At(i,i+Nyt)= 0;    
% %          At(i,i+1)  =  0;
% %           At(i,i-1)    =  -1;
% %           bt(i)=0;
%           
%         elseif p_type(i)== 13 %top(S)-right(E)
%              At(i,:)=0;
%                  At(i,i)      =  2;      
%          At(i,i-Nyt) =0;   
%          At(i,i+Nyt)= -1;    
%          At(i,i+1)  =  0;
%           At(i,i-1)    = -1;
%         bt(i)=0;
%            elseif p_type(i)== 12 %top(S)-left(E)
%                 At(i,:)=0;
%                 At(i,i)      =  2;      
%          At(i,i-Nyt) =-1;   
%          At(i,i+Nyt)= 0;    
%          At(i,i+1)  =  0;
%           At(i,i-1)    =  -1;
%              bt(i)=0;   
%        
%                 
%          elseif p_type(i)== 10 || p_type(i)==14%interior
%               At(i,:)=0;
%             At(i,i)=1;
%                At(i,i-Nyt) = 0;   
%                At(i,i+Nyt) =  0;
%                At(i,i+1)  =  0;
%                At(i,i-1)    =  0;
%              bt(i) =0;
%             end
%       end
%       
%       T=At\bt;
%       
%       figure(6);
%       fv_plotting(Xt,Yt,reshape(T,[Nyt Nxt]));
%       
%       colormap jet;
%       drawnow;
%       
%       
%   end
%  