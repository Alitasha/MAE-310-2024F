%% 备注
% Developed by GaoYinjun
% 自由度1是x方向，记为u自由度；自由度2是y方向，记为v自由度。
% 分析没有孔的结构，自定了人工解，用理论解以比较误差

clear;
close all;
clc;

%% 参数
traction = 10000.0;                    % Pa or N/m^2
LL         = 4.0;                             % m
E           = 1e9;                            % Pa or N/m^2
Poisson = 0.3;                            % unity
u = E / ( 2 * ( 1 + Poisson ) );     % Pa or N/m^2, shear Modulus
lamda = ( Poisson * E ) / ( (1 + Poisson ) * ( 1 - 2 * Poisson ) );                             % 拉梅系数
D = ( E / ( 1 - Poisson ^ 2))* [ 1, Poisson, 0; Poisson, 1, 0; 0, 0, (1 - Poisson) /2]; % 平面应力
% D = [ lamda + 2 * u, lamda, 0; lamda, lamda+2*u, 0; 0, 0, u];                            % 平面应变
exact = @(x,y) x*(1-x)*y*(1-y);
exact_x = @(x,y) (1-2*x)*y*(1-y);
exact_y = @(x,y) x*(1-x)*(1-2*y);
% 求导
exact_xx = @(x,y) (-2)*y*(1-y);
exact_yy = @(x,y) (-2)*x*(1-x);
exact_xy = @(x,y) (1-2*x)*(1-2*y); %混合导
% 应变 → 应力 → source term， 第一个箭头需要D矩阵，第二个箭头意思是(exact_xx+exact_yy)*D = -f
fx = @(x,y)  ( E / ( 1 - Poisson ^ 2))*(    2.0*x*(1-x) + 2.0*y*(1-y)*Poisson    ); % source term
fy = @(x,y)  ( E / ( 1 - Poisson ^ 2))*(    2.0*x*(1-x)*Poisson  + 2.0*y*(1-y)   ); % source term

% quadrature rule
n_int_xi      = 3;
n_int_eta   = 3;
n_int          = n_int_xi * n_int_eta;
[xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta);

% mesh generation
n_en   = 4;                      % number of nodes in an element
n_el_x = 60;                    % number of elements in x-dir
n_el_y = 60;                    % number of elements in y-dir
n_el   = n_el_x * n_el_y;  % total number of elements
n_sd = 2;                        % number of spatial dimensions

n_np_x = n_el_x + 1;      % number of nodal points in x-dir
n_np_y = n_el_y + 1;      % number of nodal points in y-dir
n_np   = n_np_x * n_np_y; % total number of nodal points
g1 = 0.03;
g2 = -0.03;
x_coor = zeros(n_np, 1);
y_coor = x_coor;

hx = 1.0 / n_el_x;        % mesh size in x-dir
hy = 1.0 / n_el_y;        % mesh size in y-dir

% %% Import mesh 不用import了
for ny = 1 : n_np_y
    for nx = 1 : n_np_x
        index = (ny-1)*n_np_x + nx; % nodal index
        x_coor(index) = (nx-1) * hx;
        y_coor(index) = (ny-1) * hy;
    end
end

% IEN array
IEN = zeros(n_el, n_en);
for ex = 1 : n_el_x
    for ey = 1 : n_el_y
        ee = (ey-1) * n_el_x + ex; % element index
        IEN(ee, 1) = (ey-1) * n_np_x + ex;
        IEN(ee, 2) = (ey-1) * n_np_x + ex + 1;
        IEN(ee, 3) =  ey    * n_np_x + ex + 1;
        IEN(ee, 4) =  ey    * n_np_x + ex;
    end
end

% ID array
ID = zeros(n_np,2);
counter = 0;
for ny = 2 : n_np_y - 1
    for nx = 2 : n_np_x - 1
        index = (ny-1)*n_np_x + nx;
        counter = counter+1;
        ID(index,1) = counter;
        counter = counter +1;% 计数器自增，以为同一全局节点的不同自由度赋唯一的值
        ID(index,2) = counter;
    end
end
ID = ID';

% LM 2自由度
n_eq = counter;
LM = zeros(n_el, 4);
counter2 = 1;
for ee = 1 : n_el
    for i = 1 : 4
        node_index = IEN(ee, i);
        for j = 1 : 2
            LM(counter2, j + (i-1)*2) = ID(j,node_index);
        end
    end
    counter2 = counter2 + 1;
end

% allocate the stiffness matrix and load vector
K = spalloc(n_eq, n_eq, 9 * n_eq); % 9 也是人工设定的
F = zeros(n_eq, 1);
% K = spalloc( n_eq , n_eq ,  n_eq * (2*(n_en-1)+1) ); % 稀疏矩阵内的信息数是多少呢，存疑

for ee = 1 : n_el   % 第一层循环-对每个单元 这个for循环将自由度修正（对比二维传热）
    x_ele = x_coor( IEN(ee, 1:n_en) );
    y_ele = y_coor( IEN(ee, 1:n_en) );
    k_ele = zeros(2*n_en, 2*n_en);  % 原因：2*n_en是因为二维线弹性，位移是两个自由度的
    f_ele = zeros(2*n_en, 1);
    
    for ll = 1 : n_int % 第二层循环-quadrature循环 K_element 和 形函数及其导
        x_l = 0.0; y_l = 0.0;
        dx_dxi = 0.0; dx_deta = 0.0;
        dy_dxi = 0.0; dy_deta = 0.0;
        for aa = 1 : n_en % 第三层循环-对每个单元的局部节点
            x_l = x_l + x_ele(aa) * Quad(aa, xi(ll), eta(ll));
            y_l = y_l + y_ele(aa) * Quad(aa, xi(ll), eta(ll));
            [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
            dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
            dx_deta = dx_deta + x_ele(aa) * Na_eta;
            dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
            dy_deta = dy_deta + y_ele(aa) * Na_eta;
        end
        
        detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;
        
        for aa = 1 : n_en
            Na = Quad(aa, xi(ll), eta(ll));
            [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
            Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
            Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;
            Ba = [ Na_x , 0 ;  0 , Na_y ; Na_y , Na_x ]; % 书本3.10.6 1 ≤ a ≤ nen
            Ba = Ba'; % 积分项计算中位于左边，需转置保证维数正确
            
            % 俩自由度，称为dof
            %             pp = 2*(aa-1);
            %             pp = 2*aa - 1;
            %             pp = 2*(aa-1)+1;
            for dof_x = 1 : 2 % Kpq 先搞定p
                pp1 = 2*(aa-1)+dof_x;
                %                 f_ele(pp) = f_ele(pp) + weight(ll) * Na * detJ * f(x_l,y_l); % source term 应补充
                % 补充完还是不对，检查维数
                % 添加对不同分量的处理：
                
                if dof_x == 1
                    f_ele( pp1 ) = f_ele( pp1 ) + weight( ll ) * detJ * fx ( x_l , y_l ) * Na;
                elseif dof_x == 2
                    f_ele( pp1 ) = f_ele( pp1 ) + weight( ll ) * detJ * fy ( x_l , y_l ) * Na;
                end
                
                for bb = 1 : n_en
                    Nb = Quad(bb, xi(ll), eta(ll));
                    [Nb_xi, Nb_eta] = Quad_grad(bb, xi(ll), eta(ll));
                    Nb_x = (Nb_xi * dy_deta - Nb_eta * dy_dxi) / detJ;
                    Nb_y = (-Nb_xi * dx_deta + Nb_eta * dx_dxi) / detJ;
                    Bb = [ Nb_x , 0 ; 0 , Nb_y ; Nb_y , Nb_x ]; % 书本3.10.6 1 ≤ a ≤ nen 位于右边，不用转置
                    extra_term = weight(ll) * detJ * Ba * D * Bb;
                    for dof_y = 1 : 2
                        % 行和列的影响
                        qq1 = 2*(bb-1)+dof_y;
                        k_ele(pp1, qq1) = k_ele(pp1, qq1) + extra_term(dof_x,dof_y);
                    end
                    
                end % end of bb loop
            end
            
        end % end of aa loop
    end % end of quadrature loop
    
    
    for aa = 1 : n_en
        Note1 = IEN(ee,aa);
        for ii = 1 : 2
            pp1 = 2*(aa-1)+ii;
            if aa == 1
                pp2 = ID(ii,Note1);
            else
                pp2 = ID(ii,Note1);
            end
            if pp2 == 0
                continue;
            end
            F(pp2) = F(pp2) + f_ele(pp1);
            
            for bb = 1 : n_en
                Note2 = IEN(ee,bb);
                for jj = 1 : 2
                    qq1 = 2*(bb-1)+jj;
                    qq2 = ID(jj, Note2);
                    if qq2 > 0
                        K(pp2, qq2) = K(pp2, qq2) + k_ele(pp1, qq1);
                    end
                end
            end
        end
    end
    
end


% solve the stiffness matrix
dn = K \ F;

% insert dn back into the vector for all nodes
disp = zeros( n_np, 2 );
% 先逐列再逐行赋予displacement值，暂未考虑纽曼条件

for ii = 1 : n_np
    index_x = ID(1,ii);
    index_y = ID(2,ii);
    countering = [index_x, index_y];
    if countering(1) > 0
        disp(ii, 1) = dn(countering(1));
    else
        disp(ii,1) = g1;
    end
    if countering(2) > 0
        disp(ii, 2) = dn(countering(2));
    else
        disp(ii,2) = g2;
    end
end


save("Q1", "disp", "n_el_x", "n_el_y");

load("Q1.mat");
hh_x = 1.0 / n_el_x;
hh_y = 1.0 / n_el_y;

% 两个自由度，分x和y呈现
figure(1)

disp_x=disp(:,1);
[X, Y] = meshgrid(0 : hh_x : 1, 0 : hh_y : 1);
Z = reshape(disp_x, n_np_x, n_np_y)';
surf(X, Y, Z);
shading interp
title('X dof: Displacement Field');
xlabel('X');
ylabel('Y');
az = -61;
el = 20;
view(az,el);

% hold on;

% y
figure(2)

disp_y=disp(:,2);
[X, Y] = meshgrid(0 : hh_x : 1, 0 : hh_y : 1);
Z = reshape(disp_y, n_np_x, n_np_y)';
surf(X, Y, Z);
shading interp
title('Y dof: Displacement Field');
xlabel('X');
ylabel('Y');
az = -61;
el = 20;
view(az,el);

% Error and Convergence
% 对每个单元（循环1）找单元坐标，每个积分点（循环2）计算换元，跑所有单元内的节点（循环3）计算dxdxi
lnL2 = zeros(1,n_el);
lnH1 = zeros(1,n_el);
% Error terms

% [xi, eta, weight] = Gauss2D(10, 10);
% L2_top = 0.0; L2_bot = 0.0; H1_top = 0.0; H1_bot = 0.0;
% for ee = 1 : n_el
%     x_ele = x_coor( IEN(ee, 1:n_en) );
%     y_ele = y_coor( IEN(ee, 1:n_en) );
%     u_ele = disp(IEN(ee,1:n_en),:);
%     for ll = 1 : 100
%         x_l = 0.0; y_l = 0.0;
%         dx_dxi = 0.0; dx_deta = 0.0;
%         dy_dxi = 0.0; dy_deta = 0.0;
%         uh = 0.0; uh_xi = 0.0; uh_eta = 0.0;
%         du_dxi = 0.0;
%         for aa = 1 : n_en % 注意u的两个分量
%             x_l = x_l + x_ele(aa) * Quad(aa, xi(ll), eta(ll));
%             y_l = y_l + y_ele(aa) * Quad(aa, xi(ll), eta(ll));
%             %             u_h = uh + u_ele(aa) * Quad(aa, xi(ll), eta(ll));
%             u_h(1) = u_h(1) + u_ele(aa,1) * Quad(aa, xi(ll), eta(ll));
%             u_h(2) = u_h(2) + u_ele(aa,2) * Quad(aa, xi(ll), eta(ll));
%             [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
%             dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
%             dx_deta = dx_deta + x_ele(aa) * Na_eta;
%             dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
%             dy_deta = dy_deta + y_ele(aa) * Na_eta;
%             du_dxi = du_dxi + u_ele(aa,:)  * Na_xi;
%             du_deta = du_deta + u_ele(aa,:)  * Na_eta;
%             
%         end
%         
%         dxi_dx = 1 / dx_dxi;
%         deta_dy = 1 / dy_deta;
%         detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;
%         %         L2_top = L2_top + weight(ll) * (uh - exact(x_l))^2 * dx_dxi;
%         %         L2_bot = L2_bot + weight(ll) * exact(x_l)^2 * dx_dxi;
%         %         from HW 4
%         %         H1_top = H1_top + weight(ll) * ( uh_xi * dxi_dx - exact_x(x_l) )^2 * dx_dxi;
%         %         H1_bot = H1_bot + weight(ll) * exact_x(x_l)^2 * dx_dxi;
%         % 2 dof for exact solution and approximated solution
%         Matrix1 = exact_u(x_l,y_l);
%         Matrix1Der = exact_uu(x_l,y_l);
%         L2_top = L2_top + weight(ll) * (uh - exact_x(x_l))^2 * dx_dxi; % 在Exact后添加对应分量
%         
%         
%         
%     end
%     
%     
%     
% end
% L2_top = sqrt(L2_top); L2_bot = sqrt(L2_bot);
% H1_top = sqrt(H1_top); H1_bot = sqrt(H1_bot);
% 
% L2_error = L2_top / L2_bot;
% H1_error = H1_top / H1_bot;
