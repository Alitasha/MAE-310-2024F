%% 备注
% Developed by GaoYinjun
% 自由度1是x方向，记为u自由度；自由度2是y方向，记为v自由度。
% 分析没有孔的结构，自定了人工解，用理论解以比较误差

clear;
close all;
clc;

%% 参数
traction = 10000.0;                    % Pa or N/m^2
% RR        = 0.5;                             % m
LL         = 4.0;                             % m
% radius   = 1;                                % m, 变量
% theta     = 0;                               % rad, 变量
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
n_int_xi  = 3;
n_int_eta = 3;
n_int     = n_int_xi * n_int_eta;
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
g = 0;
x_coor = zeros(n_np, 1);
y_coor = x_coor;

hx = 1.0 / n_el_x;        % mesh size in x-dir
hy = 1.0 / n_el_y;        % mesh size in y-dir

% %% Import mesh 不用import了
% generate the nodal coordinates
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

% LM 2自由度
n_eq = counter;
LM = zeros(n_el, 4);
counter2 = 1;
for ee = 1 : n_el
    for i = 1 : 4
        node_index = IEN(ee, i);
        for j = 1 : 2
            LM(counter2, j + (i-1)*2) = ID(node_index, j);
        end
    end
    counter2 = counter2 + 1;
end

% allocate the stiffness matrix and load vector
K = spalloc(n_eq, n_eq, 9 * n_eq); % 9 也是人工设定的
F = zeros(n_eq, 1);
% K = spalloc( n_eq , n_eq ,  n_eq * (2*(n_en-1)+1) ); % 稀疏矩阵内的信息数是多少呢，存疑
% F = zeros( n_eq , 1 );
% %类似二维传热的求解单元场的循环结构
for ee = 1 : n_el   % 第一层循环-对每个单元 这个for循环将自由度修正（对比二维传热）
    x_ele = x_coor( IEN(ee, 1:n_en) );
    y_ele = y_coor( IEN(ee, 1:n_en) );
    k_ele = zeros(2*n_en, 2*n_en);  % 原因：2*n_en是因为二维线弹性，位移是两个自由度的
    f_ele = zeros(2*n_en, 1);
    
    for ll = 1 : n_int % 第二层循环-quadrature循环 Kelement 和 形函数及其导
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
%             f_ele(aa) = f_ele(aa) + weight(ll) * detJ * f(x_l, y_l) * Na;
%             
            % 俩自由度，称为dof
            
            for dof_x = 1 : 2 % Kpq 先搞定p
                pp = 2*(aa-1)+dof_x;
                f_ele(pp) = f_ele(pp) + weight(ll) * Na * detJ; % source term 应补充
                
                for bb = 1 : n_en 
                    Nb = Quad(bb, xi(ll), eta(ll));
                    [Nb_xi, Nb_eta] = Quad_grad(bb, xi(ll), eta(ll));
                    Nb_x = (Nb_xi * dy_deta - Nb_eta * dy_dxi) / detJ;
                    Nb_y = (-Nb_xi * dx_deta + Nb_eta * dx_dxi) / detJ;
                    Bb = [ Nb_x , 0 ; 0 , Nb_y ; Nb_y , Nb_x ]; % 书本3.10.6 1 ≤ a ≤ nen 位于右边，不用转置

%                     k_ele(aa, bb) = k_ele(aa,bb) + weight(ll) * detJ * kappa * (Na_x * Nb_x + Na_y * Nb_y);
                        extra_term = weight(ll) * detJ * Ba * D * Bb;
                    for dof_y = 1 : 2
                        % 行和列的影响
                        qq = 2*(bb-1)+dof_y;
                        k_ele(pp, qq) = k_ele(pp, qq) + extra_term(dof_x,dof_y);
                    end
                    
                end % end of bb loop
            end
            
        end % end of aa loop
    end % end of quadrature loop
   
 
    
    for aa = 1 : n_en
        % 有两个自由度，则为了组装全局K，每个局部节点（aa）应该分析1方向和2方向。
        % 运用多维IEN、LM表达
        PP1 = LM( ee , 2 * aa-1 ); % 原因：空间维数是二，且二维线弹性，关心两个自由度。
        PP2 = LM( ee , 2 * aa );    % 原因：空间维数是二，且二维线弹性，关心两个自由度。
       
        % 这里先关心 PP1，再关心 PP2 。每一段内有对 QQ1 和 QQ2 的分析
        if PP1 > 0
            F( PP1 ) = F( PP1 ) + f_ele( 2 * aa - 1 );
            for bb = 1 : n_en
                % 和 PP 同样的思路，每个局部节点（aa）循环完后，由于两个自由度，对局部节点（bb）也应该判断1方向和2方向
                QQ1 = LM( ee , 2 * bb - 1 );% 对局部坐标，bb 应变更为 2 * bb - 1
                QQ2 = LM( ee , 2 * bb ); % 对局部坐标，bb 应变更为 2 * bb
                if QQ1 > 0 
                    K( PP1 , QQ1 ) = K( PP1 , QQ1 ) + k_ele( 2 * aa -1 , 2 * bb - 1 );  % 从局部（元素）k 到 全局 K，这里关心 QQ1
                else
                    % 可以在此对纽曼条件g进行补充
                end
                if QQ2 > 0 
                    K( PP1 , QQ2 ) = K( PP1 , QQ1 ) + k_ele( 2 * aa - 1, 2 * bb ); %  从局部（元素）k 到 全局 K，这里关心 QQ2
                else
                    % 可以在此对纽曼条件g进行补充
                end
            end
        end
        
        if PP2 > 0 % 关心的对象从 PP1 变成 PP2
            F( PP2 ) = F( PP2 ) + f_ele( 2 * aa );
            for bb = 1 : n_en
                QQ1 = LM( ee , 2 * bb - 1 );         % 原因：空间维数是二，且二维线弹性，关心两个自由度。 对局部坐标，bb 应变更为 2 * bb - 1
                QQ2 = LM( ee , 2 * bb );               % 原因：空间维数是二，且二维线弹性，关心两个自由度。 对局部坐标，bb 应变更为 2 * bb
                if QQ1 > 0
                    K( PP2 , QQ1 ) = K( PP2 , QQ1 ) + k_ele( 2 * aa , 2 * bb - 1 ); % 从局部（元素） k 到 全局 K，这里关心 QQ1
                else
                    % 可以在此对纽曼条件g进行补充
                end
                 if QQ2 > 0
                    K( PP2 , QQ2 ) = K( PP1 , QQ2 ) + k_ele( 2 * aa , 2 * bb ); % 从局部（元素）k 到 全局 K，这里关心 QQ2
                 else
                    % 可以在此对纽曼条件g进行补充
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
    index_x = ID(ii,1);
    index_y = ID(ii,2);
    countering = [index_x, index_y];
    if countering(1) > 0
        disp(ii, 1) = dn(countering(1));
    else
        disp(ii,1) = g;
    end
    if countering(2) > 0
        disp(ii, 2) = dn(countering(2));
    else
        disp(ii,2) = g;
    end
end


save("Q1", "disp", "n_el_x", "n_el_y");

load("Q1.mat");
hh_x = 1.0 / n_el_x;
hh_y = 1.0 / n_el_y;

% 两个自由度，分x和y呈现
figure(1)
title('X dof: Displacement Field');
xlabel('X');
ylabel('Y');
disp_x=disp(:,1);
[X, Y] = meshgrid(0 : hh_x : 1, 0 : hh_y : 1);
Z = reshape(disp_x, n_np_x, n_np_y)';
% surf(X, Y, Z);

shading interp

az = -61;
el = 20;
view(az,el);

hold on;

% y
figure(2)
title('Y dof: Displacement Field');
xlabel('X');
ylabel('Y');
disp_y=disp(:,2);
[X, Y] = meshgrid(0 : hh_x : 1, 0 : hh_y : 1);
Z = reshape(disp_y, n_np_x, n_np_y)';
% surf(X, Y, Z);

shading interp

az = -61;
el = 20;
view(az,el);










