%% 备注
% Developed by GaoYinjun
% 自由度1是x方向，记为u自由度；自由度2是y方向，记为v自由度。
clear;
close all;
clc;

%% 参数
traction = 10000.0;       % Pa or N/m^2
RR        = 0.5;                % m
LL         = 4.0;               % m
radius   = 1;                 % m, 变量
theta     = 0;                 % rad, 变量
E           = 1e9;              % Pa or N/m^2
Poisson = 0.3;               % unity
u = E / ( 2 * ( 1 + Poisson ) );     % Pa or N/m^2, shear Modulus
lamda = ( Poisson * E ) / ( (1 + Poisson ) * ( 1 - 2 * Poisson ) ); % 拉梅系数
% 模量矩阵D:
D = ( E / ( 1 - Poisson ^ 2))* [ 1, Poisson, 0; Poisson, 1, 0; 0, 0, (1 - Poisson) /2]; % 平面应力
% D = [ lamda + 2 * u, lamda, 0; lamda, lamda+2*u, 0; 0, 0, u];  % 平面应变

% %% Exact solution for stress for Plate with Hole
% sigma_rr = @( r , theta )   ( traction / 2 ) * ( 1 - RR ^ 2 / radius ^ 2 ) + ( traction / 2 ) * ( 1 - 4 * RR ^ 2 / radius ^ 2 + 3 * RR ^ 4 / radius ^ 4 ) * cos( 2 * theta );
% sigma_tt = @( r , theta )   ( traction / 2 ) * ( 1 + RR ^ 2 / radius ^ 2 ) - ( traction / 2 ) * ( 1 + 3 * RR ^ 4 / radius ^ 4 ) * cos( 2 * theta );
% sigma_rt = @( r , theta )   - ( traction / 2 ) * ( 1 + 2 * R ^ 2 / radius ^ 2 - 3 * R ^ 4 / radius ^ 4 ) * sin( 2 * theta );

%% Exact manufatured solution for Question 3 Plate without holes
% exact solution
% exact = @(x,y) x*(1-x)*y*(1-y);
% exact_x = @(x,y) (1-2*x)*y*(1-y);
% exact_y = @(x,y) x*(1-x)*(1-2*y);
% 二自由度, u''xx + u''yy = u'' = -f
% exact solution for u
exactu = @(x,y) x*(1-x)*y*(1-y);
exactu_x = @(x,y) (1-2*x)*y*(1-y);
exactu_y = @(x,y) x*(1-x)*(1-2*y);
% source term for u
% exact solution
exactv = @(x,y) x*(1-x)*y*(1-y);
exactv_x = @(x,y) (1-2*x)*y*(1-y);
exactv_y = @(x,y) x*(1-x)*(1-2*y);
% source term for v
f = @(x, y) [ 2 * E / ( Poisson ^ 2 - 1 ); 0 ];

n_en = 4;
n_int = 3;
n_qua = 3;
n_sd = 2;
%% Import mesh
Mesh1; % 运行从Gmsh导出的matlab格式的网格代码，以获取矩阵信息
[ row, col ] = size( msh.QUADS );
IEN= msh.QUADS( : , 1 : col - 1 );  % 维数：单元数 * 4 （2048*4）
nel = row;
x_coor = msh.POS(:,1);
x_coor = msh.POS(:,2);
n_np = length(x_coor); % 验证，15个全局节点


%% (or) Import mesh2
% Mesh1; % 运行从Gmsh导出的matlab格式的网格代码，以获取矩阵信息
% [ row, col ] = size( msh.QUADS );
% IEN= msh.QUADS( : , 1 : col - 1 );  % 维数：单元数 * 4 （2048*4）
% nel = row;


%% (or) Import mesh3
% Mesh1; % 运行从Gmsh导出的matlab格式的网格代码，以获取矩阵信息
% [ row, col ] = size( msh.QUADS );
% IEN= msh.QUADS( : , 1 : col - 1 );  % 维数：单元数 * 4 （2048*4）
% nel = row;

%% (or) Import mesh4
% Mesh1; % 运行从Gmsh导出的matlab格式的网格代码，以获取矩阵信息
% [ row, col ] = size( msh.QUADS );
% IEN= msh.QUADS( : , 1 : col - 1 );  % 维数：单元数 * 4 （2048*4）
% nel = row;

%% ID, IEN and LM Array
n_int_xi = 3;
n_int_eta = 3;
[xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta);

% 边界条件
ID = ones( n_np, 2 );
% 对每一行进行处理
for i = 1: size(msh.LINES,1)
    if msh.LINES(i, 3) == 11 % 底边
        ID(msh.LINES(i, 1),2) = 0; % 限制该边上所有点的ID值y自由度被限制
        ID(msh.LINES(i, 2),2) = 0;
    elseif  msh.LINES(i, 3) == 10 % 左边
        ID(msh.LINES(i, 1),1) = 0; % 限制该边上所有点的ID值x自由度被限制
        ID(msh.LINES(i, 2),1) = 0;
    end
end

% 以上经过验证，Mesh1到Mesh4都没问题

n_np_x = n_el_x + 1;
n_np_y = n_el_y + 1;
% LM = zeros( , 1 );

node_indices = IEN;
node_coords = ID;
LM = node_coords(:, node_indices);
LM = reshape(LM, 2, n_el, n_en);

n_eq = counter * 2; % 二维问题每个节点有两个自由度

% 验证LM组装正确
% ID = [1,3,0,5,7,9,0,12;2,4,0,6,8,10,11,13];
% IEN = [1,1,4,8;2,3,3,6;8,4,5,5;7,2,6,7];
% LM = [1,1,5,12;2,2,6,13;3,0,0,9;4,0,0,10;12,5,7,7;13,6,8,8;0,3,4,0;11,4,10,11];

%% 组装K和F
K = spalloc( n_eq , n_eq ,  n_eq * (2*(n_en-1)+1) ); % 稀疏矩阵内的信息数是多少呢，存疑
F = zeros( n_eq , 1 );
%类似二维传热的求解单元场的循环结构
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
    countering = ID(ii);
    if index > 0
        disp(ii, 1) = dn(countering);
        disp(ii, 2) = dn(countering + 1);
    else 
        disp(ii,1) = g;
        disp(ii,2) = g;
    end
end

save("Q1", "disp", "n_el_x", "n_el_y");

% figure;
% quiver(msh.NODES( : , 1 ), msh.NODES( : , 2 ), disp( : , 2 ), disp( : , 1 ), 'Color', 'r');

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
surf(X, Y, Z);

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
surf(X, Y, Z);

shading interp

az = -61;
el = 20;
view(az,el);