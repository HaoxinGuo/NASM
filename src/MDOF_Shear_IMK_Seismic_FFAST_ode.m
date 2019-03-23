function [u, v, a_t, fs_st, fs, T, phi, M, K, C, time, ug_up] =MDOF_Shear_IMK_Seismic_FFAST_ode...
    (h, wi, Pi, k, Xi, ug, dt, do, Vo, Fy, a_s, dcdy, a_c, tol, g, Name, StrengthLimitCheck,CMBt,ij,zeta)
%% 输入：
%        M质量矩阵 wi 层数 K 刚度矩阵 k层间刚度 C阻尼矩阵 g 加速度值
%        ug 地震动加速度序列 tol 收敛系数 Fy屈服力向量 dt  时间步长 a_s 硬化系数 d0/v0 初始位移 速度
N = length(wi);           % Number of stories 结构的层数
M = diag(wi)/g;         % Mass Matrix质量矩阵
K = ComputeK(k);        % Stiffness Matrix 刚度矩阵
if size(M,1)>1
    [C,T,phi] = C_Cal(M,K,ij,zeta);
else
    C = 2*Xi*M*sqrt(K/M);
    T = 2*pi/sqrt(K/M);
    phi = 0;
end
if dt > T(1)/30
    disp(['Warning: The time step used (dt) for the ground motion "'...
        Name '" is greater than T_1/30. This is not recommended for representing'...
        ' the response correctly.']);
end
dy = Fy./k;             % Yielding displacement of each story 每层的屈服位移
dc = dcdy.*dy;          % Capping displacement of each story 发生弱化时的位移
Fmax = (1-a_s).*Fy+a_s.*k.*dc;      % Positive Strength Limit  正屈服强度极限
Fmin = -Fmax;                       % Negative Strength Limit  负屈服强度极限
%% Check stability of the method  检查方法的稳定性
Np = length(ug); % Length of the record 包含加速度时间历程的行向量，单位为[g]
time = 0:dt:(dt*(Np-1)); % 计算时间
%% Initial Calculations 初始计算
r = ones(N,1);% Note that this assumes horizontal excitation 注意，这假设水平激励
P = -M*r*ug*g;         % Equivalent external load  等效的外部负载  得到荷载力
PP=zeros(N,Np-1);    % 恢复力1
QQ=zeros(N,Np-1);   % 恢复力2
%% 求解初始的P Q向量 对力函数插值
for i=1:N
    for j=1:Np-1
        PP(i,j)=1/2*(P(i,j)+P(i,j+1));
        QQ(i,j)=dt/12*(P(i,j)-P(i,j+1));
    end
end
%%  Initialize vectors                                            % 初始化向量
fs_st = zeros(N,Np);  % Story restoring force 每层恢复力 N 层数 Np 计算的行数
fs = fs_st;                      % Total floor restoring force  总楼层恢复力
fs(:,1) = [-diff(fs_st(:,1));fs_st(end,1)]; % diff前后两个元素的差值用来求层间剪力
u = zeros(N,Np);       % Relative displacement time history 相对位移时程曲线
v = u;                     % Relative velocity time history 相对速度时程曲线
u(:,1) = do;                               % Initial Displacement   初始位移
v(:,1) = Vo;                               % Initial Velocity       初始速度
a(:,1) = M\(P(:,1)-C*v(:,1)-fs(:,1)); % Initial Relative Acceleration初始相对加速度
% Constants  计算常数
% It(1) = 1;                        % Initialization iterations 初始化迭代变量
% Kt = K;                         % Initial tangent stiffness 总的初始切线刚度
kt = k;                                                 % 每层的初始切线刚度
U22=-dt/2*(M+dt/6*C);
deltau=zeros(N,Np);
deltav=zeros(N,Np);
deltavfs=zeros(N,Np);% 每次迭代前位移的差值
vfs=zeros(N,Np); % 每次迭代前速度的差值
deltafs2=zeros(N,Np);
deltafs=zeros(N,Np);
deltaR=zeros(N,Np);
% It = zeros(1,Np);                % Number of iterations     定义迭代次数变量
% It(1) = 1;                        % Initialization iterations 初始化迭代变量
%% 循环求解结果
for i=1:Np-1
    %% j<1;
    j=1;
    deltaPP = [dt * (PP(:,i)-fs(:,i)) ; dt*M*v(:,i)-dt^2*QQ(:,i)/2];  % 荷载不平衡力
    A = [C+1/2*dt*K,M-(dt^2)/12*K;M-(dt^2)/12*K,U22];          % 矩阵A
    delta2=A\deltaPP;                                                              % 第一次求得不平衡位移
    delta2u(:,j)=delta2(1:N,1);                                                        % 1:N 为位移
    delta2v(:,j)=delta2(N+1:end,1);                                                 % N+1:end 为速度
    deltau(:,i)=deltau(:,i)+delta2u(:,j);                                                 % 累计位移
    deltav(:,i)=deltav(:,i)+delta2v(:,j);                                                  % 累计速度
    u(:,i+1)=u(:,i)+deltau(:,i);                                                       % 得到位移
    v(:,i+1)=v(:,i)+deltav(:,i);                                                       %得到速度
    if   CMBt == 2
        [fs_st(:,i+1),~,Fmax,Fmin] = SteteDet01(u(:,i),u(:,i+1),fs_st(:,i),...
            k,Fy,a_s,Fmax,Fmin);
    elseif CMBt == 3
        [fs_st(:,i+1),~,Fmax,Fmin,LimMax,LimMin] = StateDet(h,Pi,u(:,max(i-1,1)),...
            u(:,i),u(:,i+1),fs_st(:,i),k,Fy,a_s,dc,a_c,Fmax,Fmin,LimMax,LimMin,...
            StrengthLimitCheck);
    else
        [fs_st(:,i+1),~] = StateDet02(u(:,i+1),k);
    end
    fs(:,i+1) = [-diff(fs_st(:,i+1));fs_st(end,i+1)];                %   求解层间恢复力 fs
    delta2fs(:,j)=fs(:,i+1)-fs(:,i);
    deltafs(:,i+1)=deltafs(:,i)+delta2fs(:,j);                              % delta fs
    deltavfs2(:,j) =  K*[-diff(delta2v);delta2v(end)];
    deltavfs(:,i+1) = deltavfs(:,i)+deltavfs2(:,j);                               % deltafs'
    vfs(:,i+1) = vfs(:,i)+deltafs2(:,i+1);
    %% j>1;
    while(norm(delta2,2)>tol &&  norm(deltaR,2)>tol )
        j=j+1;
        deltaR=[dt/2*(K*delta2u(:,j-1)-delta2fs(:,j-1))+dt^2/12*(deltavfs2(:,j-1)-K*delta2v(:,j-1));
            dt^2/12*(delta2fs(:,j-1)-K*delta2u(:,j-1))];
        delta2=A\deltaR;
        delta2u(:,j)=delta2(1:N,1);                                                        % 1:N 为位移
        delta2v(:,j)=delta2(N+1:end,1);                                                 % N+1:end 为速度
        deltau(:,i)=deltau(:,i)+delta2u(:,j);                                                 % 累计位移
        deltav(:,i)=deltav(:,i)+delta2v(:,j);                                                  % 累计速度
        u(:,i+1)=u(:,i)+deltau(:,i);                                                       % 得到位移
        v(:,i+1)=v(:,i)+deltav(:,i);
        if   CMBt == 2
            [fs_st(:,i+1),~,Fmax,Fmin] = SteteDet01(u(:,i),u(:,i+1),fs_st(:,i),...
                k,Fy,a_s,Fmax,Fmin);
        elseif CMBt == 3
            [fs_st(:,i+1),~,Fmax,Fmin,LimMax,LimMin] = StateDet(h,Pi,u(:,max(i-1,1)),...
                u(:,i),u(:,i+1),fs_st(:,i),k,Fy,a_s,dc,a_c,Fmax,Fmin,LimMax,LimMin,...
                StrengthLimitCheck);
        else
            [fs_st(:,i+1),~] = StateDet02(u(:,i+1),k);
        end
        fs(:,i+1) = [-diff(fs_st(:,i+1));fs_st(end,i+1)];                %   求解层间恢复力 fs
        delta2fs(:,j)=fs(:,i+1)-fs(:,i);
        deltafs(:,i+1)=deltafs(:,i)+delta2fs(:,j);                              % delta fs
        deltavfs2(:,j) =  K*[-diff(delta2v);delta2v(end)];
        deltavfs(:,i+1) = deltavfs(:,i)+deltavfs2(:,j);                               % deltafs'
        vfs(:,i+1) = vfs(:,i)+deltafs2(:,i+1);
    end
    K = ComputeK(kt);                                                 %    计算新的刚度矩阵
    a(:,i+1)=M\(P(:,i+1)-C*v(:,i+1)-fs(:,i+1));
end
a_t = a/g + r*ug; 	   % Absolute acceleration, in [g] 地面加速度+地震动加速度
ug_up=ug;
end
%%
%% 刚度计算函数
function [K] = ComputeK(k)
if length(k) > 1
    k_aux = k(2:end);
    k_aux(end+1,1) = 0;
    K = diag(k+k_aux) - diag(k(2:end),1) - diag(k(2:end),-1);
else
    K = k;
end
end
%% 恢复力曲线
function [fs,kt,Fmax2,Fmin2,LimMax2,LimMin2] = StateDet(h,Pi,u0,u1,u2,fs1,...
    k,Fy,a_s,dc,a_c,Fmax,Fmin,LimMax,LimMin,StrengthLimitCheck)
%  [fs_st(:,i+1),kt,Fmax,Fmin,LimMax,LimMin] = StateDet...
%(h,Pi,u(:,max(i-1,1)),u(:,i),u(:,i+1),fs_st(:,i),k,Fy,a_s,dc,a_c,Fmax,...
%                                   Fmin,LimMax,LimMin,StrengthLimitCheck);
%   - h/h : Vector with story heights, in [length]. 层高向量
%   - Pi/Pi: Vector with P forces for considering P-delta effects, ...
%                   in [force].If Pi = 0, no P-delta effects are considered.
%          P矢量力来考虑P-delta效应，单位 [力]。 如果Pi=0，则不考虑P-delta效应。
%   - u0/u(:,max(i-1,1)),u1/u(:,i),u2/u(:,i+1) i-1时刻 i时刻 i+1时刻的位移向量
%   - fs1/fs_st(:,i): i时刻层间剪力
%   - k/k i时刻每层刚度
%   - Fy/Fy: Vector of Yield Force of each story i, in [force]. 每层的屈服力
%   - a_s/a_s: Vector of hardening stiffness ratio of each story. 每层的硬化比向量
%   - dc/dc:   每层发生强度弱化时的位移
%   - a_c/a_c: Vector of post-capping stiffness ratio of each story. 每层的弱化比向量
%   - Fmax/Fmax  正的屈服强度极限
%   - Fmin/Fmin  负的屈服强度极限
%   - LimMax/LimMax 记录值 记录是否屈服
%   - LimMin/LimMin 记录值 记录是否屈服
%   - StrengthLimitCheck/StrengthLimitCheck:
%                           1 if strength limit is considered (recommended)
%                           0 if strength limit is NOT considered
du0 = [u0(1);diff(u0)];  % u0时刻的层间位移 u0(1)第一层（即底层）
du1 = [u1(1);diff(u1)];  % u1时刻的层间位移 u1(1)第一层（即底层）
du2 = [u2(1);diff(u2)];  % u2时刻的层间位移 u2(1)第一层（即底层）
Fmax2 = Fmax;           % 正的屈服强度极限
Fmin2 = Fmin;           % 负的屈服强度极限
LimMax2 = zeros(length(du1),1);    %
LimMin2 = zeros(length(du1),1);    %
fs = zeros(size(fs1));             % 层间剪力的大小
fs1 = fs1 + Pi./h.*du1;            % P矢量力来考虑P-delta效应
kt = k;                            % 层间刚度
for i = 1:length(du1)              % 单独计算每层层数
    fs(i) = fs1(i) + k(i)*(du2(i)-du1(i));   %计算层间剪力值
    
    if du2(i) > dc(i)                                     % 弱化阶段 第6阶段
        PosLimEnv = (1-a_s(i))*Fy(i)+a_s(i)*k(i)*dc(i)+a_c(i)...
            *k(i)*(du2(i)-dc(i)); % 第6阶段
        NegLimEnv = (a_s(i)-1)*Fy(i)+a_s(i)*k(i)*du2(i);          % 第1阶段
        ktEnv = a_c(i)*k(i);                                      % 层间刚度
    elseif du2(i) > -dc(i)                 % 2 3 阶段 负的弱化点--正的弱化点
        PosLimEnv = (1-a_s(i))*Fy(i)+a_s(i)*k(i)*du2(i);            % 第5段
        NegLimEnv = (a_s(i)-1)*Fy(i)+a_s(i)*k(i)*du2(i);            % 第2段
        ktEnv = a_s(i)*k(i);                                      % 层间刚度
    else
        PosLimEnv = (1-a_s(i))*Fy(i)+a_s(i)*k(i)*du2(i);            % 第6段
        NegLimEnv = (a_s(i)-1)*Fy(i)-a_s(i)*k(i)*dc(i)+a_c(i)*k(i)...
            *(du2(i)+dc(i));  % 第1段
        ktEnv = a_c(i)*k(i);                                      % 层间刚度
    end
    
    if fs(i) > min(Fmax(i),PosLimEnv)        % 利用层间剪力判断 发生屈服
        LimMax2(i) = 1;                      % 记录
        fs(i) = min(Fmax(i),PosLimEnv);      % 更新层间剪力值
        if fs(i) == Fmax(i)                  % 如果层间剪力等于屈服值
            kt(i) = 0;                       % 层间刚度为0
        else
            kt(i) = ktEnv;                   % 否则层间刚度为上一阶段的值
        end
    elseif fs(i) < max(Fmin(i),NegLimEnv)    % 利用层间剪力判断 发生屈服
        LimMin2(i) = 1;                      % 记录状态
        fs(i) = max(Fmin(i),NegLimEnv);      % 更新层间剪力值
        if fs(i) == Fmin(i)                  % 如果层间剪力等于屈服值
            kt(i) = 0;                       % 层间刚度为0
        else
            kt(i) = ktEnv;                   % 否则层间刚度为上一阶段的值
        end
    end
    % -  StrengthLimitCheck / StrengthLimitCheck：如果考虑强度限制，则为1（推荐）
    %                                             如果不考虑强度限制，则为0
    if StrengthLimitCheck
        if du1(i) > dc(i) && du0(i) < du1(i) && du2(i) < du1(i) && LimMax(i)
            Fmax2(i) = fs1(i);
        end
        
        if du1(i) < -dc(i) && du0(i) > du1(i) && du2(i) > du1(i) && LimMin(i)
            Fmin2(i) = fs1(i);
        end
    end
    fs(i) = fs(i) - Pi(i)/h(i)*du2(i);   % 考虑P-delta后的层间剪力
end
end
%% 双线性关系
function [fs,kt,Fmax2,Fmin2] = SteteDet01(u1,u2,fs1,k,Fy,a_s,Fmax,Fmin)
%--------------------------------------------------------------------------
%                               Input
%   [fs_st(:,i+1),kt,Fmax,Fmin,LimMax,LimMin] = StateDet...
%             (h,Pi,u(:,max(i-1,1)),u(:,i),u(:,i+1),fs_st(:,i),k,Fy,a_s
%                   ,dc,a_c,Fmax,Fmin,LimMax,LimMin,StrengthLimitCheck);
%   - u1/u(:,i),u2/u(:,i+1) i-1时刻i时刻 i+1时刻的位移向量
%   - fs1/fs_st(:,i): i时刻层间剪力
%   - k/k i时刻每层刚度
%   - Fy/Fy: Vector of Yield Force of each story i, in [force]. 每层的屈服力
%   - a_s/a_s: Vector of hardening stiffness ratio of each story. 每层的硬化
%                             比向量
%   - Fmax/Fmax  正的屈服强度极限
%   - Fmin/Fmin  负的屈服强度极限
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%                                  OUPUT
% fs 下一时刻层间剪力 LimMax2 LImMin2 记录状态
% kt 下一时刻层间刚度 Fmax2 Fmin2 屈服力
%--------------------------------------------------------------------------
du1 = [u1(1);diff(u1)]; % u1时刻的层间位移 u1(1)第一层（即底层）
du2 = [u2(1);diff(u2)]; % u2时刻的层间位移 u2(1)第一层（即底层）
Fmax2 = Fmax;           % 正的屈服强度极限
Fmin2 = Fmin;           % 负的屈服强度极限
fs = zeros(size(fs1));             % 层间剪力的大小
kt = k;                            % 层间刚度
for i = 1:length(du1)              % 单独计算每层层数
    fs(i) = fs1(i) + k(i)*(du2(i)-du1(i));                  % 计算层间剪力值
    PosLimEnv = (1-a_s(i))*Fy(i)+a_s(i)*k(i)*du2(i);            % 第5段
    NegLimEnv = (a_s(i)-1)*Fy(i)+a_s(i)*k(i)*du2(i);            % 第2段
    ktEnv = a_s(i)*k(i);                                      % 层间刚度
    if fs(i) > min(Fmax(i),PosLimEnv)        % 利用层间剪力判断 发生屈服
        fs(i) = min(Fmax(i),PosLimEnv);      % 更新层间剪力值
        if fs(i) == Fmax(i)                  % 如果层间剪力等于屈服值
            kt(i) = 0;                       % 层间刚度为0
        else
            kt(i) = ktEnv;                   % 否则层间刚度为上一阶段的值
        end
    elseif fs(i) < max(Fmin(i),NegLimEnv)    % 利用层间剪力判断 发生屈服
        fs(i) = max(Fmin(i),NegLimEnv);      % 更新层间剪力值
        if fs(i) == Fmin(i)                  % 如果层间剪力等于屈服值
            kt(i) = 0;                       % 层间刚度为0
        else
            kt(i) = ktEnv;                   % 否则层间刚度为上一阶段的值
        end
    end
end
end
%% 线弹性体系
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fs,kt] = StateDet02(u1,k)
du1 = [u1(1);diff(u1)]; % u1时刻的层间位移 u1(1)第一层（即底层）
kt=k;
fs=zeros(size(k));
for i=1:length(du1)
    %i=1为第一层
    fs(i)=k(i)*(du1(i));
    kt(i)=k(i);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%计算剪切房屋结构参数矩阵函数MCk_Cal
%%计算剪切房屋结构参数矩阵函数MCk_Cal
function[C,T,phi] = C_Cal(M,K,ij,zeta)
N = size(M,1);
i = ij(1,1);
j = ij(1,2);
zeta = zeta(:);
[phi,Ws]=eig(K,M);  %mode shape matrix
sphi = phi;
W=sqrt(diag(Ws)); %circular frequences
[w,index] = sort(W);
for i = 1:N
    sphi(:,i) = phi(:,index(i))/ phi(end,index(i));
end
phi = sphi;             % Normalized modal shapes
T = 2*pi./w;            % Undamped periods  Undamped frequencies
wi = W(i);
wj = W(j);
Wm=[1/wi wi;1/wj wj];
ab=2*(Wm\zeta);
C=ab(1)*M+ab(2)*K; %rayleigh damping matrix
end
