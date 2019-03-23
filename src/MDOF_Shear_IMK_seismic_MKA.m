function [u, v, a_t, fs_st, fs, T, phi, M, K, C, time, ug_up] = MDOF_Shear_IMK_seismic_MKA(h, wi, Pi, k, Xi, ug, dt, do, Vo, Fy, a_s, dcdy, a_c, g, Name, StrengthLimitCheck,CMBt,ij,zeta)
%% 显式HHT法求解函数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%by郭豪鑫%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
% Inputs:
%   - h : Vector with story heights, in [length].
%   - wi: Vector with story weights, from 1 to the roof (end of the
%         vector), in [force].
%   - Pi: Vector with P forces for considering P-delta effects, in [force].
%         If Pi = 0, no P-delta effects are considered.
%   - k: Vector with story stiffnesses, from 1 to the roof (end of the
%        vector), in [force/length].
%   - Xi: Vector with damping ratio of each mode.
%   - ug: Row vector containing the acceleration time history, in [g].
%   - dt: Time step. Recommendation: Use dt = min[dt_crit,T/10,sampling of ug]
%           Where   dt_crit = inf ; for Average Acceleration method
%                   dt_crir = 0.551*min(T) ; for Linear Acceleration method
%   - do: Vector with initial displacement, in [length].
%   - Vo: Vector with initial velocity, in [length/sec].
%   - Fy: Vector of Yield Force of each story i, in [force].
%   - a_s: Vector of hardening stiffness ratio of each story.
%   - dcdy: Vector of ductility capacity (dc/dy) for capping point of each story.
%   - a_c: Vector of post-capping stiffness ratio of each story.
%   - g: Standard gravity acceleration, in [length/sec^2].
%   - Name: Name of the ground motion [string].
%   - StrengthLimitCheck:   1 if strength limit is considered (recommended)
%                           0 if strength limit is NOT considered
%   - CMBt   model classes   CMBt==2 % Ibarra-Medina-Krawinkler model
%                            CMBt==3 % bilinear  model
%                            CMBt==4 % linear model
%--------------------------------------------------------------------------
% Outputs:
%   - u: Matrix with Relative Displacement time history, in [length].
%   - v: Matrix with Relative Velocity time history, in [length/sec].
%   - a_t: Matrix with Absolute Acceleration time history, in [g].
%   - fs_st: Matrix with Story Restoring Force (Shear Force) time history, in [force].
%   - fs: Matrix with Story Forces time history, in [force].
%   - T: Vector of periods of the structure, in [sec].
%   - phi: Matrix with modal shapes of the structure.
%   - M, K, C: Mass, Stiffness and Damping matrices assambled, with
%              consistent units.
%   - time: Time vector, in [sec].
%   - ug_up: Updated vector with the acceleration time history, in [g].
%
%   u, v, a_t and Fs, have N rows, with N being equal to the number of
%   stories (N = length(M) = length(K)), and the same number of columns
%   than ug, representing the time variation of the response.
%
%   All the units must be consistent between the input parameters. The
%   output parameters will be also consistent.
%--------------------------------------------------------------------------
N = length(wi);      % Number of stories 结构的层数
%% Obtain T and phi
% Note: first index corresponds to 1st floor, and last index to roof.
% 注意：第一个索引对应于第一层，最后一个索引对应于屋顶。即数据为从下往上
% M and K matrices
M = diag(wi)/g;         % Mass Matrix质量矩阵
K = ComputeK(k);        % Stiffness Matrix 刚度矩阵
if size(M,1)>1
    [C,T,phi] = C_Cal(M,K,ij,zeta);
else
    C = 2*Xi*M*sqrt(K/M);
    T = 2*pi/sqrt(K/M);
    phi = 0;
end
%% Check stability of the method     检查方法的稳定性
Np = length(ug);                    % Length of the record   包含加速度时间历程的行向量，单位为[g]
time = 0:dt:(dt*(Np-1));            % 计算时间
% If the time step is too large, display a warning. 如果时间步长太小则更新时间步长
% 并更新计算的时间序列 加速度序列以及积分步长
if dt > T(1)/30
    disp(['Warning: The time step used (dt) for the ground motion "'...
        Name '" is greater than T_1/30. This is not recommended for representing'...
        ' the response correctly.']);
    disp(['A new dt = T/30 = ' num2str(T(1)/30) ' sec is used. In GM: ' Name]);
    
    dt_ = dt/ceil(30*dt/T(1)); %积分步长 ceil 朝正无穷大方向取整
    time_ = 0:dt_:time(end);   %时间向量
    ug = interp1(time,ug,time_);
    time = time_;
    dt = dt_;
    Np = length(ug);
end
%%
%% Initial Calculations 初始计算
r = ones(N,1);          % Note that this assumes horizontal excitation 注意，这假设水平激励
P = -M*r*ug*g;          % Equivalent external load  等效的外部负载  得到荷载力

dy = Fy./k;             % Yielding displacement of each story 每层的屈服位移
dc = dcdy.*dy;          % Capping displacement of each story 发生弱化时的位移
Fmax = (1-a_s).*Fy+a_s.*k.*dc;      % Positive Strength Limit  正屈服强度极限
Fmin = -Fmax;                       % Negative Strength Limit  负屈服强度极限
INt=Fy;
Uy=dy;

LimMax = zeros(N,1);                % 初始值最大值 N 层数
LimMin = zeros(N,1);                % 初始值最小值 N 层数

% Initialize vectors        % 初始化向量
fs_st = zeros(N,Np);        % Story restoring force        每层恢复力 N 层数 Np 计算的行数
fs = fs_st;                 % Total floor restoring force  层间恢复力
fs(:,1) = [-diff(fs_st(:,1));fs_st(end,1)]; % diff 前后两个元素的差值 用来求层间恢复力
u = zeros(N,Np);            % Relative displacement time history 相对位移时程曲线
v = u;                      % Relative velocity time history 相对速度时程曲线
a = u;                      % Relative acceleration time history  相对加速度时程曲线
u(:,1) = do;                % Initial Displacement   初始位移
v(:,1) = Vo;                % Initial Velocity       初始速度
a(:,1) = M\(P(:,1)-C*v(:,1)-fs(:,1));  % Initial Relative Acceleration 初始相对加速度
%% dt-1时刻的相对加速度
rho=1;
I=eye(N);
alpha3=(2*rho^3+rho^2-1)/(rho^3+rho^2+rho+1);
alpha4=rho/(1+rho);
gamma=1/2+alpha4-alpha3;
beta=(gamma+1/2)/2;
Alpha = M+C*dt*gamma+beta*dt*dt*K;
beta1 = Alpha\M;
beta2 = (1/2+gamma)*beta1;
beta3 = Alpha\(alpha3*M+alpha4*beta*dt*dt*K+alpha4*gamma*dt*C);
M1=M*(I-beta3);
M2=M*beta3;
%% 循环计算
for i= 1: size(u,2)-1
    v(:,i+1) =  v(:,i) +dt*beta1*a(:,i);                                   %下一时刻速度 % Rel Velocity
    u(:,i+1) = u(:,i) + dt* v(:,i) + dt * dt *beta2 * a(:,i);                          %下一时刻位移
    p_ = P(:,i+1)*(1-alpha4)+alpha4* P(:,i);              %下一时刻荷载力
    if   CMBt == 2
        [fs_st(:,i+1),~,Fmax,Fmin] = SteteDet01(u(:,i),u(:,i+1),fs_st(:,i),k,Fy,a_s,Fmax,Fmin);
    elseif CMBt == 3
        [fs_st(:,i+1),~,Fmax,Fmin,LimMax,LimMin] = StateDet(h,Pi,u(:,max(i-1,1))...
            ,u(:,i),u(:,i+1),fs_st(:,i),k,Fy,a_s,dc,a_c,Fmax,Fmin,LimMax,LimMin...
            ,StrengthLimitCheck);
    else
        [fs_st(:,i+1),~] = StateDet02(u(:,i+1),k);
    end
    fs(:,i+1) = [-diff(fs_st(:,i+1));fs_st(end,i+1)];     % 求解层间恢复力
    fi = fs(:,i+1)*(1-alpha4)+fs(:,i)*alpha4;
    Ci = C*((1-alpha4)*v(:,i+1)+alpha4*v(:,i));
    a(:,i+1) = M1\(p_-fi-Ci-M2*a(:,i));           % Relative acceleration
end
a_t = a/g + r*ug; 	   % Absolute acceleration, in [g] 地面加速度+地震动加速度
ug_up = ug;
end
%% 刚度计算函数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 三线性模型-考虑p-delta效应
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fs,kt,Fmax2,Fmin2,LimMax2,LimMin2] = StateDet(h,Pi,u0,u1,u2,fs1,k,Fy,a_s,dc,a_c,Fmax,Fmin,LimMax,LimMin,StrengthLimitCheck)
%  [fs_st(:,i+1),kt,Fmax,Fmin,LimMax,LimMin] = StateDet...
%             (h,Pi,u(:,max(i-1,1)),u(:,i),u(:,i+1),fs_st(:,i),k,Fy,a_s,dc,a_c,Fmax,Fmin,LimMax,LimMin,StrengthLimitCheck);
%   - h/h : Vector with story heights, in [length]. 层高向量
%   - Pi/Pi: Vector with P forces for considering P-delta effects, in [force].If Pi = 0, no P-delta effects are considered.
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
%   - StrengthLimitCheck/StrengthLimitCheck:   1 if strength limit is considered (recommended)
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
    
    if du2(i) > dc(i)                                                              % 弱化阶段 第6阶段
        PosLimEnv = (1-a_s(i))*Fy(i)+a_s(i)*k(i)*dc(i)+a_c(i)*k(i)*(du2(i)-dc(i)); % 第6阶段
        NegLimEnv = (a_s(i)-1)*Fy(i)+a_s(i)*k(i)*du2(i);                           % 第1阶段
        ktEnv = a_c(i)*k(i);                                                       % 层间刚度
    elseif du2(i) > -dc(i)                                                 % 2 3 阶段 负的弱化点--正的弱化点
        PosLimEnv = (1-a_s(i))*Fy(i)+a_s(i)*k(i)*du2(i);                   % 第5段
        NegLimEnv = (a_s(i)-1)*Fy(i)+a_s(i)*k(i)*du2(i);                   % 第2段
        ktEnv = a_s(i)*k(i);                                               % 层间刚度
    else
        PosLimEnv = (1-a_s(i))*Fy(i)+a_s(i)*k(i)*du2(i);                             % 第6段
        NegLimEnv = (a_s(i)-1)*Fy(i)-a_s(i)*k(i)*dc(i)+a_c(i)*k(i)*(du2(i)+dc(i));   % 第1段
        ktEnv = a_c(i)*k(i);                                                         % 层间刚度
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
%%  双线性模型
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fs,kt,Fmax2,Fmin2] = SteteDet01(u1,u2,fs1,k,Fy,a_s,Fmax,Fmin)
%-----------------------------------------------------------------------------------------------------------------------------------
%   Input
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
% OUPUT
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
%% 线弹性模型
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
