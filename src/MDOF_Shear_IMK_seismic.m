function [u, v, a_t, fs_st, fs, T, phi, M, K, C, time, ug_up] = MDOF_Shear_IMK_seismic...
    (h, wi, Pi, k, Xi, ug, dt, do, Vo, Fy, a_s, dcdy, a_c, tol, g, Name, StrengthLimitCheck,CMBt,ij,zeta)
%% Newmark-AAM法求解函数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%by郭豪鑫%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [u, v, a_t, fs_st, fs, T, phi, M, K, C, time] = MDOF_ShearBiLinear_seismic...
% (h, wi, Pi, k, Xi, ug, dt, do, Vo, Fy, a_s, dcdy, a_c, tol, g, Name, ...
%                                                       StrengthLimitCheck)
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
%   - tol: Tolerance of the Newton-Raphson iterations.
%           Recommendation: Use tol = 1E-3 to 1E-8.
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
N = length(wi);                               % Number of stories 结构的层数
MaxIter = 20;                                               % 最大的迭代次数
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
% Eigenvalue analysis 特征值分析
% [phi,w2] = eig(K,M);
% w = sqrt(diag(w2));     % Undamped frequencies  Undamped frequencies
% % sort(W)若W是向量不管是列还是行向量，默认都是对A进行升序排列
% [w,index] = sort(w);
% T = 2*pi./w;            % Undamped periods  Undamped frequencies
% % Sort vectors (modal shapes) and normalize them at roof: phi_roof = 1.0；
% % 对矢量（模态形状）进行排序并在屋顶标准化：phi_roof = 1.0 振型标准化 屋顶为1；
% sphi = phi;
% for i = 1:N
%     sphi(:,i) = phi(:,index(i))/ phi(end,index(i));
% end
% phi = sphi;             % Normalized modal shapes
% C matrix  阻尼矩阵
% Mi = diag(phi'*M*phi); %归一化振型对应的广义质量
% if size(Xi) == size(Mi)
%     Ci = 2*Mi.*w.*Xi;
% else
%     Ci = 2*Mi.*w.*Xi';
% end
% C = (phi')^(-1)*diag(Ci)*phi^(-1);
%% Check stability of the method  检查方法的稳定性
Np = length(ug); % Length of the record 包含加速度时间历程的行向量，单位为[g]
time = 0:dt:(dt*(Np-1)); % 计算时间
% If the time step is too large, display a warning.
% 如果时间步长太小则更新时间步长 并更新计算的时间序列 加速度序列以及积分步长
if dt > T(1,1)/30
    disp(['Warning: The time step used (dt) for the ground motion "'...
        Name '" is greater than T_1/30. This is not recommended for representing'...
        ' the response correctly.']);
    disp(['A new dt = T/30 = ' num2str(T(1)/30) ' sec is used. In GM: ' Name]);
    dt_ = dt/ceil(30*dt/T(1));             % 积分步长 ceil 朝正无穷大方向取整
    time_ = 0:dt_:time(end);                                      % 时间向量
    ug = interp1(time,ug,time_);
    time = time_;
    dt = dt_;
    Np = length(ug);
end
%% Initial Calculations 初始计算
r = ones(N,1);% Note that this assumes horizontal excitation 注意，这假设水平激励
P = -M*r*ug*g;         % Equivalent external load  等效的外部负载  得到荷载力
dy = Fy./k;             % Yielding displacement of each story 每层的屈服位移
dc = dcdy.*dy;         % Capping displacement of each story 发生弱化时的位移
Fmax = (1-a_s).*Fy+a_s.*k.*dc;     % Positive Strength Limit  正屈服强度极限
Fmin = -Fmax;                      % Negative Strength Limit  负屈服强度极限
LimMax = zeros(N,1);                                   % 初始值最大值 N 层数
LimMin = zeros(N,1);                                   % 初始值最小值 N 层数
% Initialize vectors                                            % 初始化向量
fs_st = zeros(N,Np);  % Story restoring force 每层恢复力 N 层数 Np 计算的行数
fs = fs_st;                      % Total floor restoring force  总楼层恢复力
fs(:,1) = [-diff(fs_st(:,1));fs_st(end,1)]; % diff前后两个元素的差值用来求层间剪力
Kt = K;                         % Initial tangent stiffness 总的初始切线刚度
kt = k;                                                 % 每层的初始切线刚度
u = zeros(N,Np);       % Relative displacement time history 相对位移时程曲线
v = u;                     % Relative velocity time history 相对速度时程曲线
a = u;              % Relative acceleration time history  相对加速度时程曲线
u(:,1) = do;                               % Initial Displacement   初始位移
v(:,1) = Vo;                               % Initial Velocity       初始速度
a(:,1) = M\(P(:,1)-C*v(:,1)-fs(:,1)); % Initial Relative Acceleration初始相对加速度
% Constants Newmark 计算常数
a1 = 4/dt^2*M + 2/dt*C;
a2 = 4/dt*M + C;
dt2 = dt/10;
a1_2 = 4/dt2^2*M + 2/dt2*C;
a2_2 = 4/dt2*M + C;
R = zeros(N,Np);                   % Unbalanced force history 记录不平衡的力
It = zeros(1,Np);                % Number of iterations     定义迭代次数变量
It(1) = 1;                        % Initialization iterations 初始化迭代变量
%% Calculation for each time step  进行每一时间步的计算
kt_prev = k;  %每层的切线刚度
i = 1;        %循环变量
kk = 0;       %
while i < size(u,2)
    u(:,i+1) = u(:,i);          %把上一时间步的位移赋值给下一时刻，作为下一时刻位移的初始值
    fs_st(:,i+1) = fs_st(:,i);  %把上一时间步的抗力赋值给下一时刻，作为下一时刻抗力的初始值
    fs(:,i+1) = fs(:,i);        % 把上一时间步的层间剪力赋值给下一时刻，作为下一时刻层间剪力的初始值
    p_ = P(:,i+1) + a1*u(:,i) + a2*v(:,i) + M*a(:,i); %不平衡力
    % Newton-Raphson iterations 牛顿-拉普森迭代
    j = 0;                                   % Initialization iterations 迭代次数初始化
    R(:,i+1) = p_ - fs(:,i+1) - a1*u(:,i+1); % 用来检测收敛的不平衡力
    while sum(abs(R(:,i+1)) > tol)   &&  j < MaxIter+1
        Kt_ = Kt + a1;
        du = Kt_\R(:,i+1);
        u(:,i+1) = u(:,i+1) + du;
        % State determination 状态确定
        %% 双线性方法
        if   CMBt == 2
            [fs_st(:,i+1),kt,Fmax,Fmin] = steteDet01(u(:,i),u(:,i+1),fs_st(:,i),k,Fy,a_s,Fmax,Fmin);
            %             [Fy,Uy,fs_st(:,i+1),kt]=steteDet01(Fy,Uy,u(:,i),u(:,i+1),fs_st(:,i),k,INt,a_s);
            %         [fs_st(:,i+1),kt]=steteDet01(u(:,i),u(:,i+1),fs_st(:,i),k,Fy,a_s);
            %% 考虑强度的方法
        elseif  CMBt == 3
            [fs_st(:,i+1),kt,Fmax,Fmin,LimMax,LimMin] = StateDet...
                (h,Pi,u(:,max(i-1,1)),u(:,i),u(:,i+1),fs_st(:,i),k,Fy,a_s,dc,a_c,Fmax,Fmin,LimMax,LimMin,StrengthLimitCheck);
        else
            [fs_st(:,i+1),kt] = StateDet02(u(:,i+1),k);
        end
        % 得到的分别为层间剪力 层间刚度 正的屈服强度 负的屈服强度  屈服状态判断值
        fs(:,i+1) = [-diff(fs_st(:,i+1));fs_st(end,i+1)];    %求解层间恢复力        
        Kt = ComputeK(kt);                                 %计算新的刚度矩阵        
        R(:,i+1) = p_ - fs(:,i+1) - a1*u(:,i+1); % Unbalanced force 不平衡力
        j = j+1;                         % Increase # of iterations 迭代次数       
        if j == MaxIter
            disp(['Warning: Reached ' num2str(MaxIter) ' iterations. Convergence was not'...
                ' achieved in point i = ' num2str(i)])
        end
    end
    if [kt == kt_prev  ;  j < MaxIter+1] % 检查是否发生刚度的变化或者是否收敛
        kk = kk+1;
        % After iterations are completed, compute relative vel and acc:
        % 计算速度和加速度
        It(i+1) = j;                           % Total number of iterations
        v(:,i+1) = 2/dt*(u(:,i+1)-u(:,i)) - v(:,i);          % Rel Velocity
        a(:,i+1) = M\(P(:,i+1)-C*v(:,i+1)-fs(:,i+1)); % Relative acceleration
        kt_prev = kt;
        i = i+1;
    else  % Change in any stiffness or convergence not achieved 未达到收敛或任何刚度的变化
        dug2 = (ug(i+1)-ug(i))/10;                      % 加速度进行插值计算
        
        time_int = (time(i)+dt2):dt2:(time(i+1)-dt2); % 时间处理
        if dug2 == 0                                  % 如果加速度变化为0
            ug_int = ug(i)*ones(1,9);                 % 加速度序列
        else
            ug_int = ug(i)+dug2*(1:9);                % 加速度序列
        end
        
        time = [time(:,1:i) time_int time(:,i+1:end)]; %重新组装时间序列
        ug = [ug(:,1:i) ug_int ug(:,i+1:end)];         %重新组装加速度序列
        
        P = -M*r*ug*g;                                 % 重新组装外荷载序列
        
        u = [u(:,1:i) zeros(N,9) u(:,i+1:end)];        % 重新组装位移序列
        
        fs_st = [fs_st(:,1:i) zeros(N,9) fs_st(:,i+1:end)]; % 重新组装层间剪力
        fs = [fs(:,1:i) zeros(N,9) fs(:,i+1:end)];          % 重新组装层间恢复力
        v = [v(:,1:i) zeros(N,9) v(:,i+1:end)];             % 重新组装速度
        a = [a(:,1:i) zeros(N,9) a(:,i+1:end)];             % 重新组装加速度
        R = [R(:,1:i) zeros(N,9) R(:,i+1:end)];             % 重新组装不平衡力
        It = [It(:,1:i) zeros(1,9) It(:,i+1:end)];          % 重新组装迭代次数
        for i2 = 1:10                                       % 重新计算
            u(:,i+i2) = u(:,i+i2-1);
            fs_st(:,i+i2) = fs_st(:,i+i2-1);
            fs(:,i+i2) = fs(:,i+i2-1);
            p_ = P(:,i+i2) + a1_2*u(:,i+i2-1) + a2_2*v(:,i+i2-1) + M*a(:,i+i2-1);
            % Newton-Raphson iterations
            j = 0;
            R(:,i+i2) = p_ - fs(:,i+i2) - a1_2*u(:,i+i2);
            while sum(abs(R(:,i+i2)) > tol)   &&  j < MaxIter+1
                Kt_ = Kt + a1_2;
                du = Kt_\R(:,i+i2);
                u(:,i+i2) = u(:,i+i2) + du;
                % State determination
                % 双线性方法
                if   CMBt == 2
                    [fs_st(:,i+i2),kt,Fmax,Fmin] = steteDet01(u(:,i+i2-1)...
                        ,u(:,i+i2),fs_st(:,i+i2-1),k,Fy,a_s,Fmax,Fmin);
                elseif CMBt==3
                    % 三线性方法
                    [fs_st(:,i+i2),kt,Fmax,Fmin,LimMax,LimMin] = StateDet...
                        (h,Pi,u(:,i+i2-2),u(:,i+i2-1),u(:,i+i2),fs_st...
                        (:,i+i2-1),k,Fy,a_s,dc,a_c,Fmax,Fmin,LimMax,LimMin,...
                        StrengthLimitCheck);
				else
                    [fs_st(:,i+i2),kt] = StateDet02(u(:,i+i2),k);
                end
                fs(:,i+i2) = [-diff(fs_st(:,i+i2));fs_st(end,i+i2)];
                Kt = ComputeK(kt);
                R(:,i+i2) = p_ - fs(:,i+i2) - a1_2*u(:,i+i2); % Unbalanced force
                j = j+1;                                      % Increase # of iterations
                if j == MaxIter
                    disp(['Warning: Reached ' num2str(MaxIter) ' iterations. Convergence was not'...
                       ' achieved in point i = ' num2str(i)])
                end
            end
            % After iterations are completed, compute relative vel and acc:
            It(i+i2) = j;                                              % Total number of iterations
            v(:,i+i2) = 2/dt2*(u(:,i+i2)-u(:,i+i2-1)) - v(:,i+i2-1);   % Rel Velocity
            a(:,i+i2) = M\(P(:,i+i2)-C*v(:,i+i2)-fs(:,i+i2));          % Relative acceleration
        end
        i = i+10;
        kt_prev = kt;
    end
end
a_t = a/g + r*ug; 	                                                   % Absolute acceleration, in [g] 地面加速度+地震动加速度
ug_up = ug;
end
%% 总刚度矩阵函数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------
%INPUT k层间刚度
%OUPUT K总刚度矩阵
%-------------------------------------------------------------------------
function [K] = ComputeK(k)
if length(k) > 1
    k_aux = k(2:end);
    k_aux(end+1,1) = 0;
    K = diag(k+k_aux) - diag(k(2:end),1) - diag(k(2:end),-1);
else
    K = k;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 三线性方法
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ibarra-Medina-Krawinkler model
function [fs,kt,Fmax2,Fmin2,LimMax2,LimMin2] = StateDet(h,Pi,u0,u1,u2,fs1,...
    k,Fy,a_s,dc,a_c,Fmax,Fmin,LimMax,LimMin,StrengthLimitCheck)
%-----------------------------------------------------------------------------------------------------------------------------------
%     Input
%   [fs_st(:,i+1),kt,Fmax,Fmin,LimMax,LimMin] = StateDet...
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
%                                             0 if strength limit is NOT considered
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% OUPUT
% fs 下一时刻层间剪力 LimMax2 LImMin2 记录状态
% kt 下一时刻层间刚度 Fmax2 Fmin2 屈服力
%--------------------------------------------------------------------------
du0 = [u0(1);diff(u0)]; % u0时刻的层间位移 u0(1)第一层（即底层）
du1 = [u1(1);diff(u1)]; % u1时刻的层间位移 u1(1)第一层（即底层）
du2 = [u2(1);diff(u2)]; % u2时刻的层间位移 u2(1)第一层（即底层）
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
        LimMax2(i) = 1;                      % 记录屈服
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 双线性方法
function [fs,kt,Fmax2,Fmin2] = steteDet01(u1,u2,fs1,k,Fy,a_s,Fmax,Fmin)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%线弹性体系
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






















