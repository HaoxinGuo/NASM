function [u, v, a_t, fs_st, fs, T, phi, M, K, C, time, ug_up] = MDOF_Shear_IMK_seismic...
    (h, wi, Pi, k, Xi, ug, dt, do, Vo, Fy, a_s, dcdy, a_c, tol, g, Name, StrengthLimitCheck,CMBt,ij,zeta)
%% Newmark-AAM����⺯��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%by������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
N = length(wi);                               % Number of stories �ṹ�Ĳ���
MaxIter = 20;                                               % ���ĵ�������
%% Obtain T and phi
% Note: first index corresponds to 1st floor, and last index to roof.
% ע�⣺��һ��������Ӧ�ڵ�һ�㣬���һ��������Ӧ���ݶ���������Ϊ��������
% M and K matrices
M = diag(wi)/g;         % Mass Matrix��������
K = ComputeK(k);        % Stiffness Matrix �նȾ���
if size(M,1)>1
    [C,T,phi] = C_Cal(M,K,ij,zeta);
else
    C = 2*Xi*M*sqrt(K/M);
    T = 2*pi/sqrt(K/M);
    phi = 0;
end
% Eigenvalue analysis ����ֵ����
% [phi,w2] = eig(K,M);
% w = sqrt(diag(w2));     % Undamped frequencies  Undamped frequencies
% % sort(W)��W�������������л�����������Ĭ�϶��Ƕ�A������������
% [w,index] = sort(w);
% T = 2*pi./w;            % Undamped periods  Undamped frequencies
% % Sort vectors (modal shapes) and normalize them at roof: phi_roof = 1.0��
% % ��ʸ����ģ̬��״�������������ݶ���׼����phi_roof = 1.0 ���ͱ�׼�� �ݶ�Ϊ1��
% sphi = phi;
% for i = 1:N
%     sphi(:,i) = phi(:,index(i))/ phi(end,index(i));
% end
% phi = sphi;             % Normalized modal shapes
% C matrix  �������
% Mi = diag(phi'*M*phi); %��һ�����Ͷ�Ӧ�Ĺ�������
% if size(Xi) == size(Mi)
%     Ci = 2*Mi.*w.*Xi;
% else
%     Ci = 2*Mi.*w.*Xi';
% end
% C = (phi')^(-1)*diag(Ci)*phi^(-1);
%% Check stability of the method  ��鷽�����ȶ���
Np = length(ug); % Length of the record �������ٶ�ʱ�����̵�����������λΪ[g]
time = 0:dt:(dt*(Np-1)); % ����ʱ��
% If the time step is too large, display a warning.
% ���ʱ�䲽��̫С�����ʱ�䲽�� �����¼����ʱ������ ���ٶ������Լ����ֲ���
if dt > T(1,1)/30
    disp(['Warning: The time step used (dt) for the ground motion "'...
        Name '" is greater than T_1/30. This is not recommended for representing'...
        ' the response correctly.']);
    disp(['A new dt = T/30 = ' num2str(T(1)/30) ' sec is used. In GM: ' Name]);
    dt_ = dt/ceil(30*dt/T(1));             % ���ֲ��� ceil �����������ȡ��
    time_ = 0:dt_:time(end);                                      % ʱ������
    ug = interp1(time,ug,time_);
    time = time_;
    dt = dt_;
    Np = length(ug);
end
%% Initial Calculations ��ʼ����
r = ones(N,1);% Note that this assumes horizontal excitation ע�⣬�����ˮƽ����
P = -M*r*ug*g;         % Equivalent external load  ��Ч���ⲿ����  �õ�������
dy = Fy./k;             % Yielding displacement of each story ÿ�������λ��
dc = dcdy.*dy;         % Capping displacement of each story ��������ʱ��λ��
Fmax = (1-a_s).*Fy+a_s.*k.*dc;     % Positive Strength Limit  ������ǿ�ȼ���
Fmin = -Fmax;                      % Negative Strength Limit  ������ǿ�ȼ���
LimMax = zeros(N,1);                                   % ��ʼֵ���ֵ N ����
LimMin = zeros(N,1);                                   % ��ʼֵ��Сֵ N ����
% Initialize vectors                                            % ��ʼ������
fs_st = zeros(N,Np);  % Story restoring force ÿ��ָ��� N ���� Np ���������
fs = fs_st;                      % Total floor restoring force  ��¥��ָ���
fs(:,1) = [-diff(fs_st(:,1));fs_st(end,1)]; % diffǰ������Ԫ�صĲ�ֵ�����������
Kt = K;                         % Initial tangent stiffness �ܵĳ�ʼ���߸ն�
kt = k;                                                 % ÿ��ĳ�ʼ���߸ն�
u = zeros(N,Np);       % Relative displacement time history ���λ��ʱ������
v = u;                     % Relative velocity time history ����ٶ�ʱ������
a = u;              % Relative acceleration time history  ��Լ��ٶ�ʱ������
u(:,1) = do;                               % Initial Displacement   ��ʼλ��
v(:,1) = Vo;                               % Initial Velocity       ��ʼ�ٶ�
a(:,1) = M\(P(:,1)-C*v(:,1)-fs(:,1)); % Initial Relative Acceleration��ʼ��Լ��ٶ�
% Constants Newmark ���㳣��
a1 = 4/dt^2*M + 2/dt*C;
a2 = 4/dt*M + C;
dt2 = dt/10;
a1_2 = 4/dt2^2*M + 2/dt2*C;
a2_2 = 4/dt2*M + C;
R = zeros(N,Np);                   % Unbalanced force history ��¼��ƽ�����
It = zeros(1,Np);                % Number of iterations     ���������������
It(1) = 1;                        % Initialization iterations ��ʼ����������
%% Calculation for each time step  ����ÿһʱ�䲽�ļ���
kt_prev = k;  %ÿ������߸ն�
i = 1;        %ѭ������
kk = 0;       %
while i < size(u,2)
    u(:,i+1) = u(:,i);          %����һʱ�䲽��λ�Ƹ�ֵ����һʱ�̣���Ϊ��һʱ��λ�Ƶĳ�ʼֵ
    fs_st(:,i+1) = fs_st(:,i);  %����һʱ�䲽�Ŀ�����ֵ����һʱ�̣���Ϊ��һʱ�̿����ĳ�ʼֵ
    fs(:,i+1) = fs(:,i);        % ����һʱ�䲽�Ĳ�������ֵ����һʱ�̣���Ϊ��һʱ�̲������ĳ�ʼֵ
    p_ = P(:,i+1) + a1*u(:,i) + a2*v(:,i) + M*a(:,i); %��ƽ����
    % Newton-Raphson iterations ţ��-����ɭ����
    j = 0;                                   % Initialization iterations ����������ʼ��
    R(:,i+1) = p_ - fs(:,i+1) - a1*u(:,i+1); % ������������Ĳ�ƽ����
    while sum(abs(R(:,i+1)) > tol)   &&  j < MaxIter+1
        Kt_ = Kt + a1;
        du = Kt_\R(:,i+1);
        u(:,i+1) = u(:,i+1) + du;
        % State determination ״̬ȷ��
        %% ˫���Է���
        if   CMBt == 2
            [fs_st(:,i+1),kt,Fmax,Fmin] = steteDet01(u(:,i),u(:,i+1),fs_st(:,i),k,Fy,a_s,Fmax,Fmin);
            %             [Fy,Uy,fs_st(:,i+1),kt]=steteDet01(Fy,Uy,u(:,i),u(:,i+1),fs_st(:,i),k,INt,a_s);
            %         [fs_st(:,i+1),kt]=steteDet01(u(:,i),u(:,i+1),fs_st(:,i),k,Fy,a_s);
            %% ����ǿ�ȵķ���
        elseif  CMBt == 3
            [fs_st(:,i+1),kt,Fmax,Fmin,LimMax,LimMin] = StateDet...
                (h,Pi,u(:,max(i-1,1)),u(:,i),u(:,i+1),fs_st(:,i),k,Fy,a_s,dc,a_c,Fmax,Fmin,LimMax,LimMin,StrengthLimitCheck);
        else
            [fs_st(:,i+1),kt] = StateDet02(u(:,i+1),k);
        end
        % �õ��ķֱ�Ϊ������ ���ն� ��������ǿ�� ��������ǿ��  ����״̬�ж�ֵ
        fs(:,i+1) = [-diff(fs_st(:,i+1));fs_st(end,i+1)];    %�����ָ���        
        Kt = ComputeK(kt);                                 %�����µĸնȾ���        
        R(:,i+1) = p_ - fs(:,i+1) - a1*u(:,i+1); % Unbalanced force ��ƽ����
        j = j+1;                         % Increase # of iterations ��������       
        if j == MaxIter
            disp(['Warning: Reached ' num2str(MaxIter) ' iterations. Convergence was not'...
                ' achieved in point i = ' num2str(i)])
        end
    end
    if [kt == kt_prev  ;  j < MaxIter+1] % ����Ƿ����նȵı仯�����Ƿ�����
        kk = kk+1;
        % After iterations are completed, compute relative vel and acc:
        % �����ٶȺͼ��ٶ�
        It(i+1) = j;                           % Total number of iterations
        v(:,i+1) = 2/dt*(u(:,i+1)-u(:,i)) - v(:,i);          % Rel Velocity
        a(:,i+1) = M\(P(:,i+1)-C*v(:,i+1)-fs(:,i+1)); % Relative acceleration
        kt_prev = kt;
        i = i+1;
    else  % Change in any stiffness or convergence not achieved δ�ﵽ�������κθնȵı仯
        dug2 = (ug(i+1)-ug(i))/10;                      % ���ٶȽ��в�ֵ����
        
        time_int = (time(i)+dt2):dt2:(time(i+1)-dt2); % ʱ�䴦��
        if dug2 == 0                                  % ������ٶȱ仯Ϊ0
            ug_int = ug(i)*ones(1,9);                 % ���ٶ�����
        else
            ug_int = ug(i)+dug2*(1:9);                % ���ٶ�����
        end
        
        time = [time(:,1:i) time_int time(:,i+1:end)]; %������װʱ������
        ug = [ug(:,1:i) ug_int ug(:,i+1:end)];         %������װ���ٶ�����
        
        P = -M*r*ug*g;                                 % ������װ���������
        
        u = [u(:,1:i) zeros(N,9) u(:,i+1:end)];        % ������װλ������
        
        fs_st = [fs_st(:,1:i) zeros(N,9) fs_st(:,i+1:end)]; % ������װ������
        fs = [fs(:,1:i) zeros(N,9) fs(:,i+1:end)];          % ������װ���ָ���
        v = [v(:,1:i) zeros(N,9) v(:,i+1:end)];             % ������װ�ٶ�
        a = [a(:,1:i) zeros(N,9) a(:,i+1:end)];             % ������װ���ٶ�
        R = [R(:,1:i) zeros(N,9) R(:,i+1:end)];             % ������װ��ƽ����
        It = [It(:,1:i) zeros(1,9) It(:,i+1:end)];          % ������װ��������
        for i2 = 1:10                                       % ���¼���
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
                % ˫���Է���
                if   CMBt == 2
                    [fs_st(:,i+i2),kt,Fmax,Fmin] = steteDet01(u(:,i+i2-1)...
                        ,u(:,i+i2),fs_st(:,i+i2-1),k,Fy,a_s,Fmax,Fmin);
                elseif CMBt==3
                    % �����Է���
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
a_t = a/g + r*ug; 	                                                   % Absolute acceleration, in [g] ������ٶ�+���𶯼��ٶ�
ug_up = ug;
end
%% �ܸնȾ�����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------
%INPUT k���ն�
%OUPUT K�ܸնȾ���
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
%% �����Է���
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ibarra-Medina-Krawinkler model
function [fs,kt,Fmax2,Fmin2,LimMax2,LimMin2] = StateDet(h,Pi,u0,u1,u2,fs1,...
    k,Fy,a_s,dc,a_c,Fmax,Fmin,LimMax,LimMin,StrengthLimitCheck)
%-----------------------------------------------------------------------------------------------------------------------------------
%     Input
%   [fs_st(:,i+1),kt,Fmax,Fmin,LimMax,LimMin] = StateDet...
%             (h,Pi,u(:,max(i-1,1)),u(:,i),u(:,i+1),fs_st(:,i),k,Fy,a_s,dc,a_c,Fmax,Fmin,LimMax,LimMin,StrengthLimitCheck);
%   - h/h : Vector with story heights, in [length]. �������
%   - Pi/Pi: Vector with P forces for considering P-delta effects, in [force].If Pi = 0, no P-delta effects are considered.
%          Pʸ����������P-deltaЧӦ����λ [��]�� ���Pi=0���򲻿���P-deltaЧӦ��
%   - u0/u(:,max(i-1,1)),u1/u(:,i),u2/u(:,i+1) i-1ʱ�� iʱ�� i+1ʱ�̵�λ������
%   - fs1/fs_st(:,i): iʱ�̲�����
%   - k/k iʱ��ÿ��ն�
%   - Fy/Fy: Vector of Yield Force of each story i, in [force]. ÿ���������
%   - a_s/a_s: Vector of hardening stiffness ratio of each story. ÿ���Ӳ��������
%   - dc/dc:   ÿ�㷢��ǿ������ʱ��λ��
%   - a_c/a_c: Vector of post-capping stiffness ratio of each story. ÿ�������������
%   - Fmax/Fmax  ��������ǿ�ȼ���
%   - Fmin/Fmin  ��������ǿ�ȼ���
%   - LimMax/LimMax ��¼ֵ ��¼�Ƿ�����
%   - LimMin/LimMin ��¼ֵ ��¼�Ƿ�����
%   - StrengthLimitCheck/StrengthLimitCheck:   1 if strength limit is considered (recommended)
%                                             0 if strength limit is NOT considered
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% OUPUT
% fs ��һʱ�̲����� LimMax2 LImMin2 ��¼״̬
% kt ��һʱ�̲��ն� Fmax2 Fmin2 ������
%--------------------------------------------------------------------------
du0 = [u0(1);diff(u0)]; % u0ʱ�̵Ĳ��λ�� u0(1)��һ�㣨���ײ㣩
du1 = [u1(1);diff(u1)]; % u1ʱ�̵Ĳ��λ�� u1(1)��һ�㣨���ײ㣩
du2 = [u2(1);diff(u2)]; % u2ʱ�̵Ĳ��λ�� u2(1)��һ�㣨���ײ㣩
Fmax2 = Fmax;           % ��������ǿ�ȼ���
Fmin2 = Fmin;           % ��������ǿ�ȼ���
LimMax2 = zeros(length(du1),1);    %
LimMin2 = zeros(length(du1),1);    %

fs = zeros(size(fs1));             % �������Ĵ�С
fs1 = fs1 + Pi./h.*du1;            % Pʸ����������P-deltaЧӦ
kt = k;                            % ���ն�
for i = 1:length(du1)              % ��������ÿ�����
    fs(i) = fs1(i) + k(i)*(du2(i)-du1(i));   %���������ֵ
    
    if du2(i) > dc(i)                                                              % �����׶� ��6�׶�
        PosLimEnv = (1-a_s(i))*Fy(i)+a_s(i)*k(i)*dc(i)+a_c(i)*k(i)*(du2(i)-dc(i)); % ��6�׶�
        NegLimEnv = (a_s(i)-1)*Fy(i)+a_s(i)*k(i)*du2(i);                           % ��1�׶�
        ktEnv = a_c(i)*k(i);                                                       % ���ն�
    elseif du2(i) > -dc(i)                                                 % 2 3 �׶� ����������--����������
        PosLimEnv = (1-a_s(i))*Fy(i)+a_s(i)*k(i)*du2(i);                   % ��5��
        NegLimEnv = (a_s(i)-1)*Fy(i)+a_s(i)*k(i)*du2(i);                   % ��2��
        ktEnv = a_s(i)*k(i);                                               % ���ն�
    else
        PosLimEnv = (1-a_s(i))*Fy(i)+a_s(i)*k(i)*du2(i);                             % ��6��
        NegLimEnv = (a_s(i)-1)*Fy(i)-a_s(i)*k(i)*dc(i)+a_c(i)*k(i)*(du2(i)+dc(i));   % ��1��
        ktEnv = a_c(i)*k(i);                                                         % ���ն�
    end
    
    if fs(i) > min(Fmax(i),PosLimEnv)        % ���ò������ж� ��������
        LimMax2(i) = 1;                      % ��¼����
        fs(i) = min(Fmax(i),PosLimEnv);      % ���²�����ֵ
        if fs(i) == Fmax(i)                  % �����������������ֵ
            kt(i) = 0;                       % ���ն�Ϊ0
        else
            kt(i) = ktEnv;                   % ������ն�Ϊ��һ�׶ε�ֵ
        end
    elseif fs(i) < max(Fmin(i),NegLimEnv)    % ���ò������ж� ��������
        LimMin2(i) = 1;                      % ��¼״̬
        fs(i) = max(Fmin(i),NegLimEnv);      % ���²�����ֵ
        if fs(i) == Fmin(i)                  % �����������������ֵ
            kt(i) = 0;                       % ���ն�Ϊ0
        else
            kt(i) = ktEnv;                   % ������ն�Ϊ��һ�׶ε�ֵ
        end
    end
    % -  StrengthLimitCheck / StrengthLimitCheck���������ǿ�����ƣ���Ϊ1���Ƽ���
    %                                             ���������ǿ�����ƣ���Ϊ0
    if StrengthLimitCheck
        if du1(i) > dc(i) && du0(i) < du1(i) && du2(i) < du1(i) && LimMax(i)
            Fmax2(i) = fs1(i);
        end
        
        if du1(i) < -dc(i) && du0(i) > du1(i) && du2(i) > du1(i) && LimMin(i)
            Fmin2(i) = fs1(i);
        end
    end
    fs(i) = fs(i) - Pi(i)/h(i)*du2(i);   % ����P-delta��Ĳ�����
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ˫���Է���
function [fs,kt,Fmax2,Fmin2] = steteDet01(u1,u2,fs1,k,Fy,a_s,Fmax,Fmin)
%-----------------------------------------------------------------------------------------------------------------------------------
%   Input
%   [fs_st(:,i+1),kt,Fmax,Fmin,LimMax,LimMin] = StateDet...
%             (h,Pi,u(:,max(i-1,1)),u(:,i),u(:,i+1),fs_st(:,i),k,Fy,a_s
%                   ,dc,a_c,Fmax,Fmin,LimMax,LimMin,StrengthLimitCheck);
%   - u1/u(:,i),u2/u(:,i+1) i-1ʱ��iʱ�� i+1ʱ�̵�λ������
%   - fs1/fs_st(:,i): iʱ�̲�����
%   - k/k iʱ��ÿ��ն�
%   - Fy/Fy: Vector of Yield Force of each story i, in [force]. ÿ���������
%   - a_s/a_s: Vector of hardening stiffness ratio of each story. ÿ���Ӳ��
%                             ������
%   - Fmax/Fmax  ��������ǿ�ȼ���
%   - Fmin/Fmin  ��������ǿ�ȼ���
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% OUPUT
% fs ��һʱ�̲����� LimMax2 LImMin2 ��¼״̬
% kt ��һʱ�̲��ն� Fmax2 Fmin2 ������
%--------------------------------------------------------------------------
du1 = [u1(1);diff(u1)]; % u1ʱ�̵Ĳ��λ�� u1(1)��һ�㣨���ײ㣩
du2 = [u2(1);diff(u2)]; % u2ʱ�̵Ĳ��λ�� u2(1)��һ�㣨���ײ㣩
Fmax2 = Fmax;           % ��������ǿ�ȼ���
Fmin2 = Fmin;           % ��������ǿ�ȼ���
fs = zeros(size(fs1));             % �������Ĵ�С
kt = k;                            % ���ն�
for i = 1:length(du1)              % ��������ÿ�����
    fs(i) = fs1(i) + k(i)*(du2(i)-du1(i));                  % ���������ֵ
    PosLimEnv = (1-a_s(i))*Fy(i)+a_s(i)*k(i)*du2(i);            % ��5��
    NegLimEnv = (a_s(i)-1)*Fy(i)+a_s(i)*k(i)*du2(i);            % ��2��
    ktEnv = a_s(i)*k(i);                                      % ���ն�
    if fs(i) > min(Fmax(i),PosLimEnv)        % ���ò������ж� ��������
        fs(i) = min(Fmax(i),PosLimEnv);      % ���²�����ֵ
        if fs(i) == Fmax(i)                  % �����������������ֵ
            kt(i) = 0;                       % ���ն�Ϊ0
        else
            kt(i) = ktEnv;                   % ������ն�Ϊ��һ�׶ε�ֵ
        end
    elseif fs(i) < max(Fmin(i),NegLimEnv)    % ���ò������ж� ��������
        fs(i) = max(Fmin(i),NegLimEnv);      % ���²�����ֵ
        if fs(i) == Fmin(i)                  % �����������������ֵ
            kt(i) = 0;                       % ���ն�Ϊ0
        else
            kt(i) = ktEnv;                   % ������ն�Ϊ��һ�׶ε�ֵ
        end
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%�ߵ�����ϵ
function [fs,kt] = StateDet02(u1,k)
du1 = [u1(1);diff(u1)]; % u1ʱ�̵Ĳ��λ�� u1(1)��һ�㣨���ײ㣩
kt=k;
fs=zeros(size(k));
for i=1:length(du1)
    %i=1Ϊ��һ��
    fs(i)=k(i)*(du1(i));
    kt(i)=k(i);
end
end
%%������з��ݽṹ����������MCk_Cal
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






















