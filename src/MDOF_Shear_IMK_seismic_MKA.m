function [u, v, a_t, fs_st, fs, T, phi, M, K, C, time, ug_up] = MDOF_Shear_IMK_seismic_MKA(h, wi, Pi, k, Xi, ug, dt, do, Vo, Fy, a_s, dcdy, a_c, g, Name, StrengthLimitCheck,CMBt,ij,zeta)
%% ��ʽHHT����⺯��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%by������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
N = length(wi);      % Number of stories �ṹ�Ĳ���
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
%% Check stability of the method     ��鷽�����ȶ���
Np = length(ug);                    % Length of the record   �������ٶ�ʱ�����̵�����������λΪ[g]
time = 0:dt:(dt*(Np-1));            % ����ʱ��
% If the time step is too large, display a warning. ���ʱ�䲽��̫С�����ʱ�䲽��
% �����¼����ʱ������ ���ٶ������Լ����ֲ���
if dt > T(1)/30
    disp(['Warning: The time step used (dt) for the ground motion "'...
        Name '" is greater than T_1/30. This is not recommended for representing'...
        ' the response correctly.']);
    disp(['A new dt = T/30 = ' num2str(T(1)/30) ' sec is used. In GM: ' Name]);
    
    dt_ = dt/ceil(30*dt/T(1)); %���ֲ��� ceil �����������ȡ��
    time_ = 0:dt_:time(end);   %ʱ������
    ug = interp1(time,ug,time_);
    time = time_;
    dt = dt_;
    Np = length(ug);
end
%%
%% Initial Calculations ��ʼ����
r = ones(N,1);          % Note that this assumes horizontal excitation ע�⣬�����ˮƽ����
P = -M*r*ug*g;          % Equivalent external load  ��Ч���ⲿ����  �õ�������

dy = Fy./k;             % Yielding displacement of each story ÿ�������λ��
dc = dcdy.*dy;          % Capping displacement of each story ��������ʱ��λ��
Fmax = (1-a_s).*Fy+a_s.*k.*dc;      % Positive Strength Limit  ������ǿ�ȼ���
Fmin = -Fmax;                       % Negative Strength Limit  ������ǿ�ȼ���
INt=Fy;
Uy=dy;

LimMax = zeros(N,1);                % ��ʼֵ���ֵ N ����
LimMin = zeros(N,1);                % ��ʼֵ��Сֵ N ����

% Initialize vectors        % ��ʼ������
fs_st = zeros(N,Np);        % Story restoring force        ÿ��ָ��� N ���� Np ���������
fs = fs_st;                 % Total floor restoring force  ���ָ���
fs(:,1) = [-diff(fs_st(:,1));fs_st(end,1)]; % diff ǰ������Ԫ�صĲ�ֵ ��������ָ���
u = zeros(N,Np);            % Relative displacement time history ���λ��ʱ������
v = u;                      % Relative velocity time history ����ٶ�ʱ������
a = u;                      % Relative acceleration time history  ��Լ��ٶ�ʱ������
u(:,1) = do;                % Initial Displacement   ��ʼλ��
v(:,1) = Vo;                % Initial Velocity       ��ʼ�ٶ�
a(:,1) = M\(P(:,1)-C*v(:,1)-fs(:,1));  % Initial Relative Acceleration ��ʼ��Լ��ٶ�
%% dt-1ʱ�̵���Լ��ٶ�
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
%% ѭ������
for i= 1: size(u,2)-1
    v(:,i+1) =  v(:,i) +dt*beta1*a(:,i);                                   %��һʱ���ٶ� % Rel Velocity
    u(:,i+1) = u(:,i) + dt* v(:,i) + dt * dt *beta2 * a(:,i);                          %��һʱ��λ��
    p_ = P(:,i+1)*(1-alpha4)+alpha4* P(:,i);              %��һʱ�̺�����
    if   CMBt == 2
        [fs_st(:,i+1),~,Fmax,Fmin] = SteteDet01(u(:,i),u(:,i+1),fs_st(:,i),k,Fy,a_s,Fmax,Fmin);
    elseif CMBt == 3
        [fs_st(:,i+1),~,Fmax,Fmin,LimMax,LimMin] = StateDet(h,Pi,u(:,max(i-1,1))...
            ,u(:,i),u(:,i+1),fs_st(:,i),k,Fy,a_s,dc,a_c,Fmax,Fmin,LimMax,LimMin...
            ,StrengthLimitCheck);
    else
        [fs_st(:,i+1),~] = StateDet02(u(:,i+1),k);
    end
    fs(:,i+1) = [-diff(fs_st(:,i+1));fs_st(end,i+1)];     % �����ָ���
    fi = fs(:,i+1)*(1-alpha4)+fs(:,i)*alpha4;
    Ci = C*((1-alpha4)*v(:,i+1)+alpha4*v(:,i));
    a(:,i+1) = M1\(p_-fi-Ci-M2*a(:,i));           % Relative acceleration
end
a_t = a/g + r*ug; 	   % Absolute acceleration, in [g] ������ٶ�+���𶯼��ٶ�
ug_up = ug;
end
%% �նȼ��㺯��
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
%% �ָ�������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ������ģ��-����p-deltaЧӦ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fs,kt,Fmax2,Fmin2,LimMax2,LimMin2] = StateDet(h,Pi,u0,u1,u2,fs1,k,Fy,a_s,dc,a_c,Fmax,Fmin,LimMax,LimMin,StrengthLimitCheck)
%  [fs_st(:,i+1),kt,Fmax,Fmin,LimMax,LimMin] = StateDet...
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
%                           0 if strength limit is NOT considered
du0 = [u0(1);diff(u0)];  % u0ʱ�̵Ĳ��λ�� u0(1)��һ�㣨���ײ㣩
du1 = [u1(1);diff(u1)];  % u1ʱ�̵Ĳ��λ�� u1(1)��һ�㣨���ײ㣩
du2 = [u2(1);diff(u2)];  % u2ʱ�̵Ĳ��λ�� u2(1)��һ�㣨���ײ㣩
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
        LimMax2(i) = 1;                      % ��¼
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
%%  ˫����ģ��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fs,kt,Fmax2,Fmin2] = SteteDet01(u1,u2,fs1,k,Fy,a_s,Fmax,Fmin)
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
%% �ߵ���ģ��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
