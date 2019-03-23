function [u, v, a_t, fs_st, fs, T, phi, M, K, C, time, ug_up] =MDOF_Shear_IMK_Seismic_FFAST_ode...
    (h, wi, Pi, k, Xi, ug, dt, do, Vo, Fy, a_s, dcdy, a_c, tol, g, Name, StrengthLimitCheck,CMBt,ij,zeta)
%% ���룺
%        M�������� wi ���� K �նȾ��� k���ն� C������� g ���ٶ�ֵ
%        ug ���𶯼��ٶ����� tol ����ϵ�� Fy���������� dt  ʱ�䲽�� a_s Ӳ��ϵ�� d0/v0 ��ʼλ�� �ٶ�
N = length(wi);           % Number of stories �ṹ�Ĳ���
M = diag(wi)/g;         % Mass Matrix��������
K = ComputeK(k);        % Stiffness Matrix �նȾ���
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
dy = Fy./k;             % Yielding displacement of each story ÿ�������λ��
dc = dcdy.*dy;          % Capping displacement of each story ��������ʱ��λ��
Fmax = (1-a_s).*Fy+a_s.*k.*dc;      % Positive Strength Limit  ������ǿ�ȼ���
Fmin = -Fmax;                       % Negative Strength Limit  ������ǿ�ȼ���
%% Check stability of the method  ��鷽�����ȶ���
Np = length(ug); % Length of the record �������ٶ�ʱ�����̵�����������λΪ[g]
time = 0:dt:(dt*(Np-1)); % ����ʱ��
%% Initial Calculations ��ʼ����
r = ones(N,1);% Note that this assumes horizontal excitation ע�⣬�����ˮƽ����
P = -M*r*ug*g;         % Equivalent external load  ��Ч���ⲿ����  �õ�������
PP=zeros(N,Np-1);    % �ָ���1
QQ=zeros(N,Np-1);   % �ָ���2
%% ����ʼ��P Q���� ����������ֵ
for i=1:N
    for j=1:Np-1
        PP(i,j)=1/2*(P(i,j)+P(i,j+1));
        QQ(i,j)=dt/12*(P(i,j)-P(i,j+1));
    end
end
%%  Initialize vectors                                            % ��ʼ������
fs_st = zeros(N,Np);  % Story restoring force ÿ��ָ��� N ���� Np ���������
fs = fs_st;                      % Total floor restoring force  ��¥��ָ���
fs(:,1) = [-diff(fs_st(:,1));fs_st(end,1)]; % diffǰ������Ԫ�صĲ�ֵ�����������
u = zeros(N,Np);       % Relative displacement time history ���λ��ʱ������
v = u;                     % Relative velocity time history ����ٶ�ʱ������
u(:,1) = do;                               % Initial Displacement   ��ʼλ��
v(:,1) = Vo;                               % Initial Velocity       ��ʼ�ٶ�
a(:,1) = M\(P(:,1)-C*v(:,1)-fs(:,1)); % Initial Relative Acceleration��ʼ��Լ��ٶ�
% Constants  ���㳣��
% It(1) = 1;                        % Initialization iterations ��ʼ����������
% Kt = K;                         % Initial tangent stiffness �ܵĳ�ʼ���߸ն�
kt = k;                                                 % ÿ��ĳ�ʼ���߸ն�
U22=-dt/2*(M+dt/6*C);
deltau=zeros(N,Np);
deltav=zeros(N,Np);
deltavfs=zeros(N,Np);% ÿ�ε���ǰλ�ƵĲ�ֵ
vfs=zeros(N,Np); % ÿ�ε���ǰ�ٶȵĲ�ֵ
deltafs2=zeros(N,Np);
deltafs=zeros(N,Np);
deltaR=zeros(N,Np);
% It = zeros(1,Np);                % Number of iterations     ���������������
% It(1) = 1;                        % Initialization iterations ��ʼ����������
%% ѭ�������
for i=1:Np-1
    %% j<1;
    j=1;
    deltaPP = [dt * (PP(:,i)-fs(:,i)) ; dt*M*v(:,i)-dt^2*QQ(:,i)/2];  % ���ز�ƽ����
    A = [C+1/2*dt*K,M-(dt^2)/12*K;M-(dt^2)/12*K,U22];          % ����A
    delta2=A\deltaPP;                                                              % ��һ����ò�ƽ��λ��
    delta2u(:,j)=delta2(1:N,1);                                                        % 1:N Ϊλ��
    delta2v(:,j)=delta2(N+1:end,1);                                                 % N+1:end Ϊ�ٶ�
    deltau(:,i)=deltau(:,i)+delta2u(:,j);                                                 % �ۼ�λ��
    deltav(:,i)=deltav(:,i)+delta2v(:,j);                                                  % �ۼ��ٶ�
    u(:,i+1)=u(:,i)+deltau(:,i);                                                       % �õ�λ��
    v(:,i+1)=v(:,i)+deltav(:,i);                                                       %�õ��ٶ�
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
    fs(:,i+1) = [-diff(fs_st(:,i+1));fs_st(end,i+1)];                %   �����ָ��� fs
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
        delta2u(:,j)=delta2(1:N,1);                                                        % 1:N Ϊλ��
        delta2v(:,j)=delta2(N+1:end,1);                                                 % N+1:end Ϊ�ٶ�
        deltau(:,i)=deltau(:,i)+delta2u(:,j);                                                 % �ۼ�λ��
        deltav(:,i)=deltav(:,i)+delta2v(:,j);                                                  % �ۼ��ٶ�
        u(:,i+1)=u(:,i)+deltau(:,i);                                                       % �õ�λ��
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
        fs(:,i+1) = [-diff(fs_st(:,i+1));fs_st(end,i+1)];                %   �����ָ��� fs
        delta2fs(:,j)=fs(:,i+1)-fs(:,i);
        deltafs(:,i+1)=deltafs(:,i)+delta2fs(:,j);                              % delta fs
        deltavfs2(:,j) =  K*[-diff(delta2v);delta2v(end)];
        deltavfs(:,i+1) = deltavfs(:,i)+deltavfs2(:,j);                               % deltafs'
        vfs(:,i+1) = vfs(:,i)+deltafs2(:,i+1);
    end
    K = ComputeK(kt);                                                 %    �����µĸնȾ���
    a(:,i+1)=M\(P(:,i+1)-C*v(:,i+1)-fs(:,i+1));
end
a_t = a/g + r*ug; 	   % Absolute acceleration, in [g] ������ٶ�+���𶯼��ٶ�
ug_up=ug;
end
%%
%% �նȼ��㺯��
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
function [fs,kt,Fmax2,Fmin2,LimMax2,LimMin2] = StateDet(h,Pi,u0,u1,u2,fs1,...
    k,Fy,a_s,dc,a_c,Fmax,Fmin,LimMax,LimMin,StrengthLimitCheck)
%  [fs_st(:,i+1),kt,Fmax,Fmin,LimMax,LimMin] = StateDet...
%(h,Pi,u(:,max(i-1,1)),u(:,i),u(:,i+1),fs_st(:,i),k,Fy,a_s,dc,a_c,Fmax,...
%                                   Fmin,LimMax,LimMin,StrengthLimitCheck);
%   - h/h : Vector with story heights, in [length]. �������
%   - Pi/Pi: Vector with P forces for considering P-delta effects, ...
%                   in [force].If Pi = 0, no P-delta effects are considered.
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
%   - StrengthLimitCheck/StrengthLimitCheck:
%                           1 if strength limit is considered (recommended)
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
    
    if du2(i) > dc(i)                                     % �����׶� ��6�׶�
        PosLimEnv = (1-a_s(i))*Fy(i)+a_s(i)*k(i)*dc(i)+a_c(i)...
            *k(i)*(du2(i)-dc(i)); % ��6�׶�
        NegLimEnv = (a_s(i)-1)*Fy(i)+a_s(i)*k(i)*du2(i);          % ��1�׶�
        ktEnv = a_c(i)*k(i);                                      % ���ն�
    elseif du2(i) > -dc(i)                 % 2 3 �׶� ����������--����������
        PosLimEnv = (1-a_s(i))*Fy(i)+a_s(i)*k(i)*du2(i);            % ��5��
        NegLimEnv = (a_s(i)-1)*Fy(i)+a_s(i)*k(i)*du2(i);            % ��2��
        ktEnv = a_s(i)*k(i);                                      % ���ն�
    else
        PosLimEnv = (1-a_s(i))*Fy(i)+a_s(i)*k(i)*du2(i);            % ��6��
        NegLimEnv = (a_s(i)-1)*Fy(i)-a_s(i)*k(i)*dc(i)+a_c(i)*k(i)...
            *(du2(i)+dc(i));  % ��1��
        ktEnv = a_c(i)*k(i);                                      % ���ն�
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
%% ˫���Թ�ϵ
function [fs,kt,Fmax2,Fmin2] = SteteDet01(u1,u2,fs1,k,Fy,a_s,Fmax,Fmin)
%--------------------------------------------------------------------------
%                               Input
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
%                                  OUPUT
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
%% �ߵ�����ϵ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
