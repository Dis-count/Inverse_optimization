% ����������Լ��
% ��һ��Լ��ȷ������ r_ik (Ϊ1�����ֵ<Ϊ0 ����Сֵ )�Ĵ�С��ϵ

V_UFL = 28;   % ����������ֵ��Ŀ��ֵ��

vi = [1 ; 1 ; 0;];

uik = [ 1; 0 ; 1 ;
        0; 1 ; 0 ;
        0; 0 ; 0 ;];

m = length(vi);
n = length(uik)/m;

x_0 = [vi;uik];
x0 =[x_0;x_0*(-1)];  %�����Ҫ��double����

% ������ԭ����Costs

FC    = [5; 6; 7;];

TC    = [11; 4; 8;
         5; 7; 10;
         19; 6; 3;];

Costs =[FC;TC];

V_0 = x_0'*Costs;    %  ���������Ž����� ��Ӧ��ԭ�ɱ�����Ļ���ֵ

v1 = find(vi == 1);   % v1��¼ vi=1 �ĺ����� �� ������   Ҫ�ֱ��ҵ�viΪ 0 �� 1 ���±�

v0 = find(vi == 0);

u1 = cell(v1,2);

ui = reshape(uik,m,n)';   % ת�� �õ�uik �� �����

TCi = reshape(TC,m,n)';   % �õ�����

s = 0; %  ���ڼ�¼ vi �� rik==1 ����Ϊ1 ������

for i = v1'  % v1 ��Ҫ��������
    t=1;

    [umax,in] = max(TCi(t,:));

    u1{t,1} = i;  % ��Ԫ���� ��һ�������� viΪ 1�� ���� ��Ӧ�е� ��Сֵ�����ֵ
    u1{t,2} = find(ui(t,:) == 1);  %��Ԫ���� �ڶ����������� rik ��Ӧ���� Ϊ1 ������

    if length(find(ui(t,:) == 1))==1
        s=s+1;
    end

    t=t+1;
end


u0 = cell(m-length(v1),2);  % �˴� ���¶��� һ�����Ƶĵ�Ԫ���� ���Ϊ0 �Ĳ���

for i = v0'  % v1 ��Ҫ��������
    t=1;

    u1{t,1} = i;  % ��Ԫ���� ��һ�������� viΪ 1������
    u1{t,2} = find(ui(t,:) == 0);  %��Ԫ���� �ڶ����������� rik ��Ӧ���� Ϊ0 ������

    t=t+1;
end

fi;

% ��Ҫ���� vi Ϊ1���±�  �Լ� Ϊ0 ���±�

rik;

[h l] = max(reshape(TC,m,n), [], 2);   % ����ÿһ�е����ֵ������m �Լ� �±�����l��  ע��1 Ϊÿ�У� 2 Ϊÿ�С�

model.modelname = 'LB_Inv_UFL';
model.modelsense = 'min';

ncol = (m + m * n)*2 ;

model.lb    = zeros(ncol, 1);
model.ub    = inf(ncol, 1);

obj = ones(m + m * n ,1);

model.obj = [obj; obj];   % norm-1 c-Costs  ��Ϊ��

% model.vtype = [repmat('B', nPlants, 1); repmat('C', nPlants * nplayers, 1)];
%
% Set data for constraints and matrix

nrow = m + n + length(v1) - s + 1  ; % ǰ����Լ��������Ϊ m ����������Լ���� n ����
% ���ĸ�Լ���� vi=1 ��Ӧ���� rik==1 �������� 1 �ĸ��� �� һ����ʽ

model.A     = sparse(nrow, ncol);

model.sense = [repmat('=', 1, 1); repmat('<', nrow-1, 1)];

% Production constraints   ע������������Ҫ����  ��һ��ǳ�����

model.A(1,:) = x0;

model.rhs(1) = V_UFL-V_0;

for i=1:length(v1)
    % ����Ҫ�ҵ� i �����е� rik Ϊ 1 �����ֵ
    max()