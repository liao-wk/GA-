function funs = GA_TSP_funs
%UNTITLED2 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
funs.distance_value = @distance_value; % ����������������֮��ľ���
funs.pop_gene = @pop_gene; % ���ɳ�ʼ��Ⱥ
funs.fit = @fitness; % ��Ӧ�Ⱥ���
funs.pop_fit = @pop_fit; % ��Ⱥ��Ӧ�Ⱥ���
funs.best_so_far_fun = @best_so_far_fun; % ������Ⱥ���弰��Ӧ�ȵĴ洢
funs.normalized = @normalized; % �����һ����Ӧֵ
funs.select_operator = @select_operator; % ѡ������
funs.cross_operator = @cross_operator; % �������
funs.variation_operator = @variation_operator; % ��������
end
%% �Ӻ���1������������������֮��ľ��� %%%%%%%%%%%%%%%%%%%%
function City_dist = distance_value(city_num,City)
City_dist = zeros(city_num,city_num);
for i = 1:city_num
    for j = 1:city_num
        City_dist(i,j) = sqrt((City(i,1)-City(j,1)).^2+(City(i,2)-...
            City(j,2)).^2);
    end
end
end
%% �Ӻ���2�����ɳ�ʼ��Ⱥ %%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pop = pop_gene(NP,city_num)
pop = zeros(NP,city_num); % ��Ⱥ����������ʼֵΪ0����
for i = 1:NP
    pop(i,:) = randperm(city_num);
end
end
%% �Ӻ���3����Ӧ�Ⱥ��� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fit = fitness(city_num,chrom,City_dist)
fit = 0;
% chrom��Ⱦɫ�����˼���������
for i = 1:city_num-1
    fit = fit + City_dist(chrom(i),chrom(i+1));
end
fit = fit + City_dist(chrom(city_num),chrom(1));
end
%% �Ӻ���4����Ⱥ��Ӧ�Ⱥ���  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pop_fit = pop_fit(city_num,pop,City_dist)
pop_size = size(pop,1);
pop_fit = zeros(pop_size,1);
for i = 1:pop_size
    pop_fit(i,1) = fitness(city_num,pop(i,:),City_dist);
end
end
%% �Ӻ���5��������Ⱥ���弰��Ӧ�ȵĴ洢 %%%%%%%%%%%%%%%%%%%%
function [best_so_far,best_so_far_fit] = best_so_far_fun(Pop)
best_so_far = Pop.best_so_far;
best_so_far_fit = Pop.best_so_far_fit;
NP = size(Pop.pop,1);
for i = 1:NP
    if best_so_far_fit > Pop.pop_fit(i,1)
        best_so_far = Pop.pop(i,:);
        best_so_far_fit = Pop.pop_fit(i,1);
    end
end
end
%% �Ӻ���6,�����һ����Ӧֵ(min->max) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fit_norm = normalized(Pop,NP)
min_fit = min(Pop.pop_fit); % ��Ⱥ����С��Ӧ��
max_fit = max(Pop.pop_fit); % ��Ⱥ�������Ӧ��
fit_norm = zeros(NP,1); % ��ʼ��
fit_norm(:,1) = 1-(Pop.pop_fit(:)-min_fit)/(max_fit-min_fit+0.001);
end

%% �Ӻ���8��ѡ������,ѡ���������� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pop_fine = select_operator(NP,city_num,Pop)
pop_fine = zeros(NP,city_num); % ���������帳�ռ�
% �������ĸ���һ��Ҫ��������
p = 1;
for i = 1:NP
    if Pop.pop_fit_norm(i,1) >= rand || Pop.pop_fit(i,1) <= Pop.best_so_far_fit 
        pop_fine(p,:) = Pop.pop(i,:);
        p = p+1;
    end
end
pop_fine(p:end,:) = []; % ɾ������Ŀռ� 
end
%% �Ӻ���9���������   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Chrom1,Chrom2] = cross_operator(new_pop,Pop,city_num)
cro_dot = randperm(new_pop); % �����
Chrom1 = Pop.pop_fine(cro_dot(1),:); % ������ɵĵ�һ������
Chrom2 = Pop.pop_fine(cro_dot(2),:); % ������ɵĵڶ�������
cro_num = ceil(city_num/10); % �����ĸ������ڳ���������10
p = unidrnd(city_num-cro_num+1); % ���ѡ�񽻲淶Χ����p ��p+cro_num
for i = 1:cro_num
    x = find(Chrom1==Chrom2(p+i-1));
    y = find(Chrom2==Chrom1(p+i-1));
    %  Ⱦɫ��������н���
    temp = Chrom1(p+i-1);
    Chrom1(p+i-1) = Chrom2(p+i-1);
    Chrom2(p+i-1) = temp;
    % ��ÿ��Ⱦɫ������в�ȫ
    temp = Chrom1(x);
    Chrom1(x) = Chrom2(y);
    Chrom2(y) = temp;
end
end
%% �Ӻ���10���������ӣ��õ��µĸ��� %%%%%%%%%%%%%%%%%%%%%%%%%%
% ������ԣ��������䡣
function Chrom = variation_operator(Chrom,city_num,Pop,City_dist)
fit = fitness(city_num,Chrom,City_dist); % Ⱦɫ���Ӧ����Ӧ��
min_fit = min(Pop.pop_fine_fit); % ��Ⱥ����С��Ӧ��
max_fit = max(Pop.pop_fine_fit); % ��Ⱥ�������Ӧ��
prob = 1-(fit-min_fit)/(max_fit-min_fit+0.001); % �������,��Ӧ��Խ�ͣ�probԽ�ӽ���1.
if prob > 0.8
    prob = 0.01;
elseif prob > 0.6
    prob = 0.05;
elseif prob > 0.3
    prob = 0.1;
else
    prob = 0.2;
end
for i = 1:city_num
    if rand < prob
        r = randperm(city_num);
        temp = Chrom(i);
        Chrom(i) = Chrom(r(1));
        Chrom(r(1)) = temp;
    end
end
end











