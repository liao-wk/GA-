function funs = GA_TSP_funs2
funs.pop_fit = @pop_fit; % ��Ⱥ��Ӧ��
funs.cross_operator = @cross_operator; % ��������
funs.variation_operator = @variation_operator; % ��������
funs.variation_operator2 = @variation_operator2; % ��������2
funs.Simliar = @Simliar; % ��������Ⱦɫ��֮������ƶ�
funs.Compare_parent_child = @Compare_parent_child; % �Ƚϸ������Ӵ�
end
%% �Ӻ���1����Ӧ�Ⱥ��� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fit = fitness(city_num,chrom,City_dist)
fit = 0;
% chrom��Ⱦɫ�����˼���������
for i = 1:city_num-1
    fit = fit + City_dist(chrom(i),chrom(i+1));
end
fit = fit + City_dist(chrom(city_num),chrom(1));
end
%% �Ӻ���2����Ⱥ��Ӧ�Ⱥ���  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pop_fit = pop_fit(city_num,pop,City_dist)
pop_size = size(pop,1);
pop_fit = zeros(pop_size,1);
for i = 1:pop_size
    pop_fit(i,1) = fitness(city_num,pop(i,:),City_dist);
end
end
%% �Ӻ���3���������� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pop = cross_operator(Pop,NP,city_num)
% ���水��һ���ĸ��� Pc���н��棬
Pc = 0.85; % �������
pop_cross_size = 2*NP;
% ������Ӣ��Ⱥ�еľ�Ӣ���壬ǰ10%�ľ�Ӣ����
pool_mating = zeros(NP,city_num);  % �����
pop_fine_fit = sort(Pop.pop_fine_fit);
pop_size = size(Pop.pop_fine,1); % ��Ӣ��Ⱥ�ĸ���
fit_suspend = pop_fine_fit(ceil(pop_size*0.6),1);
q = 1;
for i = 1:pop_size
    if Pop.pop_fine_fit(i,1) <= fit_suspend
        pool_mating(q,:) = Pop.pop_fine(i,:);
        q = q+1;
    end
end
pool_mating(q:end,:) = [];
new_pop = size(pool_mating,1); 
% �Ծ�Ӣ���尴��һ���Ľ�����ʽ��н���
while new_pop <= pop_cross_size
    if rand < Pc
        cro_dot = randperm(pop_size); % �����
        Chrom1 = Pop.pop_fine(cro_dot(1),:); % ������ɵĵ�1������
        Chrom2 = Pop.pop_fine(cro_dot(2),:); % ������ɵĵ�2������
        cro_num = ceil(city_num/10); % �����ĸ������ڳ���������10
        p = unidrnd(city_num-cro_num+1); % ���ѡ�񽻲淶Χ����p ��p+cro_num
        for j = 1:cro_num
            x = find(Chrom1==Chrom2(p+j-1));
            y = find(Chrom2==Chrom1(p+j-1));
            %  Ⱦɫ��������н���,�����������µĸ���
            temp = Chrom1(p+j-1);
            Chrom1(p+j-1) = Chrom2(p+j-1);
            Chrom2(p+j-1) = temp;
            % ��ÿ��Ⱦɫ������в�ȫ
            temp = Chrom1(x);
            Chrom1(x) = Chrom2(y);
            Chrom2(y) = temp;
        end
        pool_mating = [pool_mating;Chrom1;Chrom2]; % ��佻���
        new_pop = size(pool_mating,1);
        if new_pop > pop_cross_size 
            pool_mating = pool_mating(1:pop_cross_size,:);
        end
    end
end
% �ӽ������ѡ��NP������������Ϊ�µ���Ⱥ�������᲻�������˥������
%  �������飬���������˥
pop = pool_mating(pop_cross_size-NP+1:pop_cross_size,:); % �õ��µ���Ⱥ
end
%% �Ӻ���4����������  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pop = variation_operator(Pop,city_num,NP,City_dist)
pool_mut = Pop.pop; % �����
% �������
pop_size = size(Pop.pop,1);
for i = 1:pop_size
    Chrom = Pop.pop(i,:);
    fit = fitness(city_num,Chrom,City_dist); % Ⱦɫ���Ӧ����Ӧ��
    min_fit = min(Pop.pop_fit); % ��Ⱥ����С��Ӧ��
    max_fit = max(Pop.pop_fit); % ��Ⱥ�������Ӧ��
    prob = 1-(fit-min_fit)/(max_fit-min_fit+0.001); % �������,��Ӧ��Խ�ͣ�probԽ�ӽ���1.
    if prob > 0.8
        prob = 0.01;
    elseif prob > 0.6
        prob = 0.05;
    elseif prob > 0.3
        prob = 0.1;
    else
        prob = 0.8;
    end
    for j = 1:city_num
        if rand < prob
            r = randperm(city_num);
            temp = Chrom(j);
            Chrom(j) = Chrom(r(1));
            Chrom(r(1)) = temp;
        end
    end
    pool_mut = [pool_mut;Chrom];
end

% �ӱ������ѡ��NP������������Ϊ�µ���Ⱥ�������᲻�������˥������
%  �������飬���������˥
pop = pool_mut(1:NP,:);  % �õ��µ���Ⱥ
end
%% �Ӻ���5����������  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pop = variation_operator2(Pop,city_num,NP,City_dist)
pool_mut = Pop.pop_child; % �����
% �������
pop_size = size(Pop.pop_child,1);
for i = 1:pop_size
    Chrom = Pop.pop_child(i,:);
    fit = fitness(city_num,Chrom,City_dist); % Ⱦɫ���Ӧ����Ӧ��
    min_fit = min(Pop.pop_child_fit); % ��Ⱥ����С��Ӧ��
    max_fit = max(Pop.pop_child_fit); % ��Ⱥ�������Ӧ��
    prob = 1-(fit-min_fit)/(max_fit-min_fit+0.001); % �������,��Ӧ��Խ�ͣ�probԽ�ӽ���1.
    if prob > 0.8
        prob = 0.01;
    elseif prob > 0.6
        prob = 0.05;
    elseif prob > 0.3
        prob = 0.1;
    else
        prob = 0.8;
    end
    for j = 1:city_num
        if rand < prob
            r = randperm(city_num);
            temp = Chrom(j);
            Chrom(j) = Chrom(r(1));
            Chrom(r(1)) = temp;
        end
    end
    pool_mut = [pool_mut;Chrom];
end

% �ӱ������ѡ��NP������������Ϊ�µ���Ⱥ�������᲻�������˥������
%  �������飬���������˥
pop = pool_mut(1:NP,:);  % �õ��µ���Ⱥ
end
%% �Ӻ���6����������Ⱦɫ��֮������ƶ� %%%%%%%%%%%%%%%%%%%%%%%%%%%
function sim = Simliar(Chrom1,Chrom2)
num = size(Chrom1);
sim = 0; % ��ʼ���ƶ�Ϊ0
for i = 1:num
    if Chrom1(i) == Chrom2(i)
        sim = sim + 1; % ��ͬλ�õĳ��б��һ�£�����1
    end
end
end
%% �Ӻ���7���Ƚϸ������Ӵ�  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pop = Compare_parent_child(Pop,city_num,NP,T)
pop_size_child = size(Pop.pop_child,1); % �Ӵ������ĸ���
pop_size_parent = size(Pop.pop_parent,1); % ���������ĸ���
pop = zeros(pop_size_child,city_num); 
p = 1; % ������ɵ��Ӵ������
% �ҵ����Ӵ���������ߵĸ���
for i = 1:pop_size_child
    similar_num = 0;
    similar_pop = Pop.pop_child(i,:); % ���Ӵ����Ƶĸ���
    similar_fit = 0;
    for j = 1:pop_size_parent
%         Chrom_parent = Pop.pop_parent(j,:);
        similar_num2 = Simliar(Pop.pop_child(i,:),Pop.pop_parent(j,:));
        if similar_num2 > similar_num
            similar_num = similar_num2;
            similar_pop = Pop.pop_parent(j,:);
            similar_fit = Pop.pop_parent_fit(j,:);
        end
    end
    if Pop.pop_child_fit(i,1) < similar_fit
        pop(p,:) = Pop.pop_child(i,:);
    else
        if exp(-1*abs(similar_fit-Pop.pop_child_fit(i,1))/T) > rand
            pop(p,:) = Pop.pop_child(i,:);
        else
            pop(p,:) = similar_pop;
        end
    end
    p = p+1;
end
pop(p:end,:) = [];
% ɾ���������Ӵ����ظ��ĸ��壬���������һЩ���岹��
pop = unique(pop,'rows');
new_pop_size = size(pop,1);
if new_pop_size < NP
    pop_size2 = NP - new_pop_size;
    for k = 1:pop_size2
        pop(end+1,:) = randperm(city_num);
    end
end
new_pop_size = size(pop,1);
if new_pop_size < NP
    disp("�Ӵ�����Ⱥ����С��NP��");
end
end
















