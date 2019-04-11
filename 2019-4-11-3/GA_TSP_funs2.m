function funs = GA_TSP_funs2
funs.pop_fit = @pop_fit; % ��Ⱥ��Ӧ��
funs.cross_operator = @cross_operator; % ��������
funs.variation_operator = @variation_operator; % ��������
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















