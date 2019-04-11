function funs = GA_TSP_funs
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
funs.distance_value = @distance_value; % 计算任意两个城市之间的距离
funs.pop_gene = @pop_gene; % 生成初始种群
funs.fit = @fitness; % 适应度函数
funs.pop_fit = @pop_fit; % 种群适应度函数
funs.best_so_far_fun = @best_so_far_fun; % 最优种群个体及适应度的存储
funs.normalized = @normalized; % 计算归一化适应值
funs.select_operator = @select_operator; % 选择算子
funs.cross_operator = @cross_operator; % 交叉操作
funs.variation_operator = @variation_operator; % 变异算子
end
%% 子函数1，计算任意两个城市之间的距离 %%%%%%%%%%%%%%%%%%%%
function City_dist = distance_value(city_num,City)
City_dist = zeros(city_num,city_num);
for i = 1:city_num
    for j = 1:city_num
        City_dist(i,j) = sqrt((City(i,1)-City(j,1)).^2+(City(i,2)-...
            City(j,2)).^2);
    end
end
end
%% 子函数2，生成初始种群 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pop = pop_gene(NP,city_num)
pop = zeros(NP,city_num); % 种群解向量，初始值为0向量
for i = 1:NP
    pop(i,:) = randperm(city_num);
end
end
%% 子函数3，适应度函数 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fit = fitness(city_num,chrom,City_dist)
fit = 0;
% chrom是染色体的意思，代表个体
for i = 1:city_num-1
    fit = fit + City_dist(chrom(i),chrom(i+1));
end
fit = fit + City_dist(chrom(city_num),chrom(1));
end
%% 子函数4，种群适应度函数  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pop_fit = pop_fit(city_num,pop,City_dist)
pop_size = size(pop,1);
pop_fit = zeros(pop_size,1);
for i = 1:pop_size
    pop_fit(i,1) = fitness(city_num,pop(i,:),City_dist);
end
end
%% 子函数5，最优种群个体及适应度的存储 %%%%%%%%%%%%%%%%%%%%
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
%% 子函数6,计算归一化适应值(min->max) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fit_norm = normalized(Pop,NP)
min_fit = min(Pop.pop_fit); % 种群中最小适应度
max_fit = max(Pop.pop_fit); % 种群中最大适应度
fit_norm = zeros(NP,1); % 初始化
fit_norm(:,1) = 1-(Pop.pop_fit(:)-min_fit)/(max_fit-min_fit+0.001);
end

%% 子函数8，选择算子,选出优良个体 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pop_fine = select_operator(NP,city_num,Pop)
pop_fine = zeros(NP,city_num); % 对优良个体赋空间
% 最优良的个体一定要保存下来
p = 1;
for i = 1:NP
    if Pop.pop_fit_norm(i,1) >= rand || Pop.pop_fit(i,1) <= Pop.best_so_far_fit 
        pop_fine(p,:) = Pop.pop(i,:);
        p = p+1;
    end
end
pop_fine(p:end,:) = []; % 删除多余的空间 
end
%% 子函数9，交叉操作   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Chrom1,Chrom2] = cross_operator(new_pop,Pop,city_num)
cro_dot = randperm(new_pop); % 交叉点
Chrom1 = Pop.pop_fine(cro_dot(1),:); % 随机生成的第一个个体
Chrom2 = Pop.pop_fine(cro_dot(2),:); % 随机生成的第二个个体
cro_num = ceil(city_num/10); % 交叉点的个数等于城市数除以10
p = unidrnd(city_num-cro_num+1); % 随机选择交叉范围，从p 到p+cro_num
for i = 1:cro_num
    x = find(Chrom1==Chrom2(p+i-1));
    y = find(Chrom2==Chrom1(p+i-1));
    %  染色体基因序列交换
    temp = Chrom1(p+i-1);
    Chrom1(p+i-1) = Chrom2(p+i-1);
    Chrom2(p+i-1) = temp;
    % 将每个染色体的序列补全
    temp = Chrom1(x);
    Chrom1(x) = Chrom2(y);
    Chrom2(y) = temp;
end
end
%% 子函数10，变异算子，得到新的个体 %%%%%%%%%%%%%%%%%%%%%%%%%%
% 引入策略，穷则生变。
function Chrom = variation_operator(Chrom,city_num,Pop,City_dist)
fit = fitness(city_num,Chrom,City_dist); % 染色体对应的适应度
min_fit = min(Pop.pop_fine_fit); % 种群中最小适应度
max_fit = max(Pop.pop_fine_fit); % 种群中最大适应度
prob = 1-(fit-min_fit)/(max_fit-min_fit+0.001); % 变异概率,适应度越低，prob越接近于1.
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











