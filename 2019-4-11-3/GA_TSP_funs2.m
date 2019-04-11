function funs = GA_TSP_funs2
funs.pop_fit = @pop_fit; % 种群适应度
funs.cross_operator = @cross_operator; % 交叉算子
funs.variation_operator = @variation_operator; % 变异算子
end
%% 子函数1，适应度函数 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fit = fitness(city_num,chrom,City_dist)
fit = 0;
% chrom是染色体的意思，代表个体
for i = 1:city_num-1
    fit = fit + City_dist(chrom(i),chrom(i+1));
end
fit = fit + City_dist(chrom(city_num),chrom(1));
end
%% 子函数2，种群适应度函数  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pop_fit = pop_fit(city_num,pop,City_dist)
pop_size = size(pop,1);
pop_fit = zeros(pop_size,1);
for i = 1:pop_size
    pop_fit(i,1) = fitness(city_num,pop(i,:),City_dist);
end
end
%% 子函数3，交叉算子 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pop = cross_operator(Pop,NP,city_num)
% 交叉按照一定的概率 Pc进行交叉，
Pc = 0.85; % 交叉概率
pop_cross_size = 2*NP;
% 保留精英种群中的精英个体，前10%的精英个体
pool_mating = zeros(NP,city_num);  % 交配池
pop_fine_fit = sort(Pop.pop_fine_fit);
pop_size = size(Pop.pop_fine,1); % 精英种群的个数
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
% 对精英个体按照一定的交叉概率进行交叉
while new_pop <= pop_cross_size
    if rand < Pc
        cro_dot = randperm(pop_size); % 交叉点
        Chrom1 = Pop.pop_fine(cro_dot(1),:); % 随机生成的第1个个体
        Chrom2 = Pop.pop_fine(cro_dot(2),:); % 随机生成的第2个个体
        cro_num = ceil(city_num/10); % 交叉点的个数等于城市数除以10
        p = unidrnd(city_num-cro_num+1); % 随机选择交叉范围，从p 到p+cro_num
        for j = 1:cro_num
            x = find(Chrom1==Chrom2(p+j-1));
            y = find(Chrom2==Chrom1(p+j-1));
            %  染色体基因序列交换,生出了两个新的个体
            temp = Chrom1(p+j-1);
            Chrom1(p+j-1) = Chrom2(p+j-1);
            Chrom2(p+j-1) = temp;
            % 将每个染色体的序列补全
            temp = Chrom1(x);
            Chrom1(x) = Chrom2(y);
            Chrom2(y) = temp;
        end
        pool_mating = [pool_mating;Chrom1;Chrom2]; % 填充交配池
        new_pop = size(pool_mating,1);
        if new_pop > pop_cross_size 
            pool_mating = pool_mating(1:pop_cross_size,:);
        end
    end
end
% 从交配池中选出NP个优良个体作为新的种群，这样会不会造成早衰？？？
%  经过试验，容易造成早衰
pop = pool_mating(pop_cross_size-NP+1:pop_cross_size,:); % 得到新的种群
end
%% 子函数4，变异算子  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pop = variation_operator(Pop,city_num,NP,City_dist)
pool_mut = Pop.pop; % 变异池
% 变异概率
pop_size = size(Pop.pop,1);
for i = 1:pop_size
    Chrom = Pop.pop(i,:);
    fit = fitness(city_num,Chrom,City_dist); % 染色体对应的适应度
    min_fit = min(Pop.pop_fit); % 种群中最小适应度
    max_fit = max(Pop.pop_fit); % 种群中最大适应度
    prob = 1-(fit-min_fit)/(max_fit-min_fit+0.001); % 变异概率,适应度越低，prob越接近于1.
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

% 从变异池中选出NP个优良个体作为新的种群，这样会不会造成早衰？？？
%  经过试验，容易造成早衰
pop = pool_mut(1:NP,:);  % 得到新的种群
end















