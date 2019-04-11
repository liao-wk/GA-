%%%%%%%%%%%%%%%%%%%%%%遗传算法解决TSP问题%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
City = load("City.txt"); % 导入城市坐标数据
a = load("best_so_far.txt");
city_num = size(City,1);
GATSP_funs = GA_TSP_funs; % 给函数赋句柄
%%  求任意两个城市之间的距离  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
City_dist = GATSP_funs.distance_value(city_num,City);
NP = 200; % 种群的规模，200个个体
G = 1000; % 最大遗传代数1000代
%%  生成初始种群及对应的适应度值  %%%%%%%%%%%%%%%%%%%%%%%%%%%
Pop.pop = GATSP_funs.pop_gene(NP,city_num);
Pop.pop_fit = GATSP_funs.pop_fit(city_num,Pop.pop,City_dist); 
gen = 1; % 迭代
Pop.best_so_far = Pop.pop(1,:); % 最优种群个体的初始化
Pop.best_so_far_fit = Pop.pop_fit(1);  % 最优种群个体适应度的初始化
[Pop.best_so_far,Pop.best_so_far_fit] = GATSP_funs.best_so_far_fun(Pop);
fit = zeros(G,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%遗传算法循环%%%%%%%%%%%%%%%%%%%%%%%%
while gen < G   
    fit(gen,:) = Pop.best_so_far_fit(1);
    %% 计算归一化适应值 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Pop.pop_fit_norm = GATSP_funs.normalized(Pop,NP);
    %% 选择操作，选出优良个体  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Pop.pop_fine = GATSP_funs.select_operator(NP,city_num,Pop);
    Pop.pop_fine_fit = GATSP_funs.pop_fit(city_num,Pop.pop_fine,City_dist); 
    new_pop = size(Pop.pop_fine,1); % 优良个体的数量 
    while new_pop < NP 
        %% 交叉操作,得到子代  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 交叉策略引入保留精英个体的策略
        [Chrom1,Chrom2] = GATSP_funs.cross_operator(new_pop,Pop,city_num);
        %% 变异操作，得到子代  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Chrom1 = GATSP_funs.variation_operator(Chrom1,city_num,Pop,City_dist);
        Chrom2 = GATSP_funs.variation_operator(Chrom2,city_num,Pop,City_dist);
        % 在交配池中添加子代
        Pop.pop_fine = [Pop.pop_fine;Chrom1;Chrom2];
        new_pop = size(Pop.pop_fine,1);
        % 保留优良个体
        if new_pop >= NP
            index = find(min(Pop.pop_fit));
            Pop.pop(index,:) = Pop.best_so_far;
            if index > NP
                Pop.pop_fine = Pop.pop_fine(new_pop-NP:end,:);
            else
                Pop.pop_fine = Pop.pop_fine(1:NP,:);
            end
        end
    end
    Pop.pop = Pop.pop_fine; % 更新种群
    Pop.pop_fit = GATSP_funs.pop_fit(city_num,Pop.pop,City_dist);
    [Pop.best_so_far,Pop.best_so_far_fit] = GATSP_funs.best_so_far_fun(Pop);
    clear Pop.pop_fine;
    gen = gen + 1;
end
%% 绘图
fit(gen:end,:) = []; 
figure
plot(fit);
title(["优化最短距离：",num2str(Pop.best_so_far_fit)]);
xlabel("迭代次数");
ylabel("目标函数值");


