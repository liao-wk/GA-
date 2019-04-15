%%%%%%%%%%%%%%%%%%遗传算法引入模拟退火算法解决TSP问题%%%%%%%%%%%%%%%%%%
City = load("City.txt"); % 导入城市坐标数据
city_num = size(City,1);
GATSP_funs = GA_TSP_funs; % 给函数赋句柄
GATSP_funs2 = GA_TSP_funs2; % 给函数赋句柄
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
fit = zeros(100000,1);
%%%%%%%%%%%%%%%%%%%%%%遗传算法嵌套模拟退火算法循环SSB1算法%%%%%%%%%%%%%%%%%
T = 100; % 初始温度为100度
L = 100; % 内层循环次数
K = 0.89; % 衰减参数
p = 1;
while T > 0.01
    fit(p,:) = Pop.best_so_far_fit;
    p = p+1;
    %% 创建内循环  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:20
        %% 计算归一化适应值 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Pop.pop_fit_norm = GATSP_funs.normalized(Pop,NP);
        %% 选择操作，选出优良个体  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Pop.pop_fine = GATSP_funs.select_operator(NP,city_num,Pop);
        Pop.pop_fine_fit = GATSP_funs.pop_fit(city_num,Pop.pop_fine,City_dist); 
        Pop.pop_parent = Pop.pop; % 将当前代个体赋给父代
        Pop.pop_parent_fit = Pop.pop_fit; % 父代的适应度
        %% 交叉操作,得到子代   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Pop.pop_child = GATSP_funs2.cross_operator(Pop,NP,city_num);
        Pop.pop_child_fit = GATSP_funs.pop_fit(city_num,Pop.pop,City_dist);
        %% 变异操作  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Pop.pop_child = GATSP_funs2.variation_operator2(Pop,city_num,NP,City_dist);
        Pop.pop_child_fit = GATSP_funs.pop_fit(city_num,Pop.pop,City_dist);
        %% 模拟退火算法接受劣解策略  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Pop.pop = GATSP_funs2.Compare_parent_child(Pop,city_num,NP,T); 
        Pop.pop_fit = GATSP_funs.pop_fit(city_num,Pop.pop,City_dist); 
        [Pop.best_so_far,Pop.best_so_far_fit] = GATSP_funs.best_so_far_fun(Pop);
    end
    T = T*K;
end
%% 绘图
fit(p:end,:) = []; 
figure
plot(fit);
title(["优化最短距离：",num2str(Pop.best_so_far_fit)]);
xlabel("迭代次数");
ylabel("目标函数值");














