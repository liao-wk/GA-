%%%%%%%%%%%%%%%%%%%%%%�Ŵ��㷨���TSP����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
City = load("City.txt"); % ���������������
a = load("best_so_far.txt");
city_num = size(City,1);
GATSP_funs = GA_TSP_funs; % �����������
GATSP_funs2 = GA_TSP_funs2; % �����������
%%  ��������������֮��ľ���  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
City_dist = GATSP_funs.distance_value(city_num,City);
NP = 200; % ��Ⱥ�Ĺ�ģ��200������
G = 500; % ����Ŵ�����1000��
%%  ���ɳ�ʼ��Ⱥ����Ӧ����Ӧ��ֵ  %%%%%%%%%%%%%%%%%%%%%%%%%%%
Pop.pop = GATSP_funs.pop_gene(NP,city_num);
Pop.pop_fit = GATSP_funs.pop_fit(city_num,Pop.pop,City_dist); 
gen = 1; % ����
Pop.best_so_far = Pop.pop(1,:); % ������Ⱥ����ĳ�ʼ��
Pop.best_so_far_fit = Pop.pop_fit(1);  % ������Ⱥ������Ӧ�ȵĳ�ʼ��
[Pop.best_so_far,Pop.best_so_far_fit] = GATSP_funs.best_so_far_fun(Pop);
fit = zeros(G,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%�Ŵ��㷨ѭ��%%%%%%%%%%%%%%%%%%%%%%%%
while gen < G
    fit(gen,:) = Pop.best_so_far_fit(1);
    %% �����һ����Ӧֵ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Pop.pop_fit_norm = GATSP_funs.normalized(Pop,NP);
    %% ѡ�������ѡ����������  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Pop.pop_fine = GATSP_funs.select_operator(NP,city_num,Pop);
    Pop.pop_fine_fit = GATSP_funs.pop_fit(city_num,Pop.pop_fine,City_dist); 
    %% �������   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Pop.pop = GATSP_funs2.cross_operator(Pop,NP,city_num);
    Pop.pop_fit = GATSP_funs.pop_fit(city_num,Pop.pop,City_dist);
    %% �����һ����Ӧֵ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Pop.pop_fit_norm = GATSP_funs.normalized(Pop,NP);
    %% ѡ�������ѡ����������  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Pop.pop_fine = GATSP_funs.select_operator(NP,city_num,Pop);
    Pop.pop_fine_fit = GATSP_funs.pop_fit(city_num,Pop.pop_fine,City_dist); 
    %% �������  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Pop.pop = GATSP_funs2.variation_operator(Pop,city_num,NP,City_dist);
    Pop.pop_fit = GATSP_funs.pop_fit(city_num,Pop.pop,City_dist);
    [Pop.best_so_far,Pop.best_so_far_fit] = GATSP_funs.best_so_far_fun(Pop);
    gen = gen+1;
end
%% ��ͼ
fit(gen:end,:) = []; 
figure
plot(fit);
title(["�Ż���̾��룺",num2str(Pop.best_so_far_fit)]);
xlabel("��������");
ylabel("Ŀ�꺯��ֵ");














