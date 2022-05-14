clear;
clc;

%初始化种群
factory=1;   %加工车间得数量
retailer=10;  %输入看顾客的数目

%%设计10*10的矩阵存储不同顾客之间的距离矩阵
load('retailer10');
zuobiao=c1011;
h=pdist(zuobiao);
dist=squareform(h);  
realdist=zeros(11,11);
%设置标准差因子来控制随机性的大小
e=1; %%e=0.5
%%采用蒙特卡洛方法
for i=1:11
    for j=1:11
        w=dist(i,j);
        r=e*w; %%标准差的大小
        y=normrnd(w,r,[1 100]);
        realdist(i,j)=mean(y);
    end
    
end

size=[3 5 6 4 2 2 3 4 3 7];  %不同零售上对工件的尺寸大小
%生产的过程种的先后次序的的惩罚生产成本、排在【123456789】越靠前的惩罚系数越低
%不同零售商的单位时间生产惩罚成本【】


time=[4 7 6 3 5 5 3 4 6 2];  


delaytime=[2 1 2 3 1 1 3 4 2 3]; %不同订单的延期惩罚系数以及早到的库存惩罚系数
expectime=[19 28 10 18 25 5 4 10 19 15];  %不同零售商的预期到达时间

vehicle=3;  %配送车两的数目
vehicleSize=[9 6 8];  %配送的车辆的对应的最大的车辆容量
vehicleCost=[3 2 3];  %不同车辆的的单位时间运行成本
%采用双层的编码、第一层表示零售上的配送次序、第二层代表选择的配送车辆
Chrom=zeros(50,20);%预定义零矩阵，用于存数20个染色体

%k=seed(randperm(numel(seed)));
tic
NIND=200;
for i=1:NIND
    hou=zeros(1,retailer);%预定义零矩阵，
    chu=zeros(1,3);
    Chrom(i,1:retailer)=(randperm(retailer));%生成染色体第一层
   
    for j=1:retailer
        hou(j)=unidrnd(vehicle);    %生成染色体的第二层
    end

    chu(1)=unidrnd(4);
    chu(2)=unidrnd(4);
    chu(3)=retailer-chu(2)-chu(1);
    Chrom(i,11:20)=hou;
    Chrom(i,21:23)=chu;
    
end

 WNumber=retailer*2;
 XOVR=0.8;
 MUTR=0.2;
 time_opt=zeros(NIND,200);% 预定义NIND*100的矩阵存储100代种群中的各个个体时间
 Obj=zeros(200,1);
for generation=1:200
%%设计存储矩阵存储行使的车辆所走过的路过的零售商—
 
%%重新生成零售商的订单的满足的先后时间 例子 ：第九个零售上的需求先到的、到达的时间和预期到达的时间做对比时间差、
%% 第一行为不同零售商车辆的到达时间（不包含等待时间）、第二行为包含等待时间、第三行为订单实际到达零售商的时间、第四行为实际和预期差
%%第一行位不同顾客订单的生产加工完胜时间、第二行为不同订单的配送到达的时间（不包含等待生产的时间）、第三实际到达的时间、第四位时间床惩罚



%%设计矩阵存储在生产的过程中的等待时间、是否产生了等待、等待时间为该车辆的出发时间减去该零售商的订单以及之前的订单加工时间计算出来的库存的成本为该订单的等待时间*单位等待库存成本

%%以下进行解码、计算适应度、选择————————————

 T_qunti=zeros(NIND,1);
%%不同的配送次数所路过的零售商
%该函数用来1 :计算种群中各个染色体的适应度，
%          2 :通过选择生成新的种群，其中，最优个体直接保留
%  将该时间放到NIND行100矩阵中，每一列代表一代中的NIND个个体的时间
for i=1:NIND         %%群体中的各个染色体得适应度
        some=Chrom(i,:);
F=calfit(retailer,some,time,realdist,expectime,delaytime);
    T_qunti(i,1)=F;   %将该染色体对目标函数值统计

     
time_add=sum(T_qunti); %计算出种群中各个染色体总时间和
%%计算种群染色体适应度的平均值
avg=time_add/i;
time_indiv=T_qunti/time_add; %计算每一个个体与总时间的比值
min_time=min(T_qunti);%该代种群中时间最短的个体时间
max_time=max(T_qunti);%该代种群中时间最短的个体时间
end %%截至该句，计算出每个染色体的时间
% 
% 
% disp(['第几条',num2str(1)])
%        disp(Chrom(1,:))
%以下执行*选择*操作 ---.>稳态复制的方法
next_pop=Chrom;%初始化群体
best_flag=0;
  for tt=1:NIND
    if T_qunti(tt,1)==min_time %如果该个体为截止到当前代最好的个体，则保留
       best_flag=best_flag+1;
       next_pop(best_flag,:)=Chrom(tt,:);
%              disp(Chrom(tt,:));
%        disp(['适应度',num2str(min_time)]);

       time_indiv(tt,:)=2;% 如果该个体为该代中最优，则将该个体直接复制到下一代，然后将time_indiv赋值为2，以避开交叉和变异
    end
  end

flag=best_flag;
while flag<NIND
%     for z=1:20 %下一代群体的前flag个个体直接取上一带的最佳个体，剩下的个体用随即方法选择
%         sj=rand;
%         if  time_indiv(z,1)<sj %如果随机数大于第i个染色体的概率，认为该染色体较好，保留
%             next_pop(flag+1,:)=Chrom(z,:);
%             flag=flag+1;
%         
%         end
%     end
     Pos1=unidrnd(NIND);%
     Pos2=unidrnd(NIND);
%二元锦标赛选择的第几条染色体
  while Pos1==Pos2      
        Pos2=unidrnd(NIND);
  end  
  if  time_indiv(Pos1,1)<time_indiv(Pos2,1)
      next_pop(flag+1,:)=Chrom(Pos1,:);
  else
      next_pop(flag+1,:)=Chrom(Pos2,:);
  end
 flag=flag+1;
end
Chrom=next_pop;



SelNum=randperm(NIND);
 %交叉个体组个数
%  Num=NIND/2;
%  Num=2*fix(Num);
ChromNew=Chrom;
 
 for mm=(NIND-best_flag):NIND-1
     rd=rand;
             S1=Chrom(SelNum(mm),:);
             S2=Chrom(SelNum(mm+1),:);
%                         %取交换的个体;
%                         S1=Chrom(SelNum(1),:);
%                         S2=Chrom(SelNum(2),:);
                        
            if best_flag==NIND||SelNum(mm)==1 % 如果该种群最优个体达到最大，跳出交叉
                break
            end
            if T_qunti(SelNum(mm),:)>avg
                bijiao=(T_qunti(SelNum(mm),:)-avg)/(max_time-avg);
                XOVR=0.9/(0.3+exp(bijiao));
            else
                XOVR=0.8;
            end
     if XOVR>rd;
         

            

            A=S1;
            B=S2;
             %交叉点--postion
            n=fix(retailer*rand); 
            C=[A(1:retailer+n) B(n+retailer+1:2*retailer) A(21:23) ];
            D=[B(1:retailer+n) A(n+retailer+1:2*retailer) B(21:23)];

%             c1=sum(C(n+10:18)==1);
%             d1=sum(D(n+10:18)==1);
%             c2=sum(C(n+10:18)==2);
%             d2=sum(D(n+10:18)==2);
%             c3=sum(C(n+10:18)==3);
%             d3=sum(D(n+10:18)==3);
            
%             E=[ones(1,3-c1) 2*ones(1,3-c2) 3*ones(1,3-c3) C(n+1:9)];%重排C
%             F=[ones(1,3-d1) 2*ones(1,3-d2) 3*ones(1,3-d3) D(n+1:9)];%重排C
            %取出交叉位及其之前的元素放到ex1和ex2中
%             ex1=E(1,1:(9-c1-c2-c3));
%             ex2=F(1,1:(9-c1-c2-c3));
%             ex1=ex1(randperm(numel(ex1)));
%             ex2=ex2(randperm(numel(ex2)));
%             E=[ex1 E(n+1:9)];
%             F=[ex1 F(n+1:9)];
%            %放入新群
             ChromNew(SelNum(mm),:)=C;
             ChromNew(SelNum(mm+1),:)=D;
%              Chrom(SelNum(mm),:)=ChromNew(SelNum(mm),:);
%              Chrom(SelNum(mm+1),:)=ChromNew(SelNum(mm+1),:);
%             
     end
 end% % %截止到该行语句 实现了交叉，以下进行变异操作





% % %--------------------------------------------------------------------------
% % %截止到该行语句 实现了交叉，以下进行变异操作
% % %-启发式变异算子  直接更新该染色体、将染色体按照生产的完成时间排序

Chrom=ChromNew;

for i=best_flag+1:NIND  %是否变异
            if T_qunti(i)>avg
                bijiao=(T_qunti(i)-avg)/(max_time-avg);
                XOVR=0.8/(0.1+exp(bijiao));
            else
                XOVR=0.2;
            end
    mt=rand;
  if MUTR>mt;     

 
   S=Chrom(i,:);
   if  generation<50
   prTime=zeros(retailer,1);
 %此处进行两阶段启发式变异
  for x=21:23
       plant=S(x);
       ftime=0;
       for y=1:plant
%            if y==1
%                ftime=time(S(y));
%            else
%                ftime=time(S(y))+time(S(y-1));
%            end
                      if x==21;
                            ftime=ftime+time(S(y));
                             prTime(S(y),1)=ftime;
                      elseif x==22;
                          ftime=ftime+time(S(y+S(21)));
                          prTime(S(y+S(21)),1)=ftime;
                      else
                          ftime=ftime+time(S(y+retailer-S(23)));
                          prTime(S(y+retailer-S(23)),1)=ftime; 
                      end;
           
           
       end
  end
 [PT,PTdex]=sortrows(prTime); 
 for z=1:retailer
      aaaa=find(S(1:10)==PTdex(z))+10;
      if z<4
          S(aaaa)=1;
      elseif z<7
           S(aaaa)=2;
      elseif z<10
           S(aaaa)=3;
      end
      
 end

   else
     %交换
     

         Pos1=unidrnd(WNumber/2);%变异位置
    Pos2=unidrnd(WNumber/2);
     
%变异位置不相同
  while Pos1==Pos2      
        Pos2=unidrnd(WNumber/2);
  end  
 %取数据
   temp=S(Pos1);
    S(Pos1)=S(Pos2);
    S(Pos2)=temp;
     end


   ChromNew(i,:)=S;
 end
end
Chrom=ChromNew;

time_opt(:,generation)=T_qunti;%时间矩阵，用来存储N代种群中各个染色体的时间

% % %--------------------------------------------------------------------------
% % %截止到该行语句 实现变异
% % %--------------------------------------------------------------------------
%pp=T_qunti==min_time;%pp记录T_qunti中与最小时间个体相同的个体数量（为20*1的矩阵）
%%引入禁忌搜索的操作
% gt=unidrnd(20);
% while gt==1
%     gt=unidrnd(20);
% end

for gt=NIND-NIND/2:NIND

some=Chrom(gt,:);
TabuList=zeros(retailer,23);                      % (tabu list)
TabuLength=round((retailer*(retailer-1)/2)^0.5);%禁忌表长度(tabu length)
Candidates=10;                               %候选集的个数 (全部领域解个数)
CandidateNum=zeros(Candidates,23);       %候选解集合
BSF=some;  
F=calfit(retailer,BSF,time,realdist,expectime,delaytime);%best so far;
BestL=F;                                    %当前最佳解距离
p=1;                                         %记录迭代次数
StopL=10;                                  %最大迭代次数



while p<StopL
    tabufit=zeros(Candidates,2);
    for q=1:Candidates
             Pos1=unidrnd(retailer);%
     Pos2=unidrnd(retailer);
          while Pos1==Pos2      
        Pos2=unidrnd(retailer);
          end  

        %交换
        temp=BSF(Pos1);
        BSF(Pos1)=BSF(Pos2);
        BSF(Pos2)=temp;
        CandidateNum(q,:)=BSF;

        tabufit(q,1)=calfit(retailer,BSF,time,realdist,expectime,delaytime);
        tabufit(q,2)=q;
% if   tabufit(q,1)<BestL
%     BestL=tabufit(q,1);
%     biao=q;
% end
        
    end
    
    tabufit=sortrows(tabufit,1);
    if tabufit(1,1)<=BestL
            BSF=CandidateNum(tabufit(1,2),:);
    BestL=tabufit(1,1);

    else 
        flag=0;
        for f=1:6
            if TabuList(f,:)==CandidateNum(tabufit(1,2),:);
            flag=1;
            end
        end
        if  flag==1
             BSF=CandidateNum(tabufit(2,2),:);
    BestL=tabufit(2,1);
        else
                        BSF=CandidateNum(tabufit(1,2),:);
    BestL=tabufit(1,1);
        end
    end
    for d=1:5
        TabuList(d+1,:)=TabuList(d,:);
    end
    TabuList(1,:)=BSF;
    

    
        p=p+1;                                                          %迭代次数加1
end

for g=1:6
    
    fit=calfit(retailer,TabuList(g,:),time,realdist,expectime,delaytime);
    if BestL>fit
        BestL=fit;
        BSF=TabuList(g,:);
    end
end
  
Chrom(gt,:)=BSF;

end



Obj(generation)=min(T_qunti);


end

figure(1)
plot(1:generation,Obj,'-')
hold on
xlabel('迭代次数')
ylabel('适应度')
title('优化过程')


display('最优个体为：')
Chrom_best=Chrom(1,:) %显示最优个体的染色体编码，每一代进化时都将最优个体放置在群体的第一行
display('最优函数适应度值为:')
min_time
toc
% display('解码后对应最优个体的加工顺序为：')
% P_best=P
% figure(2)
% for sl=0:999
%     plot((sl+0.05):0.05:(1+sl),(time_opt(:,sl+1))','r*')
%     hold on
% end



