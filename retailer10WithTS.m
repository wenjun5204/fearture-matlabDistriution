clear;
clc;

%��ʼ����Ⱥ
factory=1;   %�ӹ����������
retailer=10;  %���뿴�˿͵���Ŀ

%%���10*10�ľ���洢��ͬ�˿�֮��ľ������
load('retailer10');
zuobiao=c1011;
h=pdist(zuobiao);
dist=squareform(h);  
realdist=zeros(11,11);
%���ñ�׼����������������ԵĴ�С
e=1; %%e=0.5
%%�������ؿ��巽��
for i=1:11
    for j=1:11
        w=dist(i,j);
        r=e*w; %%��׼��Ĵ�С
        y=normrnd(w,r,[1 100]);
        realdist(i,j)=mean(y);
    end
    
end

size=[3 5 6 4 2 2 3 4 3 7];  %��ͬ�����϶Թ����ĳߴ��С
%�����Ĺ����ֵ��Ⱥ����ĵĳͷ������ɱ������ڡ�123456789��Խ��ǰ�ĳͷ�ϵ��Խ��
%��ͬ�����̵ĵ�λʱ�������ͷ��ɱ�����


time=[4 7 6 3 5 5 3 4 6 2];  


delaytime=[2 1 2 3 1 1 3 4 2 3]; %��ͬ���������ڳͷ�ϵ���Լ��絽�Ŀ��ͷ�ϵ��
expectime=[19 28 10 18 25 5 4 10 19 15];  %��ͬ�����̵�Ԥ�ڵ���ʱ��

vehicle=3;  %���ͳ�������Ŀ
vehicleSize=[9 6 8];  %���͵ĳ����Ķ�Ӧ�����ĳ�������
vehicleCost=[3 2 3];  %��ͬ�����ĵĵ�λʱ�����гɱ�
%����˫��ı��롢��һ���ʾ�����ϵ����ʹ��򡢵ڶ������ѡ������ͳ���
Chrom=zeros(50,20);%Ԥ������������ڴ���20��Ⱦɫ��

%k=seed(randperm(numel(seed)));
tic
NIND=200;
for i=1:NIND
    hou=zeros(1,retailer);%Ԥ���������
    chu=zeros(1,3);
    Chrom(i,1:retailer)=(randperm(retailer));%����Ⱦɫ���һ��
   
    for j=1:retailer
        hou(j)=unidrnd(vehicle);    %����Ⱦɫ��ĵڶ���
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
 time_opt=zeros(NIND,200);% Ԥ����NIND*100�ľ���洢100����Ⱥ�еĸ�������ʱ��
 Obj=zeros(200,1);
for generation=1:200
%%��ƴ洢����洢��ʹ�ĳ������߹���·���������̡�
 
%%�������������̵Ķ�����������Ⱥ�ʱ�� ���� ���ھŸ������ϵ������ȵ��ġ������ʱ���Ԥ�ڵ����ʱ�����Ա�ʱ��
%% ��һ��Ϊ��ͬ�����̳����ĵ���ʱ�䣨�������ȴ�ʱ�䣩���ڶ���Ϊ�����ȴ�ʱ�䡢������Ϊ����ʵ�ʵ��������̵�ʱ�䡢������Ϊʵ�ʺ�Ԥ�ڲ�
%%��һ��λ��ͬ�˿Ͷ����������ӹ���ʤʱ�䡢�ڶ���Ϊ��ͬ���������͵����ʱ�䣨�������ȴ�������ʱ�䣩������ʵ�ʵ����ʱ�䡢����λʱ�䴲�ͷ�



%%��ƾ���洢�������Ĺ����еĵȴ�ʱ�䡢�Ƿ�����˵ȴ����ȴ�ʱ��Ϊ�ó����ĳ���ʱ���ȥ�������̵Ķ����Լ�֮ǰ�Ķ����ӹ�ʱ���������Ŀ��ĳɱ�Ϊ�ö����ĵȴ�ʱ��*��λ�ȴ����ɱ�

%%���½��н��롢������Ӧ�ȡ�ѡ�񡪡���������������������

 T_qunti=zeros(NIND,1);
%%��ͬ�����ʹ�����·����������
%�ú�������1 :������Ⱥ�и���Ⱦɫ�����Ӧ�ȣ�
%          2 :ͨ��ѡ�������µ���Ⱥ�����У����Ÿ���ֱ�ӱ���
%  ����ʱ��ŵ�NIND��100�����У�ÿһ�д���һ���е�NIND�������ʱ��
for i=1:NIND         %%Ⱥ���еĸ���Ⱦɫ�����Ӧ��
        some=Chrom(i,:);
F=calfit(retailer,some,time,realdist,expectime,delaytime);
    T_qunti(i,1)=F;   %����Ⱦɫ���Ŀ�꺯��ֵͳ��

     
time_add=sum(T_qunti); %�������Ⱥ�и���Ⱦɫ����ʱ���
%%������ȺȾɫ����Ӧ�ȵ�ƽ��ֵ
avg=time_add/i;
time_indiv=T_qunti/time_add; %����ÿһ����������ʱ��ı�ֵ
min_time=min(T_qunti);%�ô���Ⱥ��ʱ����̵ĸ���ʱ��
max_time=max(T_qunti);%�ô���Ⱥ��ʱ����̵ĸ���ʱ��
end %%�����þ䣬�����ÿ��Ⱦɫ���ʱ��
% 
% 
% disp(['�ڼ���',num2str(1)])
%        disp(Chrom(1,:))
%����ִ��*ѡ��*���� ---.>��̬���Ƶķ���
next_pop=Chrom;%��ʼ��Ⱥ��
best_flag=0;
  for tt=1:NIND
    if T_qunti(tt,1)==min_time %����ø���Ϊ��ֹ����ǰ����õĸ��壬����
       best_flag=best_flag+1;
       next_pop(best_flag,:)=Chrom(tt,:);
%              disp(Chrom(tt,:));
%        disp(['��Ӧ��',num2str(min_time)]);

       time_indiv(tt,:)=2;% ����ø���Ϊ�ô������ţ��򽫸ø���ֱ�Ӹ��Ƶ���һ����Ȼ��time_indiv��ֵΪ2���Աܿ�����ͱ���
    end
  end

flag=best_flag;
while flag<NIND
%     for z=1:20 %��һ��Ⱥ���ǰflag������ֱ��ȡ��һ������Ѹ��壬ʣ�µĸ������漴����ѡ��
%         sj=rand;
%         if  time_indiv(z,1)<sj %�����������ڵ�i��Ⱦɫ��ĸ��ʣ���Ϊ��Ⱦɫ��Ϻã�����
%             next_pop(flag+1,:)=Chrom(z,:);
%             flag=flag+1;
%         
%         end
%     end
     Pos1=unidrnd(NIND);%
     Pos2=unidrnd(NIND);
%��Ԫ������ѡ��ĵڼ���Ⱦɫ��
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
 %������������
%  Num=NIND/2;
%  Num=2*fix(Num);
ChromNew=Chrom;
 
 for mm=(NIND-best_flag):NIND-1
     rd=rand;
             S1=Chrom(SelNum(mm),:);
             S2=Chrom(SelNum(mm+1),:);
%                         %ȡ�����ĸ���;
%                         S1=Chrom(SelNum(1),:);
%                         S2=Chrom(SelNum(2),:);
                        
            if best_flag==NIND||SelNum(mm)==1 % �������Ⱥ���Ÿ���ﵽ�����������
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
             %�����--postion
            n=fix(retailer*rand); 
            C=[A(1:retailer+n) B(n+retailer+1:2*retailer) A(21:23) ];
            D=[B(1:retailer+n) A(n+retailer+1:2*retailer) B(21:23)];

%             c1=sum(C(n+10:18)==1);
%             d1=sum(D(n+10:18)==1);
%             c2=sum(C(n+10:18)==2);
%             d2=sum(D(n+10:18)==2);
%             c3=sum(C(n+10:18)==3);
%             d3=sum(D(n+10:18)==3);
            
%             E=[ones(1,3-c1) 2*ones(1,3-c2) 3*ones(1,3-c3) C(n+1:9)];%����C
%             F=[ones(1,3-d1) 2*ones(1,3-d2) 3*ones(1,3-d3) D(n+1:9)];%����C
            %ȡ������λ����֮ǰ��Ԫ�طŵ�ex1��ex2��
%             ex1=E(1,1:(9-c1-c2-c3));
%             ex2=F(1,1:(9-c1-c2-c3));
%             ex1=ex1(randperm(numel(ex1)));
%             ex2=ex2(randperm(numel(ex2)));
%             E=[ex1 E(n+1:9)];
%             F=[ex1 F(n+1:9)];
%            %������Ⱥ
             ChromNew(SelNum(mm),:)=C;
             ChromNew(SelNum(mm+1),:)=D;
%              Chrom(SelNum(mm),:)=ChromNew(SelNum(mm),:);
%              Chrom(SelNum(mm+1),:)=ChromNew(SelNum(mm+1),:);
%             
     end
 end% % %��ֹ��������� ʵ���˽��棬���½��б������





% % %--------------------------------------------------------------------------
% % %��ֹ��������� ʵ���˽��棬���½��б������
% % %-����ʽ��������  ֱ�Ӹ��¸�Ⱦɫ�塢��Ⱦɫ�尴�����������ʱ������

Chrom=ChromNew;

for i=best_flag+1:NIND  %�Ƿ����
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
 %�˴��������׶�����ʽ����
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
     %����
     

         Pos1=unidrnd(WNumber/2);%����λ��
    Pos2=unidrnd(WNumber/2);
     
%����λ�ò���ͬ
  while Pos1==Pos2      
        Pos2=unidrnd(WNumber/2);
  end  
 %ȡ����
   temp=S(Pos1);
    S(Pos1)=S(Pos2);
    S(Pos2)=temp;
     end


   ChromNew(i,:)=S;
 end
end
Chrom=ChromNew;

time_opt(:,generation)=T_qunti;%ʱ����������洢N����Ⱥ�и���Ⱦɫ���ʱ��

% % %--------------------------------------------------------------------------
% % %��ֹ��������� ʵ�ֱ���
% % %--------------------------------------------------------------------------
%pp=T_qunti==min_time;%pp��¼T_qunti������Сʱ�������ͬ�ĸ���������Ϊ20*1�ľ���
%%������������Ĳ���
% gt=unidrnd(20);
% while gt==1
%     gt=unidrnd(20);
% end

for gt=NIND-NIND/2:NIND

some=Chrom(gt,:);
TabuList=zeros(retailer,23);                      % (tabu list)
TabuLength=round((retailer*(retailer-1)/2)^0.5);%���ɱ���(tabu length)
Candidates=10;                               %��ѡ���ĸ��� (ȫ����������)
CandidateNum=zeros(Candidates,23);       %��ѡ�⼯��
BSF=some;  
F=calfit(retailer,BSF,time,realdist,expectime,delaytime);%best so far;
BestL=F;                                    %��ǰ��ѽ����
p=1;                                         %��¼��������
StopL=10;                                  %����������



while p<StopL
    tabufit=zeros(Candidates,2);
    for q=1:Candidates
             Pos1=unidrnd(retailer);%
     Pos2=unidrnd(retailer);
          while Pos1==Pos2      
        Pos2=unidrnd(retailer);
          end  

        %����
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
    

    
        p=p+1;                                                          %����������1
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
xlabel('��������')
ylabel('��Ӧ��')
title('�Ż�����')


display('���Ÿ���Ϊ��')
Chrom_best=Chrom(1,:) %��ʾ���Ÿ����Ⱦɫ����룬ÿһ������ʱ�������Ÿ��������Ⱥ��ĵ�һ��
display('���ź�����Ӧ��ֵΪ:')
min_time
toc
% display('������Ӧ���Ÿ���ļӹ�˳��Ϊ��')
% P_best=P
% figure(2)
% for sl=0:999
%     plot((sl+0.05):0.05:(1+sl),(time_opt(:,sl+1))','r*')
%     hold on
% end



