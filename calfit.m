function F=calfit(retailer,some,time,realdist,expectime,delaytime)


sumvhicous=0;
vhicous=[3,0.8,5,0.6,9];
for shc=1:3
    sumvhicous=sumvhicous+vhicous(shc);
end

sumrawcous=0;
rawcous=[0.3,0.2,0.5,0.6,0.9,0.1,0.3,0.2,0.7,0.2];
for shc=1:10
    sumrawcous=sumrawcous+rawcous(shc);
end

cousatsa=zeros(3,10);

cousatsa(1,:)=[70 64 90 65 59 111 70 80 72 107];
cousatsa(2,:)=[87 68 80 56 76 89 79 76 62 113];
cousatsa(3,:)=[80 62 65 66 71 121 84 84 68 88];


arriveT=zeros(5,retailer); 
   s1=[];   s2=[];   s3=[];  
       
      proT=0;

            for a=1:3
                s=some(2*retailer+a);
             c=0;
            for b=1:s
                      if a==1;
                            c=c+time(some(b));
                      else
                          c=c+time(some(b+some(a-1+retailer)));
                      end;
                if a==1;
                          arriveT(1,some(b))=c;
                     else if a==2
                              arriveT(1,some(b+some(2*retailer+1)))=c;
                          else
                              arriveT(1,some(b+retailer-some(23)))=c;
                          end
                    
                   
                end
            end
            end
   proT=max(arriveT(1,:));

   
    chenum=zeros(1,3);
    for j=retailer+1:2*retailer
        if some(j)==1
             chenum(1)=chenum(1)+1;
                    s1(chenum(1))=some(j-retailer);
        elseif some(j)==2
             chenum(2)=chenum(2)+1;
                    s2(chenum(2))=some(j-retailer);
        elseif some(j)==3
             chenum(3)=chenum(3)+1;
                    s3(chenum(3))=some(j-retailer);
        end
    end  
    

ss(1,1:chenum(1))=s1;
ss(2,1:chenum(2))=s2;
ss(3,1:chenum(3))=s3;    
    
for o=1:3
    if chenum(o)>3   %某列车带的工件数量大于3的情况
        if chenum(o)==7
            arriveT(5,ss(o,3))=1;
            arriveT(5,ss(o,6))=1;
            arriveT(5,ss(o,7))=1;
                                yx=0;
         for r=1:3
            if r==1
                yx=yx++realdist(ss(o,r),11);
            else
                yx=yx++realdist(ss(o,r),ss(o,r-1));
            end
            arriveT(2,ss(o,r))=yx;
            
         end
         
        yx=0;
          for r=4:6
            if r==4
                yx=yx++realdist(ss(o,r),10);
            else
                yx=yx++realdist(ss(o,r),ss(o,r-1));
            end
            arriveT(2,ss(o,r))=yx;
            
          end
          arriveT(2,ss(o,7))=realdist(ss(o,7),10);
         
        else
             arriveT(5,ss(o,3))=1;
            arriveT(5,ss(o,chenum(o)))=1;
                    yx=0;
         for r=1:3
            if r==1
                yx=yx++realdist(ss(o,r),retailer+1);
            else
                yx=yx++realdist(ss(o,r),ss(o,r-1));
            end
            arriveT(2,ss(o,r))=yx;
            
         end
         
        yx=0;
          for r=4:chenum(o)
            if r==4
                yx=yx++realdist(ss(o,r),retailer+1);
            else
                yx=yx++realdist(ss(o,r),ss(o,r-1));
            end
            arriveT(2,ss(o,r))=yx;
            
         end
        end

        
    elseif chenum(o)==0
   
    else   %%如果车辆的带的工件数量小于3
        
        yx=0;
        for r=1:chenum(o)
            if r==1
                yx=yx++realdist(ss(o,r),11);
            else
                yx=yx++realdist(ss(o,r),ss(o,r-1));
            end
            arriveT(2,ss(o,r))=yx;
            
        end
        arriveT(5,ss(o,chenum(o)))=1;
    end
end

  %%计算该染色体的缺货成本和超时间窗成本
       sumpunish=0;
%        sumkucun=0;
       sumde=0;
       sumpro=0;
    for l=1:retailer
        arriveT(3,l)=arriveT(1,l)+arriveT(2,l);    %%第三行 配送到达时间
        arriveT(4,l)=expectime(1)-arriveT(1,l)-arriveT(2,l);  %%第四行预期达到时间和实际到达时间的插值
  
       cost=arriveT(4,l);
       punish=abs(cost*delaytime(some(j-retailer)));
       sumpunish=sumpunish+punish;

%        kucun=arriveT(2,l)*stockcost(Chrom(i,j-9));
%               sumkucun=sumkucun+kucun;
                             de=arriveT(1,l)*arriveT(5,l);   %%计算的是旅途成本
               sumde=sumde+de;
    end
    
    for mm=1:retailer
        getipro=time(some(retailer))*cousatsa(some(retailer+10),some(retailer));
        sumpro=sumpro+getipro;
    end
    

       t=sumde+sumpunish*0.1+sumpro*0.001;      %得出该染色体对应的加工时间
       F=t+sumrawcous+sumvhicous;