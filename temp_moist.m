function y=temp_moist(Nt,nput,yearr_f,yearr,mscut,R10,Q10);
moist_1l    =   mois5_tao2(5840);
moist_1    =   moist_1l(yearr_f:yearr,1:5);
temp_1l     =   tmp_tao(5840);
temp_1     =   temp_1l(yearr_f:yearr,1:5);
temp       =   temp_1(:,nput);
moist      =   moist_1(:,nput);
for i = 1:(Nt)
   sumtemp = 0;
   tmp=temp(i);
   if (i > 10)
      for j = i-9:i
         sumtemp=sumtemp+temp(j);
      end
      tmp=sumtemp/10;
   end
   tmp=R10*Q10^((tmp-10)/10); % Yuanhe Yang revised
   moisture=1;
   if (moist(i)<mscut)
       moisture=1.0-1/mscut*(mscut-moist(i)); % Yuanhe Yang revised
   end
   product(i)=tmp*moisture;
%     if (product(i)>=4)
%         product(i)=4;
%     end       
end
y=product;
  