%This program was written by Chang Gyo Jung, which based on Xuhui Zhou's
%version to study SOM decomposition using MCMC

%clear all
for rep=1:3
for nput=2:5;%nput=2; % 2-control, 3-warming, 4-clipping, 5-warming+clipping
for year=1:16;
format long e;
RandSeed=clock;
rng('default');
rng(RandSeed(6));

resp_type = 1;%%%root-exclusion-1; intact-soil-2

%start
if resp_type==1
    disp('root-exclusion')
    fn_parameterr = 'paramrter_hetero';
    soil_respration =   soil_resp_hetero(160);
    %soil carbon
    if year==1
        if nput == 2||nput ==3
            soilCarbonn = soil_UnClip_(2);
            soilCarbon = soilCarbonn(1,:);
            soilCTime        =   soilCarbon(:,1);
            soilCValue       =   soilCarbon(:,nput);
        end
    elseif year==9
        if nput == 2||nput == 3
            soilCarbonn = soil_UnClip_(2);
            soilCarbon = soilCarbonn(2,:);
            soilCTime        =   soilCarbon(:,1);
            soilCValue       =   soilCarbon(:,nput);
        end
    end
else
    disp('intact-soil') 
     fn_parameterr = 'paramrter_total';
    soil_respration =   soil_resp_total(189);
    %soil carbon
    if year==1
        if nput == 2||nput ==3
            soilCarbonn = soil_UnClip_total(2);
            soilCarbon = soilCarbonn(1:2,:);
            soilCTime        =   soilCarbon(:,1);
            soilCValue       =   soilCarbon(:,nput);
        end
    end
    if nput == 2||nput ==3
        if year==2||year==3||year==4||year==6||year==7||year==11||year==14
            soilCarbonn = soil_UnClip_total(18);
            soilCarbon = soilCarbonn(365*year==soilCarbonn(:,1),:);
            soilCTime        =   soilCarbon(:,1);
            soilCValue       =   soilCarbon(:,nput);
        end
    else
        if year==11||year==14
            soilCarbonn = soil_Clip_total(2);
            soilCarbon = soilCarbonn(365*year==soilCarbonn(:,1),:);
            soilCTime        =   soilCarbon(:,1);
            soilCValue       =   soilCarbon(:,nput);
        end
    end
    
    end


nsim     =  50000;         % control running time
Nt       =   365;   %five compartments
HS       =   10;
m        =   5;   
mscut    =   0.2; 
cbnScale =   1.0; 

fn_parameter = strcat("year",num2str(year),"_", fn_parameterr,"_", num2str(rep),"_",date,".xlsx");

%*********************prior c--the transfer coefficients********************
if resp_type==1 %hetero
cmin(1) = 0.5e-5; cmin(2) = 0.1e-4; cmin(3) = 0.1e-6;  cmin(4) = 1e-9;cmin(12) = 0.6;cmin(13) = 0.05;
cmax(1) = 6.5e-3; cmax(2) = 5.4e-3; cmax(3) = 4.5e-4;  cmax(4) = 1.0e-5;cmax(12) = 5;cmax(13) = 1;
     
cmin(5) = 1e-1;cmin(6) = 5e-2;   cmin(7) = 2e-1; cmin(8) = 0.1e-3; cmin(9) = 1e-1;cmin(10) = 1e-2;cmin(11) = 3e-1;
cmax(5) = 7e-1;cmax(6) = 1.5e-1; cmax(7) = 7e-1; cmax(8) = 8.4e-3; cmax(9) = 6e-1;cmax(10) = 2e-1;cmax(11) = 7e-1;

else %total
cmin(1) = 1e-6;cmin(2) = 1e-7; cmin(3) = 0.5e-5; cmin(4) = 0.1e-4; cmin(5) = 0.1e-6;  cmin(6) = 1e-9;cmin(16) = 0.6;cmin(17) = 0.05;
cmax(1) = 1.5e-2;cmax(2) = 8.7e-4; cmax(3) = 6.5e-3; cmax(4) = 5.4e-3; cmax(5) = 4.5e-4;  cmax(6) = 1.0e-5;cmax(16) = 5;cmax(17) = 1;

cmin(7) = 1e-1;cmin(8) = 5e-2;   cmin(9) = 2e-1; cmin(10) = 0.1e-3; cmin(11) = 1e-1;  cmin(12) = 1e-2;cmin(13) = 3e-1;cmin(14) = 0.1;cmin(15) = 0.1;
cmax(7) = 7e-1;cmax(8) = 1.5e-1; cmax(9) = 7e-1; cmax(10) = 8.4e-3; cmax(11) = 6e-1;  cmax(12) = 2e-1;cmax(13) = 7e-1;cmax(14) = 0.5;cmax(15) = 0.5;

end

%*********************set input matrix and input carbon******************
yearr = 365 * year;

%select 365 days
if year==1
    yearr_f = 1;
else yearr_f = 365 * (year -1) + 1;
end
carbon_input = carbon2;
u=carbon_input(yearr_f:yearr,nput);

%data sets below; soil respiration, biomass, and soil carbon
bioMasss        =   biom(16);
bioMass         =   bioMasss(year,:);
biomTime        =   bioMass(:,1);
biomValue       =   bioMass(:,nput);  

%bnpp
rootbioMasss        =   bnpp(16);
rootbioMass         =   rootbioMasss(year,:);
rootbiomTime        =   rootbioMass(:,1);
rootbiomValue       =   rootbioMass(:,nput);  

%soil respiration
soil_resprationn = soil_respration(yearr_f<soil_respration(:,1),:);
soil_resprationnn = soil_resprationn(yearr>soil_resprationn(:,1),:);
soilTime        =   soil_resprationnn(:,1);
soilResValue    =   soil_resprationnn(:,nput);

cdif  = (cmax-cmin)';

%simulation starts
c_op=cmin'+rand*cdif;
J_last = 300000;

record_index=1;
DJ1=2*var(soilResValue);
if resp_type==2
    if year ==1
        if nput==2||nput==3
            DJ4=2*var(soilCValue);
        end
    end
end

   J(1)  =   0;
   J(2)  =   0;
   J(3)  =   0;
   J(4)  =   0;
   J(5)  =   0;
%Simulation starts
upgraded=0;
for simu=1:nsim
    counter=simu
    upgradedd=upgraded
    d = 7; %default 7
     if resp_type==2 %root-exclusion-1; intact soil-2; 
         c_new=Generate_tot2(c_op,cmin,cmax,d);  %generate pars
     else
         c_new=Generate_het2(c_op,cmin,cmax,d);  %generate pars
     end
    
if resp_type==1 %hetero
    if nput==2
        x00 = x0_UC_tot;
        x000 = x0_UC_het;
        tot_mle = mlepars_UC_tot;
    elseif nput==3
        x00 = x0_UW_tot;
        x000 = x0_UW_het;
        tot_mle = mlepars_UW_tot;
    elseif nput==4
        x00 = x0_CC_tot;
        x000 = x0_CC_het;
        tot_mle = mlepars_CC_tot;
    elseif nput==5
        x00 = x0_CW_tot;
        x000 = x0_CW_het;
        tot_mle = mlepars_CW_tot;
    end
    
    tau=temp_moist(Nt,nput,yearr_f,yearr,mscut,c_new(13),c_new(12)); %het
    tau_tot=temp_moist(Nt,nput,yearr_f,yearr,mscut,tot_mle(17,year),tot_mle(16,year));%tot

    phi_slResp = [0 0 (1-c_new(5)-c_new(6))*c_new(1) (1-c_new(7)-c_new(8))*c_new(2) (1-c_new(9)-c_new(10))*c_new(3) (1-c_new(11))*c_new(4)];
    phi_foilage  = [1 0 0 0 0 0];
    phi_root  = [0 1 0 0 0 0];    
    phi_soilC1 = [0 0 0 1 1 1]; 

    A=[-1      0       0          0          0           0
        0     -1       0          0          0           0
        1      0       -1         0          0           0
        0      0   c_new(5)      -1       c_new(9)   c_new(11)
        0      0   c_new(6)    c_new(7)     -1           0
        0      0       0       c_new(8)  c_new(10)     -1];

    K=diag([tot_mle(1,year) tot_mle(2,year) c_new(1) c_new(2) c_new(3) c_new(4)]);
    x0  = [x00(1,year) x00(2,year) x000(:,year)']'; %initial value of states]';
    b   =  [tot_mle(14,year) tot_mle(15,year) 0 0 0 0]'; 
    x_last = x0;
    for i = 1:Nt
        x_present = (eye(6)+A*K*tau(i))*x_last + b*u(i);
        x(:,i) = x_present;
        x_last = x_present;
    end
    
else %total
    phi_slResp = [0 0 (1-c_new(7)-c_new(8))*c_new(3) (1-c_new(9)-c_new(10))*c_new(4) (1-c_new(11)-c_new(12))*c_new(5) (1-c_new(13))*c_new(6)];
    phi_foilage  = [1 0 0 0 0 0];
    phi_root  = [0 1 0 0 0 0];    
    phi_soilC1 = [0 0 0 1 1 1]; 
    
    A=[-1      0       0          0          0           0
        0     -1       0          0          0           0
        1      1       -1         0          0           0
        0      0   c_new(7)      -1       c_new(11)   c_new(13)
        0      0   c_new(8)    c_new(9)     -1           0
        0      0       0       c_new(10)  c_new(12)     -1];

    K=diag(c_new(1:6));
    tau=temp_moist(Nt,nput,yearr_f,yearr,mscut,c_new(17),c_new(16));

if nput==2
    x00 = x0_UC_tot;
elseif nput==3
    x00 = x0_UW_tot;
elseif nput==4
    x00 = x0_CC_tot;
elseif nput==5
    x00 = x0_CW_tot;
end
    
    x0 = x00(:,year);    
    x_last=x0;
    b   =  [c_new(14) c_new(15) 0 0 0 0]';  

 for i = 1:365 
    x_present = (eye(6)+A*K*tau(i))*x_last + b*u(i);
    x(:,i) = x_present;
    x_last = x_present;
 end
end
   
if resp_type==1    
     if year==1
        for i=1:length(soilTime)       
            soilResp_simu(i)    = tau(soilTime(i))*phi_slResp*x(:,soilTime(i));
        end
    else
        for i=1:length(soilTime)      
            soilResp_simu(i)    = tau(soilTime(i)-365*(year-1))*phi_slResp*x(:,(soilTime(i)-365*(year-1)));
        end
    end

    if nput == 2||nput == 3
        if year==1
            for i=1:length(soilCTime)
                soilC1_simu(i)  = phi_soilC1*x(:,soilCTime(i));
            end
            
        elseif year==9
            for i=1:length(soilCTime)
                soilC1_simu  = phi_soilC1*x(:,soilCTime(i)-365*(year-1));
            end
        end
                
    end 
else
    if year==1
        for i=1:length(soilTime)       
            soilResp_simu(i)    = tau(soilTime(i))*phi_slResp*x(:,soilTime(i))+0.25*(1-b(1)-b(2))*u(soilTime(i));%
        end
        Biom_simu    = phi_foilage*x(:,biomTime(year));
    else
        for i=1:length(soilTime)      
            soilResp_simu(i)    = tau([soilTime(i)-365*(year-1)])*phi_slResp*x(:,(soilTime(i)-365*(year-1)))+0.25*(1-b(1)-b(2))*u(soilTime(i)-365*(year-1));%
        end 
        Biom_simu    = phi_foilage*x(:,biomTime-365*(year-1));
        rootBiom_simu    = phi_root*x(:,rootbiomTime-365*(year-1));
        
    end
    
    if nput==2||nput==3
        if year==1
            for i=1:length(soilCTime)
                soilC1_simu(i)  = phi_soilC1*x(:,soilCTime(i));
            end
        elseif year==2||year==3||year==4||year==6||year==7||year==11||year==14
                soilC1_simu  = phi_soilC1*x(:,365);
        end
     elseif nput==4||nput==5 
        if year==11||year==14
           soilC1_simu  = phi_soilC1*x(:,365);
        end
    end
end
                           
%cost function J
   if resp_type==2
       j1 = 5;
       j2 = 0.001;
       j3 = 0.01;
       j5 = 100;
       
       J(1)  = (norm(soilResp_simu-soilResValue'))^2*j1;
       J(2)  = (norm(Biom_simu-biomValue))^2*j2;
       J(5)  = (norm(mean(soilResp_simu)-mean(soilResValue)))^2*j5;
       
       if year==4||year==5||year==6||year==7||year==8||year==9||year==11||year==12||year==13||year==14||year==15||year==16;
           J(3)  =   (norm(rootBiom_simu-rootbiomValue))^2*j3;  
           
       end
   else
       j1 = 100;
       j5 = 800;
  
       J(1)  = (norm(soilResp_simu-soilResValue'))^2*j1;
       J(5)  = (norm(mean(soilResp_simu)-mean(soilResValue)))^2*j5;
   end
      
   %soil carbon
   if resp_type==1
       if year==1||year==9
           if nput==2||nput==3
             J(4)  =   (norm(soilC1_simu-soilCValue'))^2*0.00001;
           end
       end
   else
       if year==1||year==2||year==3||year==4||year==6||year==7
           if nput==2||nput==3
               J(4)  =   (norm(soilC1_simu-soilCValue'))^2*0.00001;
           end
       elseif year==11||year==14
               J(4)  =   (norm(soilC1_simu-soilCValue'))^2*0.00001;
       end
   end
 
   if resp_type==1
     if year==1||year==9
        if nput==2||nput==3
            J_new=(J(1)/DJ1+J(4)+J(5));
        else
            J_new=(J(1)/DJ1+J(5));
        end
       else
        J_new=(J(1)/DJ1+J(5));
       end
        
   else
       if year==1
        if nput==2||nput==3
            J_new=(J(1)/DJ1+J(2)+J(4)/DJ4+J(5));
        elseif nput==4||nput==5
            J_new=(J(1)/DJ1+J(2)+J(5));
        end             
       else
           J_new=(J(1)/DJ1+J(2)+J(3)+J(4)+J(5));
       end
   end
           
   delta_J = J_new-J_last;
	
	if min(1, exp(-delta_J)) >rand
        c_op=c_new;
        J_last=J_new;
        upgraded=upgraded+1
        c_upgraded(:,upgraded)=c_op;
        J_upgraded(:,upgraded)=J_last; 
    end
    
    c_rec(:,record_index)=c_op;
    J_record(record_index)=J_last;
    record_index=record_index+1;
end

%save posterior distribution
if resp_type==1
    c1=c_upgraded(1,end/2:end);
    c2=c_upgraded(2,end/2:end);
    c3=c_upgraded(3,end/2:end);
    c4=c_upgraded(4,end/2:end);
    c5=c_upgraded(5,end/2:end);
    c6=c_upgraded(6,end/2:end);
    c7=c_upgraded(7,end/2:end);
    c8=c_upgraded(8,end/2:end);
    c9=c_upgraded(9,end/2:end);
    c10=c_upgraded(10,end/2:end);
    c11=c_upgraded(11,end/2:end);
    c12=c_upgraded(12,end/2:end);
    c13=c_upgraded(13,end/2:end);

if nput == 2
xlswrite(fn_parameter, c1','UC','A'); 
xlswrite(fn_parameter, c2','UC','B');
xlswrite(fn_parameter, c3','UC','C');
xlswrite(fn_parameter, c4','UC','D');
xlswrite(fn_parameter, c5','UC','E');
xlswrite(fn_parameter, c6','UC','F');
xlswrite(fn_parameter, c7','UC','G');
xlswrite(fn_parameter, c8','UC','H'); 
xlswrite(fn_parameter, c9','UC','I');
xlswrite(fn_parameter, c10','UC','J');
xlswrite(fn_parameter, c11','UC','K');
xlswrite(fn_parameter, c12','UC','L');
xlswrite(fn_parameter, c13','UC','M');
elseif nput == 3
        xlswrite(fn_parameter, c1','UW','A'); 
        xlswrite(fn_parameter, c2','UW','B');
        xlswrite(fn_parameter, c3','UW','C');
        xlswrite(fn_parameter, c4','UW','D');
        xlswrite(fn_parameter, c5','UW','E');
        xlswrite(fn_parameter, c6','UW','F');
        xlswrite(fn_parameter, c7','UW','G');
        xlswrite(fn_parameter, c8','UW','H'); 
        xlswrite(fn_parameter, c9','UW','I');
        xlswrite(fn_parameter, c10','UW','J');
        xlswrite(fn_parameter, c11','UW','K');
        xlswrite(fn_parameter, c12','UW','L');
        xlswrite(fn_parameter, c13','UW','M');
elseif nput == 4
            xlswrite(fn_parameter, c1','CC','A'); 
            xlswrite(fn_parameter, c2','CC','B');
            xlswrite(fn_parameter, c3','CC','C');
            xlswrite(fn_parameter, c4','CC','D');
            xlswrite(fn_parameter, c5','CC','E');
            xlswrite(fn_parameter, c6','CC','F');
            xlswrite(fn_parameter, c7','CC','G');
            xlswrite(fn_parameter, c8','CC','H'); 
            xlswrite(fn_parameter, c9','CC','I');
            xlswrite(fn_parameter, c10','CC','J');
            xlswrite(fn_parameter, c11','CC','K');
            xlswrite(fn_parameter, c12','CC','L');
            xlswrite(fn_parameter, c13','CC','M');
elseif nput ==5
             xlswrite(fn_parameter, c1','CW','A'); 
             xlswrite(fn_parameter, c2','CW','B');
             xlswrite(fn_parameter, c3','CW','C');
             xlswrite(fn_parameter, c4','CW','D');
             xlswrite(fn_parameter, c5','CW','E');
             xlswrite(fn_parameter, c6','CW','F');
             xlswrite(fn_parameter, c7','CW','G');
             xlswrite(fn_parameter, c8','CW','H'); 
             xlswrite(fn_parameter, c9','CW','I');
             xlswrite(fn_parameter, c10','CW','J');
             xlswrite(fn_parameter, c11','CW','K');
             xlswrite(fn_parameter, c12','CW','L');
             xlswrite(fn_parameter, c13','CW','M');
end
    
else
   c1=c_upgraded(1,end/2:end);
   c2=c_upgraded(2,end/2:end);
   c3=c_upgraded(3,end/2:end);
   c4=c_upgraded(4,end/2:end);
   c5=c_upgraded(5,end/2:end);
   c6=c_upgraded(6,end/2:end);
   c7=c_upgraded(7,end/2:end);
   c8=c_upgraded(8,end/2:end);
   c9=c_upgraded(9,end/2:end);
   c10=c_upgraded(10,end/2:end);
   c11=c_upgraded(11,end/2:end);
   c12=c_upgraded(12,end/2:end);
   c13=c_upgraded(13,end/2:end);
   c14=c_upgraded(14,end/2:end);
   c15=c_upgraded(15,end/2:end);
   c16=c_upgraded(16,end/2:end);
   c17=c_upgraded(17,end/2:end);
    
if nput == 2
xlswrite(fn_parameter, c1','UC','A'); 
xlswrite(fn_parameter, c2','UC','B');
xlswrite(fn_parameter, c3','UC','C');
xlswrite(fn_parameter, c4','UC','D');
xlswrite(fn_parameter, c5','UC','E');
xlswrite(fn_parameter, c6','UC','F');
xlswrite(fn_parameter, c7','UC','G');
xlswrite(fn_parameter, c8','UC','H'); 
xlswrite(fn_parameter, c9','UC','I');
xlswrite(fn_parameter, c10','UC','J');
xlswrite(fn_parameter, c11','UC','K');
xlswrite(fn_parameter, c12','UC','L');
xlswrite(fn_parameter, c13','UC','M');
xlswrite(fn_parameter, c14','UC','N');
xlswrite(fn_parameter, c15','UC','O');
xlswrite(fn_parameter, c16','UC','P');
xlswrite(fn_parameter, c17','UC','Q');

elseif nput == 3
        xlswrite(fn_parameter, c1','UW','A'); 
        xlswrite(fn_parameter, c2','UW','B');
        xlswrite(fn_parameter, c3','UW','C');
        xlswrite(fn_parameter, c4','UW','D');
        xlswrite(fn_parameter, c5','UW','E');
        xlswrite(fn_parameter, c6','UW','F');
        xlswrite(fn_parameter, c7','UW','G');
        xlswrite(fn_parameter, c8','UW','H'); 
        xlswrite(fn_parameter, c9','UW','I');
        xlswrite(fn_parameter, c10','UW','J');
        xlswrite(fn_parameter, c11','UW','K');
        xlswrite(fn_parameter, c12','UW','L');
        xlswrite(fn_parameter, c13','UW','M');
        xlswrite(fn_parameter, c14','UW','N');
        xlswrite(fn_parameter, c15','UW','O');
        xlswrite(fn_parameter, c16','UW','P');
        xlswrite(fn_parameter, c17','UW','Q');

elseif nput == 4
            xlswrite(fn_parameter, c1','CC','A'); 
            xlswrite(fn_parameter, c2','CC','B');
            xlswrite(fn_parameter, c3','CC','C');
            xlswrite(fn_parameter, c4','CC','D');
            xlswrite(fn_parameter, c5','CC','E');
            xlswrite(fn_parameter, c6','CC','F');
            xlswrite(fn_parameter, c7','CC','G');
            xlswrite(fn_parameter, c8','CC','H'); 
            xlswrite(fn_parameter, c9','CC','I');
            xlswrite(fn_parameter, c10','CC','J');
            xlswrite(fn_parameter, c11','CC','K');
            xlswrite(fn_parameter, c12','CC','L');
            xlswrite(fn_parameter, c13','CC','M');
            xlswrite(fn_parameter, c14','CC','N');
            xlswrite(fn_parameter, c15','CC','O');
            xlswrite(fn_parameter, c16','CC','P');
            xlswrite(fn_parameter, c17','CC','Q');

elseif nput ==5
             xlswrite(fn_parameter, c1','CW','A'); 
             xlswrite(fn_parameter, c2','CW','B');
             xlswrite(fn_parameter, c3','CW','C');
             xlswrite(fn_parameter, c4','CW','D');
             xlswrite(fn_parameter, c5','CW','E');
             xlswrite(fn_parameter, c6','CW','F');
             xlswrite(fn_parameter, c7','CW','G');
             xlswrite(fn_parameter, c8','CW','H'); 
             xlswrite(fn_parameter, c9','CW','I');
             xlswrite(fn_parameter, c10','CW','J');
             xlswrite(fn_parameter, c11','CW','K');
             xlswrite(fn_parameter, c12','CW','L');
             xlswrite(fn_parameter, c13','CW','M');
             xlswrite(fn_parameter, c14','CW','N');
             xlswrite(fn_parameter, c15','CW','O');
             xlswrite(fn_parameter, c16','CW','P');
             xlswrite(fn_parameter, c17','CW','Q');
end
end
  clearvars -except year rep nput;        
end
end
end