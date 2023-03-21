n = 1000;
be_1 = -0.5;
be_bm = 1;
be_mm = 0.5;
be_mb = 0.8;
be_metro= 1;
be_u= 1.01;
be_auto= -0.6;
lambda_base= 1;
sig_1=1;
sig_2=1;

T = readtable('base1000data.xlsx');

ID=T.ID;
Time_bb= T.Time_bb;
Time_bm= T.Time_bm;
Time_mb= T.Time_mb;
Time_mm= T.Time_mm;
Time_auto= T.Time_auto;
eta_bb = T.Error_bb;
eta_bm= T.Error_bm;
eta_mb = T.Error_mb;
eta_mm= T.Error_mm;
Error_auto = T.Error_auto;

P= [-1:0.1:1];

s = width(P);
lambda = zeros([1 s]);
corrmat= table(transpose(P), transpose(lambda), transpose(lambda), transpose(lambda));
corrmat.Properties.VariableNames = ["p", 'lambda0', 'lambda1', 'lambda2'];
p=0.2;
for i = 1:s
        p= table2array(corrmat(i, 1)); %takes the first row value for p, generates lambda given that p
        s11= sqrt(sig_1.^2);
        s22= sqrt(1- (p^2))*sig_2;
        s21= p*sig_2;
        S= [s11, 0; s21, s22]; %creates the cholesky matrix
        eta_m= [eta_mb, eta_mm]; %eta columns from train book 
        eta_b= [eta_bb, eta_bm];
        error_m= S*transpose(eta_m);
        error_b= S*transpose(eta_b); %generating joint dist
        Error_bb= error_b(1,:);
        Error_bm= error_b(2,:);
        Error_mb= error_m(1,:);
        Error_mm= error_m(2,:);
        df = table(ID, Time_bb, Time_bm, Time_mb, Time_mm, Time_auto, ...
            transpose(Error_bb), transpose(Error_bm), transpose(Error_mm), transpose(Error_mb), Error_auto);
        df.Properties.VariableNames = ["ID","Time_bb", "Time_bm", "Time_mb" , "Time_mm", "Time_auto", "Error_bb", "Error_bm", "Error_mm", "Error_mb", "Error_auto"];
        
        bb = (be_1.*transpose(Time_bb)) + Error_bb;
        bm = be_bm + (be_1*transpose(Time_bm)) + Error_bm;
        mm = be_mm+ (be_1*transpose(Time_mm)) + Error_mm;
        mb = be_mb + (be_1*transpose(Time_mb)) + Error_mb;
        auto= be_auto + (be_1*transpose(Time_auto))+ transpose(Error_auto);
        df.bb= transpose(bb);
        df.bm= transpose(bm);
        df.mb= transpose(mb);
        df.mm= transpose(mm);
        df.auto= transpose(auto);
   for r = 1:n    
        if df.bb(r,:)> df.bm(r,:)
            df.choice_bb(r,:)= 1;
            df.choice_bm(r,:) = 0;
        else df.choice_bb(r,:)= 0;
             df.choice_bm(r,:) = 1;
        end        
   end   
        
          for r = 1:n    
                if df.mb(r,:)> df.mm(r,:)
                    df.choice_mb(r,:)= 1;
                    df.choice_mm(r,:) = 0;
                else df.choice_mb(r,:)= 0;
                     df.choice_mm(r,:) = 1;
                end        
          end 

        E_bb= (be_1.*transpose(Time_bb));
        E_bm= be_bm+ (be_1.*transpose(Time_bm));
        E_mb= be_mb + (be_1.*transpose(Time_mb));
        E_mm= be_mm + (be_1.*transpose(Time_mm));
        df.E_bb= transpose(E_bb);
        df.E_bm= transpose(E_bm);
        df.E_mb= transpose(E_mb);
        df.E_mm= transpose(E_mm);
        df.E_sumbus = exp(df.E_bb)+ exp(df.E_bm); 
        df.E_summetro = exp(df.E_mb)+ exp(df.E_mm);
        df.Ibus= log(df.E_sumbus);
        df.v_bus= be_u.*df.Ibus; 
        df.Imetro= log(df.E_summetro);
        df.v_metro= be_u.*df.Imetro; 
       
         for r = 1:n    
                if df.v_bus(r,:)> df.v_metro(r,:)
                    df.choice_bus(r,:)= 1;
                    df.choice_metro(r,:) = 0;
                else df.choice_bus(r,:)= 0;
                     df.choice_metro(r,:) = 1;
                end        
         end 
        I_bus= transpose(df.Ibus);
        I_metro= transpose(df.Imetro);
        E_bus= (be_u.*I_bus);
        df.E_bus= transpose(E_bus);
        E_metro= be_metro + (be_1.*I_metro);
        df.E_metro = transpose(E_metro);
        df.E_sumtransit = exp(df.E_bus)+ exp(df.E_metro); 
        df.Ik= log(df.E_sumtransit);
        df.v_transit= lambda_base.*df.Ik;
        
         for r = 1:n    
                if df.v_transit(r,:)> df.auto(r,:)
                    df.choice_transit(r,:)= 1;
                    df.choice_auto(r,:) = 0;
                else df.choice_transit(r,:)= 0;
                     df.choice_auto(r,:) = 1;
                end        
         end 
        
        choice_bb= df.choice_bb;
        choice_bm= df.choice_bm;
        choice_mm= df.choice_mm;
        choice_mb= df.choice_mb;
        choice_bus= df.choice_bus;
        choice_metro= df.choice_metro;
        choice_transit= df.choice_transit;
        choice_auto= df.choice_auto;
        b_1 = -0.1;
        b_bm= 1;
        b_mm= 1;
        b_mb=1;
        beta0   = [b_1, b_bm, b_mm, b_mb];
        
        options = optimset('Display','iter','Maxiter',1000,'MaxFunEvals',10000,'TolFun',1e-50);
        [betahat3,value,exit,output] = fminunc('nestbm',beta0,options,...
               Time_bb, Time_bm, Time_mb , Time_mm, Time_auto, choice_bb, choice_bm, choice_mm, choice_mb, choice_bus, choice_metro, choice_transit, choice_auto);
        
        beta_1 = betahat3(1,1);
        beta_bm= betahat3(1,2);
        beta_mm= betahat3(1,3);
        beta_mb= betahat3(1,4);
        lambda_2= array2table(betahat3(1,1));
        corrmat(i,4)= lambda_2;

        
        bb = (beta_1.*Time_bb);
        bm = beta_bm + (beta_1.*Time_bm);
        mm = beta_mm+ (beta_1.*Time_mm);
        mb = beta_mb + (beta_1.*Time_mb);
        sumb = exp(bb)+ exp(bm); 
        summ= exp(mm) + exp(mb);
        I_b= log(sumb);
        I_m= log(summ);
        
        b_metro = 1.5;
        b_u= 0.5;
        beta0transit= [b_metro, b_u];
        
        options = optimset('Display','iter','Maxiter',1000,'MaxFunEvals',10000,'TolFun',1e-50);
        [betahat2,value,exit,output] = fminunc('nesttransit',beta0transit,options,...
               I_b, I_m,  choice_bus, choice_metro);
        
        beta_metro = betahat2(1,1);
        beta_u= betahat2(1,2);
        lambda_1= array2table(betahat2(1,2));
        corrmat(i,3)= lambda_1;
        
        bus =  (beta_u.*I_b);
        metro = beta_metro + (beta_u.*I_m);
        sumtransit= exp(bus)+ exp(metro); 
        I_transit= log(sumtransit);
        
        b_auto = 0.5;
        lambda_0= 0.5;
        beta_2= -0.1;
        beta0mode= [lambda_0, b_auto, beta_2];
        
        options = optimset('Display','iter','Maxiter',1000,'MaxFunEvals',10000,'TolFun',1e-50);
        [betahat1, value,exit,output] = fminunc('nestmode',beta0mode,options,...
              I_transit, Time_auto, choice_transit, choice_auto);

        lambda_0= array2table(betahat1(1,1));
       corrmat(i,2)= lambda_0;

end

writetable(corrmat,'n1000nests3.xlsx');
