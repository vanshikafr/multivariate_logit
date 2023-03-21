n = 1000;
be_metro= 1;
be_auto= 1;
be_1= -0.5;
sig_1=1;
sig_2=1;
lambda_base= 1;
T = readtable('base1000data.xlsx');

ID=T.ID;
Time_bus= T.Time_bus;
Time_metro= T.Time_metro;
Time_auto= T.Time_auto;
eta_bus = T.Error_bus;
eta_metro = T.Error_metro;
Error_auto = T.Error_auto;

P= [-1:0.1:1];

s = width(P);
lambda = zeros([1 s]);
corrmat= table(transpose(P), transpose(lambda), transpose(lambda));
corrmat.Properties.VariableNames = ["p", 'lambda0', 'lambda1'];

for i = 1:s
        p= table2array(corrmat(i, 1)); %takes the first row value for p, generates lambda given that p
        s11= sqrt(sig_1.^2);
        s22= sqrt(1- (p^2))*sig_2;
        s21= p*sig_2;
        S= [s11, 0; s21, s22]; %creates the cholesky matrix
        eta_transit= [eta_bus, eta_metro]; %eta columns from train book
        error_transit= S*transpose(eta_transit); %generating joint dist
        Error_bus= error_transit(1,:);
        Error_metro= error_transit(2,:);
        df = table(ID, Time_bus, Time_metro, Time_auto, ...
            transpose(Error_bus), transpose(Error_metro), Error_auto);
        df.Properties.VariableNames = ["ID","Time_bus", "Time_metro", "Time_auto", "Error_bus", "Error_metro", "Error_auto"];
        
        bus = (be_1.*transpose(Time_bus)) + Error_bus;
        metro = be_metro + (be_1*transpose(Time_metro)) + Error_metro;
        auto= be_auto + (be_1*transpose(Time_auto))+ transpose(Error_auto);
        df.bus= transpose(bus);
        df.metro= transpose(metro);
        df.auto= transpose(auto);

         for r = 1:n    
                if df.bus(r,:)> df.metro(r,:)
                    df.choice_bus(r,:)= 1;
                    df.choice_metro(r,:) = 0;
                else df.choice_bus(r,:)= 0;
                     df.choice_metro(r,:) = 1;
                end        
         end 
      
        E_bus= (be_1.*transpose(Time_bus));
        df.E_bus= transpose(E_bus);
        E_metro= be_metro + (be_1.*transpose(Time_metro));
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
        
        choice_bus= df.choice_bus;
        choice_metro= df.choice_metro;
        choice_transit= df.choice_transit;
        choice_auto= df.choice_auto;
        b_1 = 1;
        b_metro= 2;
        beta1transit = [b_1, b_metro];
        
        options = optimset('Display','iter','Maxiter',1000,'MaxFunEvals',10000,'TolFun',1e-50);
        [betahat1,value,exit,output] = fminunc('nest2transit',beta1transit,options,...
               Time_bus, Time_metro, Time_auto, choice_bus, choice_metro, choice_transit, choice_auto);
        
        beta_1 = betahat1(1,1);
        beta_metro= betahat1(1,2);
        lambda_1= array2table(betahat1(1,1));
        corrmat(i,3)= lambda_1;

        
        bus = (beta_1.*transpose(Time_bus));
        metro = beta_metro + (beta_1.*transpose(Time_metro));
        sumtransit = exp(bus)+ exp(metro); 
        I_transit= log(sumtransit);
        
        b_0= 0.5;
        b_auto = 0.5;
        beta_2= -0.1;
        beta0mode= [b_0, b_auto, beta_2];
        
        options = optimset('Display','iter','Maxiter',1000,'MaxFunEvals',10000,'TolFun',1e-50);
        [betahat1, value,exit,output] = fminunc('nest2mode',beta0mode,options,...
              I_transit, Time_auto, choice_transit, choice_auto);

        lambda_0= array2table(betahat1(1,1));
       corrmat(i,2)= lambda_0;

end

writetable(corrmat,'n1000nests2.xlsx');
