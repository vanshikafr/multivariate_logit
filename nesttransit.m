function [LL]  = nesttransit(beta,...
    I_b, I_m, choice_bus, choice_metro)


beta_metro = beta(1,1);
beta_u= beta(1,2);


bus =  (beta_u.*I_b);
metro = beta_metro + (beta_u.*I_m);

sumtransit = exp(bus)+ exp(metro); 
p_bus= exp(bus)./(sumtransit);
p_metro= exp(metro)./sumtransit;
sum_bus = sum(transpose(choice_bus).*log(p_bus));
sum_metro = sum(transpose(choice_metro).*log(p_metro));

sum_all= [sum_bus, sum_metro];
ll= sum(sum_all);
LL= -ll;

end