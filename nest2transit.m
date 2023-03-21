function [LL]  = nest2transit(beta,...
    Time_bus, Time_metro, Time_auto, choice_bus, choice_metro, choice_auto, choice_transit)

lambda_1 = beta(1,1);
beta_metro= beta(1,2);

bus = (lambda_1.*transpose(Time_bus));
metro = beta_metro + (lambda_1.*transpose(Time_metro));

%nest transit
sumtransit = exp(bus)+ exp(metro); 
p_bus= exp(bus)./(sumtransit);
p_metro= exp(metro)./sumtransit;
logbus= log(p_bus);
logmetro= log(p_metro);

sum_bus = sum(transpose(choice_bus).*logbus);
sum_metro = sum(transpose(choice_metro).*logmetro);

sum_all= [sum_bus, sum_metro];
ll= sum(sum_all);
LL= -ll;

end