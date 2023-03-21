function [LL]  = nestmode(beta,...
    I_transit, Time_auto, choice_transit, choice_auto)

lambda_0= beta(1,1);
beta_auto= beta(1,2);
beta_2= beta(1,3);

transit =  lambda_0.*transpose(I_transit);
auto= beta_auto + (beta_2.*Time_auto);

summode = exp(transit)+ exp(auto); 
p_transit= exp(transit)./(summode);
p_auto= exp(auto)./summode;
sum_transit = sum(transpose(choice_transit).*log(p_transit));
sum_auto = sum(transpose(choice_auto).*log(p_auto));

sum_all= [sum_transit, sum_auto];
ll= sum(sum_all);
LL= -ll;

end

