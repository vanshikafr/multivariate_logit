function [LL]  = nestbm(beta,...
    Time_bb, Time_bm, Time_mb , Time_mm, Time_auto, choice_bb, choice_bm, choice_mm, choice_mb, choice_bus, choice_metro, choice_transit, choice_auto)

beta_1 = beta(1,1);
beta_bm= beta(1,2);
beta_mm= beta(1,3);
beta_mb= beta(1,4);



bb = (beta_1.*transpose(Time_bb));
bm = beta_bm + (beta_1.*transpose(Time_bm));
mm = beta_mm+ (beta_1.*transpose(Time_mm));
mb = beta_mb + (beta_1.*transpose(Time_mb));


%nest bus
sumb = exp(bb)+ exp(bm); 
p_bb= exp(bb)./(sumb);
p_bm= exp(bm)./sumb;
sum_bb = sum(transpose(choice_bb).*log(p_bb));
sum_bm = sum(transpose(choice_bm).*log(p_bm));


%nest metro
summ= exp(mm) + exp(mb);
p_mm = exp(mm)./summ;
p_mb= exp(mb)./summ;
sum_mm = sum(transpose(choice_mm).*log(p_mm));
sum_mb = sum(transpose(choice_mb).*log(p_mb));

sum_all= [sum_bb, sum_bm, sum_mm, sum_mb];
ll= sum(sum_all);
LL= -ll;

end