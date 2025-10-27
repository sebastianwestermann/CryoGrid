% c1=999.842594;
% c2=6.793952e-2; 
% c3=9.095290e-3;
% c4=1.001685e-4;
% c5=1.120083e-6;
% c6=6.536332e-9;
% 
% t2 = T.^2;
% t3 = t2.*T;
% t4 = t3.*T;
% t5 = t4.*T;
% 
% term1  =  c1;
% term2  =  c2 .* T;
% term3  = -c3 .* t2;
% term4  =  c4 .* t3;
% term5  = -c5 .* t4;
% term6  =  c6 .* t5;
% 
% density_water  =  term6  + term5  + term4  + term2 + term3 + term1;



T=[0:0.1:30]';
%in g salt/kg water
salt_conc = 0.*T+35;

c1=999.842594;
c2=6.793952e-2; 
c3=9.095290e-3;
c4=1.001685e-4;
c5=1.120083e-6;
c6=6.536332e-9;
d1=8.24493e-1;
d2=4.0899e-3;
d3=7.6438e-5;
d4=8.2467e-7;
d5=5.3875e-9; 
d6=5.72466e-3;
d7=1.0227e-4;
d8=1.6546e-6; 
d9=4.8314e-4;


t2 = T.^2;
t3 = t2.*T;
t4 = t3.*T;
t5 = t4.*T;
s2 = salt_conc.*salt_conc;
s32 = salt_conc.^1.5;

term1  =  c1;
term2  =  c2 .* T;
term3  = -c3 .* t2;
term4  =  c4 .* t3;
term5  = -c5 .* t4;
term6  =  c6 .* t5;
term7  =  d1;
term8  = -d2 .* T;
term9  =  d3 .* t2;
term10  = -d4 .* t3;
term11 =  d5 .* t4;
term12 = -d6;
term13 =  d7 .* T;
term14 = -d8 .* t2;
term15 =  d9;

dpure  =  term6  + term5  + term4  + term2 + term3 + term1;
csal1  = (term11 + term10  + term9  + term8 + term7) .* salt_conc;
csal32 = (term14 + term13 + term12) .* s32;
csal2  =  term15 .* s2;

density_water = dpure + csal1 + csal32 + csal2;