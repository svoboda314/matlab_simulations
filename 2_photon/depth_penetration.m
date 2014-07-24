NA = 1;
no = 1.3;

lambda_800=0.800;
lambda_1280=1.280;
lambda_1700=1.700;

le_800 = 140;
le_1280 = 285;
le_1700 = 400; 

z=400:1000;

sbr2p_800=(6*NA^2*z.^2/(lambda_800*le_800)).*exp(-2*z/le_800);
semilogy(z, sbr2p_800, 'k')

hold

z=400:2000;

sbr2p_1280=(6*NA^2*z.^2/(lambda_1280*le_1280)).*exp(-2*z/le_1280);
plot(z, sbr2p_1280, 'b')

z=1000:3000;

sbr3p_1280=(14.7*z.^4*NA^6/(lambda_1280^3*le_1280)).*exp(-3*z/le_1280)
plot(z, sbr3p_1280, 'g')

z=1000:4000;

sbr3p_1700=(14.7*z.^4*NA^6/(lambda_1700^3*le_1700)).*exp(-3*z/le_1700)
plot(z, sbr3p_1700, 'r')

plot(400:4000, 1, 'k')

hold off