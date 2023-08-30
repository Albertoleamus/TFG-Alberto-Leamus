function [Vmp, Imp, Pmp] = MPP_Calculator(IL, I0, Rs, Rsh, nNsVth, Ns)
%Algortimo para el cálculo del MPP con los parámetros hallados, basado en el
%propuesto en el TFG de Alvaro Benito Oliva, sin utilizar la función de Lambert, con los mismos valores iniciales
MaxError = 0.000005; 
VmpA = 0.4*Ns*0.7; 
VmpB = 0.9*Ns*0.7; 
for i=1:1:500
VmpC = (VmpA+VmpB)/2;
I1 = IL+I0;
I2 = 0; 
for j=1:1:200
Ih = (I1+I2)/2;
Err1 = IL-((VmpA+I1*Rs)/Rsh)-I0*(exp((VmpA+I1*Rs)/nNsVth)-1)-I1; 
Err2 = IL-((VmpA+I2*Rs)/Rsh)-I0*(exp((VmpA+I2*Rs)/nNsVth)-1)-I2;
Errh = IL-((VmpA+Ih*Rs)/Rsh)-I0*(exp((VmpA+Ih*Rs)/nNsVth)-1)-Ih; 
if (Err1*Errh < 0)
I2 = Ih;
elseif (Err2*Errh <0) 
I1 = Ih;
end
if(abs(I1-I2) < MaxError)
break;
end
end
ImpA=(I1+I2)/2; 
dP_dV_A=-VmpA*(I0/nNsVth*exp((VmpA+ImpA*Rs)/nNsVth)+1/Rsh)/(1+I0*Rs/nNsVth*exp((VmpA+ImpA*Rs)/nNsVth)+Rs/Rsh)+ImpA;
I1 = IL+I0; 
I2 = 0;  
for j=1:1:200
Ih = (I1+I2)/2; 
Err1 = IL-((VmpB+I1*Rs)/Rsh)-I0*(exp((VmpB+I1*Rs)/nNsVth)-1)-I1; Err2 = IL-((VmpB+I2*Rs)/Rsh)-I0*(exp((VmpB+I2*Rs)/nNsVth)-1)-I2; Errh = IL-((VmpB+Ih*Rs)/Rsh)-I0*(exp((VmpB+Ih*Rs)/nNsVth)-1)-Ih; 
if (Err1*Errh < 0)
I2 = Ih;
elseif (Err2*Errh <0) 
I1 = Ih;
end
if(abs(I1-I2) < MaxError) 
break;
end
end
ImpB = (I1+I2)/2;
dP_dV_B=-VmpB*(I0/nNsVth*exp((VmpB+ImpB*Rs)/nNsVth)+1/Rsh)/(1+I0*Rs/nNsVth*exp((VmpB+ImpB*Rs)/nNsVth)+Rs/Rsh)+ImpB; %dP/dV for VmpB
I1 = IL+I0; 
I2 = 0; 
for j=1:1:200
Ih = (I1+I2)/2; 
Err1 = IL-((VmpC+I1*Rs)/Rsh)-I0*(exp((VmpC+I1*Rs)/nNsVth)-1)-I1; Err2 = IL-((VmpC+I2*Rs)/Rsh)-I0*(exp((VmpC+I2*Rs)/nNsVth)-1)-I2; Errh = IL-((VmpC+Ih*Rs)/Rsh)-I0*(exp((VmpC+Ih*Rs)/nNsVth)-1)-Ih; if (Err1*Errh < 0)
I2 = Ih;
elseif (Err2*Errh <0) I1 = Ih;
end
if(abs(I1-I2) < MaxError) 
break;
end
end
ImpC = (I1+I2)/2;
dP_dV_C=-VmpC*(I0/nNsVth*exp((VmpC+ImpC*Rs)/nNsVth)+1/Rsh)/(1+I0*Rs/nNsVth*exp((VmpC+ImpC*Rs)/nNsVth)+Rs/Rsh)+ImpC; %dP/dV for VmpC
if (dP_dV_C*dP_dV_B < 0) 
VmpA = VmpC;
else
VmpB = VmpC; 
end
if(abs(VmpA-VmpB)<MaxError) 
break
end
end
Vmp = (VmpB+VmpA)/2; 
Imp = (ImpB+ImpA)/2; 
Pmp = Vmp*Imp;
end