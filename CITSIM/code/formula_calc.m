k = 1.3806505e-23 ;      % Boltzmann constant (J/K)
R = 8.3145 ;     % Ideal gas constant (J/(mol*K))
kg_amu = 1.6605402e-27; %(kg/amu) conversion factor
eV_J = 6.2415095e+18;   % (eV/J) conversion factor
c= 29979245800  ;        % cm.s^-1
h= 6.6260695729e-34  ;  % J.s^-1

%Mean free path and Knudsen number
T=300;
r=1.8e-10; % nitrogen VDW radius
P_mtorr=1;
P_Pa=P_mtorr/7.50061683
P_mbar=P_Pa/100;
lambda=k*T/(sqrt(2)*pi*(2*r)^2*P_Pa)
L=1e-4;
Kn=lambda/L

%% p_ratio
gamma=1+2/5
p_ratio=(2/(gamma+1))^(gamma/(gamma-1))

%% Reynolds number and viscosity
%http://www.mhtl.uwaterloo.ca/old/onlinetools/airprop/airprop.html
gamma=1.4;
nu_wiki = @(m,d,P_Pa) (1/3)*sqrt(8*k*T/(pi*m*kg_amu))*k*T/(sqrt(2)*pi*(2*r)^2*P_Pa);
nu_HS= @(m,d) (5/16)*(1/((2*r)^2))*sqrt(k*m*kg_amu*T/pi)*(k*T/P_Pa);
cs = @(m,gamma) sqrt(gamma*k*T/(m*kg_amu));
Mflow = @(gamma,p1,p2,T1,m) (p1/(R/(m*1e-3))/T1)*sqrt(2*gamma/(gamma-1)*(R/(m*1e-3))*T1*(p2/p1)^(2/gamma)*(1-(p2/p1)^((gamma-1)/gamma)));
massflow = Mflow(1.66,133.3,120,300,4)
%%
%speed of sound
cs_N2=cs(28,1.4); %m.s^-1, gamma heat cap. ratio N2 at 20C (not much T dep.)
cs_He=cs(4,1.66);%m.s^-1, gamma at 20C

%kinematic viscosities
length_trap=1e-3;% meter
%N2 VDW radius 2.25e-10m, other value tabulated 2.16A
nu_wiki_N2=nu_wiki(28,2.25e-10,P_Pa);%m^2.s-1 
nu_HS_N2=nu_HS(28,2.25e-10);
%He VDW radius 2.1e-10m, other value tabulated 1.4A
nu_wiki_He=nu_wiki(4,2.1e-10,P_Pa);%m^2.s-1 
nu_HS_He=nu_HS(4,2.1e-10);

%Reynolds number
Re_wiki_N2=cs_N2*length_trap/nu_wiki_N2
Re_HS_N2=cs_N2*length_trap/nu_HS_N2
Re_wiki_He=cs_He*length_trap/nu_wiki_He;
Re_HS_He=cs_He*length_trap/nu_HS_He;

%lattice spacings
dt_lbm=0.0005;
dx_lbm=1/192;

%lattice viscosity and tau
nu_lbm_wiki_N2=dt_lbm/(dx_lbm^2*Re_wiki_N2);
nu_lbm_HS_N2=dt_lbm/(dx_lbm^2*Re_HS_N2);
nu_lbm_wiki_He=dt_lbm/(dx_lbm^2*Re_wiki_He);
nu_lbm_HS_He=dt_lbm/(dx_lbm^2*Re_HS_He);
tau_wiki_N2=0.5*(6*nu_lbm_wiki_N2+1)
tau_HS_N2=0.5*(6*nu_lbm_HS_N2+1)
tau_wiki_He=0.5*(6*nu_lbm_wiki_He+1);
tau_HS_He=0.5*(6*nu_lbm_HS_He+1);

%% plot viscosity
P_mTorr_array=1000:1000:760000;
P_Pa_array=P_mTorr/7.50061683;
i=1;
for p=P_Pa_array 
nu_wiki_N2_array(i)=nu_wiki(28,2.25e-10,p);
i=i+1;
end 
plot(P_Pa_array,nu_wiki_N2_array,'-g','LineWidth',1.5)
xlabel('Pressure(Pa)','fontsize',14)
ylabel('viscosity(m^2.s^-1)','fontsize',14)
%%
figure
R=8.3144621;%J.k^-1.mol^-1
gamma=1.4;
M_air=29e-3; %Kg.mol^-1
p_ratio=1:-0.01:0.53;
q_ratio=sqrt(p_ratio.^(2/gamma).*(1-p_ratio.^((gamma-1)/gamma)));
q_ratio=q_ratio/(max(q_ratio));
plot(p_ratio,q_ratio,'-r','LineWidth',1.5)
set(gca,'xdir','rev','FontSize',18)
xlabel('P_{2}/P_{1}','fontsize',26)
ylabel('j/jmax','fontsize',26)
filetowrite='C:\Users\coupier\Documents\Ongoing Work\Flow\bernouilli_choked_flow.png';
print('-dpng',filetowrite)
%% energy transfer between collisions toy model
kg_amu = 1.6605402e-27; %(kg/amu) conversion factor
eV_J = 6.2415095e+18;   % (eV/J) conversion factor
k = 1.3806505e-23 ;      % Boltzmann constant (J/K)
T=300; %temp. in Kelvin
cor=0.1; %coef. of restitution
KETh=1.5*k*T*eV_J;
mi=10000; %amu
gas_type='N2';
mg=28; %amu
KE_ion_init=0.01;
KE_gas_init=KETh;
KEs_init=0.01:0.05:10;
KE_ion_final=[];
KE_gas_final=[];
KE_ion_init_ratio=[];
KE_ion_final_ratio=[];
for KE_ion_init=KEs_init
ui=sqrt(2*KE_ion_init/(eV_J*kg_amu*mi)); %m.s-1
ug=sqrt(2*KE_gas_init/(eV_J*kg_amu*mg)); %m.s-1
vi=(ui*(mi-cor*mg)+mg*(1+cor)*ug)/(mi+mg);
vg=(ug*(mg-cor*mi)+mi*(1+cor)*ui)/(mi+mg);
KE_ion_final=[KE_ion_final 0.5*mi*vi^2*kg_amu*eV_J]; %eV
KE_gas_final=[KE_gas_final 0.5*mg*vg^2*kg_amu*eV_J]; %eV
KE_ion_init_ratio=[KE_ion_init_ratio,  KE_ion_init/(KE_ion_init+KE_gas_init)];
%KE_gas_init_ratio=KE_gas_init/(KE_ion_init+KE_gas_init)
KE_ion_final_ratio=[KE_ion_final_ratio,  KE_ion_final/(KE_ion_final+KE_gas_final)];
%KE_gas_final_ratio=KE_gas_final/(KE_ion_final+KE_gas_final)
end
% figure
% [hAx,hLine1,hLine2]=plotyy([KEs_init,KEs_init],[KE_ion_final,KE_gas_final],[KEs_init,KEs_init],[KE_ion_init_ratio,KE_ion_final_ratio]);
% hLine1.LineStyle = '-';
% hLine2.LineStyle = '--';
% ylabel(hAx(1),'Final KE (eV)') % left y-axis
% ylabel(hAx(2),'Ion KE / Total KE') % right y-axis
% xlabel('Initial Ion KE (eV)','fontsize',12)
% legend('Ion Final KE (eV)','Initial Ion KE/Initial Total KE','Gas final KE (eV)','Final Ion KE/Final Total KE','Location','NorthWest')
% title('Kinetic Energy sharing Toy collision model')
figure
hold on 
plot(KEs_init,KE_ion_final,'r-')
plot(KEs_init,KE_gas_final,'b--')
ylabel('Final KE (eV)','fontsize',12)
xlabel('Initial Ion KE (eV)','fontsize',12)
legend('Ion Final KE (eV)','Gas final KE (eV)','Location','best')
title({'1D collision model - Kinetic Energy sharing';[ gas_type ' Gas, COR = ' num2str(cor) ', T=' num2str(T) 'K']})
filetowrite=['C:\Users\coupier\Documents\Ongoing Work\my paper\KEsharing_' gas_type '_cor_' num2str(cor) '_T_' num2str(T) '.png'];
%print('-dpng',filetowrite)
figure
hold on 
plot(KEs_init,KE_ion_init_ratio,'r-')
plot(KEs_init,KE_ion_final_ratio,'b--')
ylim([0,1])
ylabel('KEion/(KEion+KEgas) ratio','fontsize',12)
xlabel('Initial Ion KE (eV)','fontsize',12)
legend('Ion initial KE ratio','Ion final KE ratio','Location','best')
title({'1D collision model - Ion Kinetic Energy';[ gas_type ' Gas, COR = ' num2str(cor) ', T=' num2str(T) 'K']})
filetowrite=['C:\Users\coupier\Documents\Ongoing Work\my paper\IonKE_' gas_type '_cor_' num2str(cor) '_T_' num2str(T) '.png'];
%print('-dpng',filetowrite)
%% collision cross sections comparison
% Langevin
Langevin_prob = 1.6963021e2*sqrt(alpha/(ion_mass*_gas_mass_amu/(ion_mass+_gas_mass_amu)))*(_pressure_pa/_temperature_k)*(ion_time_step)
Langevin_MFP_mm=speed_ion*collision_prob/(ion_time_step)
% Harsphere
-- Compute mean gas speed (mm/us)
				local c_bar_gas = sqrt(8*k*_temperature_k/pi/(_gas_mass_amu * kg_amu)) / 1000

				-- Compute median gas speed (mm/us)
				local c_star_gas = sqrt(2*k*_temperature_k/(_gas_mass_amu * kg_amu)) / 1000

				-- Compute mean relative speed (mm/us) between gas and ion.
				local s = speed_ion / c_star_gas
				local c_bar_rel = c_bar_gas * (
					(s + 1/(2*s)) * 0.5 * sqrt(pi) * erf(s) + 0.5 * exp(-s*s))
HS_effective_MFP_mm = 1000 * k * _temperature_k *(speed_ion / c_bar_rel) / (_pressure_pa * _sigma_m2)
HS_prob= 1 - exp(- speed_ion * ion_time_step / HS_effective_MFP_mm)
% Hybrid


%% collision cross section E, vr dependence
k=1.3806505e-23 ;
T=300;
P_mtorr=500;
P_Pa=P_mtorr/7.50061683
E=100/0.25e-3;% V.m-1
EN=E*k*T/P_Pa
m_ion=78;
m_gas=28;
d_ion_and_gas=428e-12; %m2
lowFieldThreshold=(m_ion/(m_ion+m_gas))^0.5*d_ion_and_gas^2

%% collision cross section langevin hardsphere comparison
k=1.3806505e-23 ;
kg_amu = 1.6605402e-27;
e=1.60217657e-19;
epsilon0=8.854187817e-12;
alpha_N2=1.7403;
alpha_He=0.2050;
r_N2=180e-12;
r_He=140e-12;
r_benzene=250e-12;%?
r_toluene=390e-12;
P_mtorr=500;
P_Pa=P_mtorr/7.50061683;
T=300;
vrel=300; %m/s
sigma_hs_N2=pi*(r_N2+r_toluene)^2;
lang_const=(e/(2*epsilon0))*sqrt(4*pi*epsilon0*1e-30/kg_amu);
sigma_lang_N2=lang_const*sqrt(alpha_N2/(92*28/(92+28)))/vrel;
vrels=1:5000;
sigma_lang_N2_array=lang_const*sqrt(alpha_N2/(92*28/(92+28)))./vrels;
sigma_hs_N2_array=ones(1,5000)*sigma_hs_N2;
sigma_tot_N2_array=sigma_lang_N2_array+sigma_hs_N2_array;
figure
hold on 
plot(vrels,sigma_hs_N2_array,'b-','LineWidth',1)
plot(vrels,sigma_lang_N2_array,'r-','LineWidth',1)
plot(vrels,sigma_tot_N2_array,'k-','LineWidth',2)
ylim([0,0.05e-16])
xlim([100,5000])
ylabel('cs (m^{2})','fontsize',16)
xlabel('v_{rel} (m/s)','fontsize',16)
legend({'Hardsphere ccs','Langevin ccs','Hybrid ccs','Location','best'},'fontsize',16)
title('Hybrid collision cross section','fontsize',16)
filetowrite=['C:\Users\coupier\Documents\Ongoing Work\my paper\discussion with spangler\Hybrid_CCS.png'];
print('-dpng',filetowrite)
figure
hold on
K0_array=2*(92+28)^2/92^2/28./(sigma_hs_N2_array.*vrels);
plot(vrels,K0_array,'k-','LineWidth',1)
xlim([1,5000])
ylabel('K0 (m^{2})','fontsize',12)
xlabel('v_{rel} (m/s)','fontsize',12)
%% Collision center of mass 
% We place ourselves in the reference frame of the gas particle
%vxIon2=vxIon-vxGas
%vyIon2=vyIon-vyGas
%vzIon2=vzIon-vzGas
%Calculate spheric coordinate in the center of mass frame
% speedionr=sqrt(vxIon2^2+vyIon2^2+vzIon2^2) %radial velocity in the  
% azionr=-atan2(vzIon2,vxIon2)
% elionr=(pi/2)-acos(vyIon2/speedionr)
% speedGasr=sqrt(vxGas2^2+vyGas2^2+vzGas2^2) %radial velocity in the  
% azGasr=-atan2(vzGas2,vxGas2)
% elGasr=(pi/2)-acos(vyGas2/speedGasr)
% [VazIon,VelIon,VrIon] = cart2sph(vxIon2,vyIon2,vzIon2);  
% [VazGas,VelGas,VrGas] = cart2sph(vxGas2,vyGas2,vzGas2); 
%equivalent to standard cart2sph transformation
% azionr2=atan2(vyIon2,vxIon2)
% elionr2=(pi/2)-acos(vzIon2/speedionr)
% pick 2 angles for take into account the fact that we collide extended spheres 
% impacttheta=5*2*pi/360;%asin(sqrt(0.999999999*rU(2)));
% impactphi=0*2*pi/360;%2*pi*rU(1);
%determine the position of collision on the gas molecule
%transform minus relative velocity -(vion-vgas) in spherical coordinate:
% nvrel=sqrt((vxIon-vxGas)^2+(vyIon-vyGas)^2+(vzIon-vzGas)^2);
% [VrelAz,VrelEl,VrelR] = cart2sph((vxIon-vxGas)/nvrel,(vyIon-vyGas)/nvrel,(vzIon-vzGas)/nvrel);
%in the ref frame aligned with the rel velocity the rel velocity
%is(0,90,VrelR), then once it is correctly oriented to the impact point 
%in this frame, it becomes (impactphi,90-impacttheta,VrelR). Then,
%we transform the vector back into the original ref frame
%(impactphi+VrelAz,90-impacttheta-90+VrelEl,VrelR). 
% [Xref2x,Yref2x,Zref2x]=sph2cart(impactphi+VrelAz,impacttheta+VrelEl,1);
% [XImpact,YImpact,ZImpact]=sph2cart(impactphi,pi/2-impacttheta,VrelR);
% 
% XImpact2=ZImpact*sin(-VrelEl)+XImpact*cos(-VrelEl);
% YImpact2=YImpact;
% ZImpact2=ZImpact*cos(-VrelEl)-XImpact*sin(-VrelEl);
% XImpact3=XImpact2*cos(VrelAz)-YImpact2*sin(VrelAz);
% YImpact3=XImpact2*sin(VrelAz)+YImpact2*cos(VrelAz);
% ZImpact3=ZImpact2;
% vrion=speedionr*cos(impacttheta);
% vtion=speedionr*sin(impacttheta);
% vrion2=vrion*(mIon-mGas*eRestitution)/(mIon+mGas);
% cangle=cos(twopi/4-impacttheta);sangle=sin(twopi/4-impacttheta);
% vxIon3=vrion2*cangle-vtion*sangle;
% vyIon3=vrion2*sangle+vtion*cangle;
% vzIon3=0;
% ctheta = COS(impactphi); stheta = SIN(impactphi);
% vxIon3=vxIon3*ctheta-vyIon3*stheta;
% vyIon3=vxIon3*stheta+vyIon3*ctheta;
% celback= COS(-twopi/4+elionr); selback = SIN(-twopi/4+elionr);
% vxIon3=vxIon3*celback-vyIon3*selback;
% vyIon3=vxIon3*selback+vyIon3*celback;
% cazback= COS(azionr); sazback = SIN(azionr);
% vxIon3=vxIon3*cazback-vyIon3*sazback;
% vyIon3=vxIon3*sazback+vyIon3*cazback;
% vxIon3=vxIon3+vxGas;
% vyIon3=vyIon3+vyGas;
% vzIon3=vzIon3+vzGas;
% figure
% hold on
% quiver3(0,0,0,100,0,0,0,'k')
% quiver3(0,0,0,0,100,0,0,'k')
% quiver3(0,0,0,0,0,100,0,'k')
% quiver3(0,0,0,vxIon,vyIon,vzIon,'r--')
% quiver3(0,0,0,vxGas,vyGas,vzGas,'b--')
% quiver3(0,0,0,vxIon2,vyIon2,vzIon2,'r-')
% quiver3(0,0,0,vxGas2,vyGas2,vzGas2,'b-')
% xlabel('x')
% ylabel('y')
% zlabel('z')
% text(0.3,0.2,0.5,{['vxIon= ' num2str(vxIon) ', vyIon= ' num2str(vyIon) ',vzIon=' num2str(vzIon)];...
%     ['vxGas= ' num2str(vxGas) ', vyGas= ' num2str(vyGas) ',vzGas=' num2str(vzGas)];...
%     ['uxIon= ' num2str(vxIon2) ', uyIon= ' num2str(vyIon2) ',uzIon=' num2str(vzIon2)];...
%     ['urIon= ' num2str(VrIon) ', uazIon= ' num2str(VazIon*360/(2*pi)) ',uelIon=' num2str(VelIon*360/(2*pi))];...
%     ['uxGas= ' num2str(vxGas2) ', uyGas= ' num2str(vyGas2) ',uzGas=' num2str(vzGas2)];...
%     ['urGas= ' num2str(VrGas) ', uazGas= ' num2str(VazGas*360/(2*pi)) ',uelGas=' num2str(VelGas*360/(2*pi))]},'Units','normalized')
% view(90,0)


vxIon=100; vyIon=0; vzIon=0;
nvIon=sqrt(vxIon^2+vyIon^2+vzIon^2);
vxGas=20; vyGas=150; vzGas=-76; 
mGas=28; mIon=128; cor=1;
k=1.3806505e-23 ;
T=300;
kg_amu = 1.6605402e-27;
iterations=2000;
voutIon=zeros(iterations,6);
for i=1:iterations 

vGas = random('Normal',0,sqrt(k*T/(mGas*kg_amu)),1,3);
rU=random('unif',0,1,1,2);
% We place ourselves in the reference frame of the center of mass
vxIon2=vxIon-(mGas*vxGas+mIon*vxIon)/(mIon+mGas);
vyIon2=vyIon-(mGas*vyGas+mIon*vyIon)/(mIon+mGas);
vzIon2=vzIon-(mGas*vzGas+mIon*vzIon)/(mIon+mGas);
vxGas2=vxGas-(mGas*vxGas+mIon*vxIon)/(mIon+mGas);
vyGas2=vyGas-(mGas*vyGas+mIon*vyIon)/(mIon+mGas);
vzGas2=vzGas-(mGas*vzGas+mIon*vzIon)/(mIon+mGas);

save=false;
impacttheta=asin(sqrt(0.999999999*rU(2))); %30*2*pi/360;
impactphi=2*pi*rU(1); %0*2*pi/360;
nvrel=sqrt((vxIon-vxGas)^2+(vyIon-vyGas)^2+(vzIon-vzGas)^2);
[mVrelAz,mVrelEl,mVrelR]=cart2sph((vxGas-vxIon)/nvrel,(vyGas-vyIon)/nvrel,(vzGas-vzIon)/nvrel);
[XImpact,YImpact,ZImpact]=sph2cart(impactphi,pi/2-impacttheta,mVrelR);
% XImpact=1;
% YImpact=1;
% ZImpact=0;
% nImpact=sqrt((XImpact)^2+(YImpact)^2+(ZImpact)^2);
% XImpact=XImpact/nImpact;
% YImpact=YImpact/nImpact;
% ZImpact=ZImpact/nImpact;
%mVrelEl=-90*2*pi/360;
%mVrelAz=0*2*pi/360;

XImpact2=ZImpact*sin(pi/2-mVrelEl)+XImpact*cos(pi/2-mVrelEl);
YImpact2=YImpact;
ZImpact2=ZImpact*cos(pi/2-mVrelEl)-XImpact*sin(pi/2-mVrelEl);
XImpact3=XImpact2*cos(mVrelAz)-YImpact2*sin(mVrelAz);
YImpact3=XImpact2*sin(mVrelAz)+YImpact2*cos(mVrelAz);
ZImpact3=ZImpact2;
[AzImpact3,ElImpact3,RImpact3]=cart2sph(XImpact3,YImpact3,ZImpact3);
[XIon,YIon,ZIon]=sph2cart(AzImpact3,ElImpact3,RImpact3+2);
%at this point we have determined the collision
%We will now calculate the radial force exerted by each ball on each other
%1 collision is treated using velocities in the center of mass:
%XImpact3 is chosen as norm vector along radial direction
%calculating radial part of velocities in CM
VIonCMrad=-(vxIon2*XImpact3+vyIon2*YImpact3+vzIon2*ZImpact3);
VGasCMrad=-(vxGas2*XImpact3+vyGas2*YImpact3+vzGas2*ZImpact3);
%calculating transverse normalized vector
XtransVector=vxIon2+VIonCMrad*XImpact3;
YtransVector=vyIon2+VIonCMrad*YImpact3;
ZtransVector=vzIon2+VIonCMrad*ZImpact3;
normtransvector=sqrt(XtransVector^2+YtransVector^2+ZtransVector^2);
XtransnormVector=XtransVector/normtransvector;
YtransnormVector=YtransVector/normtransvector;
ZtransnormVector=ZtransVector/normtransvector;
%determining transverse components of velocities in CM
VIonCMtrans=vxIon2*XtransnormVector+vyIon2*YtransnormVector+vzIon2*ZtransnormVector;
VGasCMtrans=vxGas2*XtransnormVector+vyGas2*YtransnormVector+vzGas2*ZtransnormVector;
%Calculating resulting velocities after collision using radial components
%of the force only:
VIonCMradaftercol=(cor*mGas*(VGasCMrad-VIonCMrad)+mIon*VIonCMrad+mGas*VGasCMrad)/(mGas+mIon);
VGasCMradaftercol=(cor*mIon*(-VGasCMrad+VIonCMrad)+mIon*VIonCMrad+mGas*VGasCMrad)/(mGas+mIon);
%Recomposing ion velocities in CM after collision
VxIonCMaftercol=-VIonCMradaftercol*XImpact3+VIonCMtrans*XtransnormVector;
VyIonCMaftercol=-VIonCMradaftercol*YImpact3+VIonCMtrans*YtransnormVector;
VzIonCMaftercol=-VIonCMradaftercol*ZImpact3+VIonCMtrans*ZtransnormVector;
%not necessary to calculate normalized transverse vector since XtransVector
%can be used to recompose the ion vector after collision
%VIonCMtrans*XtransnormVector
%VIonCMtrans*YtransnormVector
%VIonCMtrans*ZtransnormVector
%normalization of ion velocity in CM after collision for visualization:
VIonCMaftercolnorm=sqrt(VxIonCMaftercol^2+VyIonCMaftercol^2+VzIonCMaftercol^2);
VxIonCMaftercolnorm=VxIonCMaftercol/VIonCMaftercolnorm;
VyIonCMaftercolnorm=VyIonCMaftercol/VIonCMaftercolnorm;
VzIonCMaftercolnorm=VzIonCMaftercol/VIonCMaftercolnorm;
%Recomposing gas velocities in CM after collision:
VxGasCMaftercol=-VGasCMradaftercol*XImpact3+VGasCMtrans*XtransnormVector;
VyGasCMaftercol=-VGasCMradaftercol*YImpact3+VGasCMtrans*YtransnormVector;
VzGasCMaftercol=-VGasCMradaftercol*ZImpact3+VGasCMtrans*ZtransnormVector;
%normalization of gas velocity in CM after collision for visualization:
VGasCMaftercolnorm=sqrt(VxGasCMaftercol^2+VyGasCMaftercol^2+VzGasCMaftercol^2);
VxGasCMaftercolnorm=VxGasCMaftercol/VGasCMaftercolnorm;
VyGasCMaftercolnorm=VyGasCMaftercol/VGasCMaftercolnorm;
VzGasCMaftercolnorm=VzGasCMaftercol/VGasCMaftercolnorm;
%Transformation back into the LAB reference frame of the velocities after
%collision:
%Ion
VxIonLABaftercol=VxIonCMaftercol+(mGas*vxGas+mIon*vxIon)/(mIon+mGas);
VyIonLABaftercol=VyIonCMaftercol+(mGas*vyGas+mIon*vyIon)/(mIon+mGas);
VzIonLABaftercol=VzIonCMaftercol+(mGas*vzGas+mIon*vzIon)/(mIon+mGas);
VIonLABaftercolnorm=sqrt(VxIonLABaftercol^2+VyIonLABaftercol^2+VzIonLABaftercol^2);
VxIonLABaftercolnorm=VxIonLABaftercol/VIonLABaftercolnorm;
VyIonLABaftercolnorm=VyIonLABaftercol/VIonLABaftercolnorm;
VzIonLABaftercolnorm=VzIonLABaftercol/VIonLABaftercolnorm;
%Gas molecule
VxGasLABaftercol=VxGasCMaftercol+(mGas*vxGas+mIon*vxIon)/(mIon+mGas);
VyGasLABaftercol=VyGasCMaftercol+(mGas*vyGas+mIon*vyIon)/(mIon+mGas);
VzGasLABaftercol=VzGasCMaftercol+(mGas*vzGas+mIon*vzIon)/(mIon+mGas);
VGasLABaftercolnorm=sqrt(VxGasLABaftercol^2+VyGasLABaftercol^2+VzGasLABaftercol^2);
VxGasLABaftercolnorm=VxGasLABaftercol/VGasLABaftercolnorm;
VyGasLABaftercolnorm=VyGasLABaftercol/VGasLABaftercolnorm;
VzGasLABaftercolnorm=VzGasLABaftercol/VGasLABaftercolnorm;

%SIMION's algorithm
% We place ourselves in the reference frame of the gas molecule
vxIon2sim=vxIon-vxGas;
vyIon2sim=vyIon-vyGas;
vzIon2sim=vzIon-vzGas;
%transformation in spheric coordinates:
speedionr=sqrt(vxIon2sim^2+vyIon2sim^2+vzIon2sim^2);
azionr=-atan2(vzIon2sim,vxIon2sim);
elionr=(pi/2)-acos(vyIon2sim/speedionr);
%pick 2 angles for collision:
% impactangle=90*2*pi/360;%asin(sqrt(0.999999999*rU(2)));
% impacttheta=0*2*pi/360;%2*pi*rU(1);
%projection of radial component onto line of collision:
vrion=speedionr*cos(impacttheta);
vtion=speedionr*sin(impacttheta);
%collision
vrion2=vrion*(mIon-cor*mGas)/(mGas+mIon);
vrgas2=(cor+1)*mIon*vrion/(mGas+mIon);
%reference frames transfo to randomize the collision impact
cangle=cos(pi/2-impacttheta);sangle=sin(pi/2-impacttheta);
vxIon3=vrion2*cangle-vtion*sangle;
vyIon3=vrion2*sangle+vtion*cangle;
vzIon3=0;
ctheta=cos(impactphi);stheta=sin(impactphi);
vxIon4=vxIon3*ctheta+vzIon3*stheta;
vyIon4=vyIon3;
vzIon4=-vxIon3*stheta+vzIon3*ctheta;
celback= cos(-pi/2+elionr); selback = sin(-pi/2+elionr);
vxIon5=vxIon4*celback-vyIon4*selback;
vyIon5=vxIon4*selback+vyIon4*celback;
vzIon5=vzIon4;
cazback= cos(azionr); sazback = sin(azionr);
vxIon6=vxIon5*cazback+vzIon5*sazback;
vyIon6=vyIon5;
vzIon6=-vxIon5*sazback+vzIon5*cazback;
vxIon7=vxIon6+vxGas;
vyIon7=vyIon6+vyGas;
vzIon7=vzIon6+vzGas;
nvIon7=sqrt(vxIon7^2+vyIon7^2+vzIon7^2);

voutIon(i,1:6)=[VxIonLABaftercol, VyIonLABaftercol, VzIonLABaftercol,vxIon7,vyIon7,vzIon7];

end
square=voutIon.^2;
norms=[sqrt(square(:,1)+square(:,2)+square(:,3)),sqrt(square(:,4)+square(:,5)+square(:,6))];
ratios=norms(:,1)./norms(:,2);
angles=[acos(voutIon(:,1:3)*[vxIon; vyIon; vzIon]./(nvIon*norms(:,1))).*(360/(2*pi)),...
    acos(voutIon(:,4:6)*[vxIon; vyIon; vzIon]./(nvIon*norms(:,2))).*(360/(2*pi))];

nbins=20;
figure
[n,xout]=hist(norms(:,1),nbins);
%binsize=(max(eject_times)-min(eject_times))/nbins;
%histnorm=n/(binsize*nb_of_times);
bar(xout,n,'b')

figure
[n,xout]=hist(norms(:,2),nbins);
bar(xout,n,'r')


[x,y,z] = sphere;
x2 = 2*x; y2 = 2*y; z2 = 2*z;
figure
hold on
surf(x,y,z,'FaceAlpha',0.2) % centered at (3,-2,0)
surf(x2+XIon,y2+YIon,z2+ZIon,'FaceAlpha',0.2)
nvIon=sqrt(vxIon^2+vyIon^2+vzIon^2);
nvIon2=sqrt(vxIon2^2+vyIon2^2+vzIon2^2);
nvGas=sqrt(vxGas^2+vyGas^2+vzGas^2);
nvGas2=sqrt(vxGas2^2+vyGas2^2+vzGas2^2);
quiver3(128*XIon/156,128*YIon/156,128*ZIon/156,2,0,0,0,'k','LineWidth',2,'MaxHeadSize',0.4)
quiver3(128*XIon/156,128*YIon/156,128*ZIon/156,0,2,0,0,'k','LineWidth',2,'MaxHeadSize',0.4)
quiver3(128*XIon/156,128*YIon/156,128*ZIon/156,0,0,2,0,'k','LineWidth',2,'MaxHeadSize',0.4)
quiver3(XIon,YIon,ZIon,vxIon/nvIon,vyIon/nvIon,vzIon/nvIon,'r-','LineWidth',2,'MaxHeadSize',0.8)
quiver3(0,0,0,vxGas/nvGas,vyGas/nvGas,vzGas/nvGas,'b-','LineWidth',2,'MaxHeadSize',0.8)
%quiver3(0,0,0,(vxIon-vxGas)/nvrel,(vyIon-vyGas)/nvrel,(vzIon-vzGas)/nvrel,'m--','LineWidth',2,'MaxHeadSize',0.8)
quiver3(128*XIon/156,128*YIon/156,128*ZIon/156,vxIon2/nvIon2,vyIon2/nvIon2,vzIon2/nvIon2,'-','Color',[255/255,128/255,0],'LineWidth',2,'MaxHeadSize',0.8)
quiver3(128*XIon/156,128*YIon/156,128*ZIon/156,vxGas2/nvGas2,vyGas2/nvGas2,vzGas2/nvGas2,'-','Color',[0/255,128/255,255/255],'LineWidth',2,'MaxHeadSize',0.8)
%quiver3(3,0,0,XImpact/nImpact,YImpact/nImpact,ZImpact/nImpact,'c-','LineWidth',2,'MaxHeadSize',0.8)
%quiver3(0,0,0,-(vxIon-vxGas)/nvrel,-(vyIon-vyGas)/nvrel,-(vzIon-vzGas)/nvrel,'y--','LineWidth',2,'MaxHeadSize',0.8)
quiver3(XImpact3,YImpact3,ZImpact3,(vxIon-vxGas)/nvrel,(vyIon-vyGas)/nvrel,(vzIon-vzGas)/nvrel,'g-.','LineWidth',2,'MaxHeadSize',0.8)
%quiver3(0,0,0,XImpact,YImpact,ZImpact,'m-','LineWidth',2,'MaxHeadSize',0.8)
quiver3(0,0,0,XIon,YIon,ZIon,'k:','LineWidth',2,'MaxHeadSize',0.8,'ShowArrowHead','off')
%quiver3(XImpact3,YImpact3,ZImpact3,XtransnormVector,YtransnormVector,ZtransnormVector,'y--','LineWidth',2,'MaxHeadSize',0.8)
quiver3(128*XIon/156,128*YIon/156,128*ZIon/156,VxIonCMaftercolnorm,VyIonCMaftercolnorm,VzIonCMaftercolnorm,'--','Color',[255/255,128/255,0],'LineWidth',2,'MaxHeadSize',0.8)
quiver3(128*XIon/156,128*YIon/156,128*ZIon/156,VxGasCMaftercolnorm,VyGasCMaftercolnorm,VzGasCMaftercolnorm,'--','Color',[0/255,128/255,255/255],'LineWidth',2,'MaxHeadSize',0.8)
quiver3(XIon,YIon,ZIon,VxIonLABaftercolnorm,VyIonLABaftercolnorm,VzIonLABaftercolnorm,'r--','LineWidth',2,'MaxHeadSize',0.9)
quiver3(0,0,0,VxGasLABaftercolnorm,VyGasLABaftercolnorm,VzGasLABaftercolnorm,'b--','LineWidth',2,'MaxHeadSize',0.9)
quiver3(XIon,YIon,ZIon,vxIon7/nvIon7,vyIon7/nvIon7,vzIon7/nvIon7,'m--','LineWidth',2,'MaxHeadSize',0.9)
xlabel('x')
ylabel('y')
zlabel('z')
view(105,28)
if save
filetowrite=['C:\Users\coupier\Documents\Ongoing Work\my paper\kinematics\sphere_theta_' num2str(impacttheta*360/(2*pi))...
    '_phi_' num2str(impactphi*360/(2*pi)) '_cor_' num2str(cor) '.png'];
print('-dpng',filetowrite)
end
%second figure with real velocities (not normalized)
[x,y,z] = sphere;
x2 = 2*x; y2 = 2*y; z2 = 2*z;
figure
hold on
surf(x,y,z,'FaceAlpha',0.2) % centered at (3,-2,0)
surf(x2+XIon,y2+YIon,z2+ZIon,'FaceAlpha',0.2)
quiver3(128*XIon/156,128*YIon/156,128*ZIon/156,2,0,0,0,'k','LineWidth',2,'MaxHeadSize',0.4)
quiver3(128*XIon/156,128*YIon/156,128*ZIon/156,0,2,0,0,'k','LineWidth',2,'MaxHeadSize',0.4)
quiver3(128*XIon/156,128*YIon/156,128*ZIon/156,0,0,2,0,'k','LineWidth',2,'MaxHeadSize',0.4)
quiver3(XIon,YIon,ZIon,vxIon,vyIon,vzIon,'r-','LineWidth',3,'MaxHeadSize',0.8)
quiver3(0,0,0,vxGas,vyGas,vzGas,'b-','LineWidth',3,'MaxHeadSize',0.8)
%quiver3(0,0,0,(vxIon-vxGas),(vyIon-vyGas),(vzIon-vzGas),'m--','LineWidth',2,'MaxHeadSize',0.8)
quiver3(128*XIon/156,128*YIon/156,128*ZIon/156,vxIon2,vyIon2,vzIon2,'-','Color',[255/255,128/255,0],'LineWidth',2,'MaxHeadSize',0.8)
quiver3(128*XIon/156,128*YIon/156,128*ZIon/156,vxGas2,vyGas2,vzGas2,'-','Color',[0/255,128/255,255/255],'LineWidth',2,'MaxHeadSize',0.8)
%quiver3(3,0,0,XImpact/nImpact,YImpact/nImpact,ZImpact/nImpact,'c-','LineWidth',2,'MaxHeadSize',0.8)
%quiver3(0,0,0,-(vxIon-vxGas)/nvrel,-(vyIon-vyGas)/nvrel,-(vzIon-vzGas)/nvrel,'y--','LineWidth',2,'MaxHeadSize',0.8)
quiver3(XImpact3,YImpact3,ZImpact3,(vxIon-vxGas),(vyIon-vyGas),(vzIon-vzGas),'g-.','LineWidth',2,'MaxHeadSize',0.8)
%quiver3(0,0,0,XImpact,YImpact,ZImpact,'m-','LineWidth',2,'MaxHeadSize',0.8)
quiver3(0,0,0,XIon,YIon,ZIon,'k:','LineWidth',2,'MaxHeadSize',0.8,'ShowArrowHead','off')
%quiver3(XImpact3,YImpact3,ZImpact3,XtransnormVector,YtransnormVector,ZtransnormVector,'y--','LineWidth',2,'MaxHeadSize',0.8)
quiver3(128*XIon/156,128*YIon/156,128*ZIon/156,VxIonCMaftercol,VyIonCMaftercol,VzIonCMaftercol,'--','Color',[255/255,128/255,0],'LineWidth',2,'MaxHeadSize',0.8)
quiver3(128*XIon/156,128*YIon/156,128*ZIon/156,VxGasCMaftercol,VyGasCMaftercol,VzGasCMaftercol,'--','Color',[0/255,128/255,255/255],'LineWidth',2,'MaxHeadSize',0.8)
quiver3(XIon,YIon,ZIon,VxIonLABaftercol,VyIonLABaftercol,VzIonLABaftercol,'r--','LineWidth',3,'MaxHeadSize',0.8)
quiver3(0,0,0,VxGasLABaftercol,VyGasLABaftercol,VzGasLABaftercol,'b--','LineWidth',3,'MaxHeadSize',0.8)
quiver3(XIon,YIon,ZIon,vxIon7,vyIon7,vzIon7,'m--','LineWidth',2,'MaxHeadSize',0.9)
xlabel('x')
ylabel('y')
zlabel('z')
view(105,28)
if save
filetowrite=['C:\Users\coupier\Documents\Ongoing Work\my paper\kinematics\realspeed_theta_' num2str(impacttheta*360/(2*pi))...
    '_phi_' num2str(impactphi*360/(2*pi)) '_cor_' num2str(cor) '.png'];
print('-dpng',filetowrite)
end
%% SIMION collision model 
vxIon=-10; vyIon=-30; vzIon=84; 
vxGas=20; vyGas=150; vzGas=-76; 
mGas=28; mIon=128; cor=1;
rU=random('unif',0,1,1,2);
save=false;
% We place ourselves in the reference frame of the gas molecule
vxIon2=vxIon-vxGas
vyIon2=vyIon-vyGas
vzIon2=vzIon-vzGas
%transformation in spheric coordinates:
speedionr=sqrt(vxIon2^2+vyIon2^2+vzIon2^2)
azionr=-atan2(vzIon2,vxIon2)
elionr=(pi/2)-acos(vyIon2/speedionr)
%pick 2 angles for collision:
impactangle=90*2*pi/360;%asin(sqrt(0.999999999*rU(2)));
impacttheta=0*2*pi/360;%2*pi*rU(1);
%projection of radial component onto line of collision:
vrion=speedionr*cos(impactangle)
vtion=speedionr*sin(impactangle)
%collision
vrion2=vrion*(mIon-cor*mGas)/(mGas+mIon)
vrgas2=(cor+1)*mIon*vrion/(mGas+mIon);
%reference frames transfo to randomize the collision impact
cangle=cos(pi/2-impactangle);sangle=sin(pi/2-impactangle);
vxIon3=vrion2*cangle-vtion*sangle;
vyIon3=vrion2*sangle+vtion*cangle;
vzIon3=0;
ctheta=cos(impacttheta);stheta=sin(impacttheta);
vxIon4=vxIon3*ctheta+vzIon3*stheta;
vyIon4=vyIon3;
vzIon4=-vxIon3*stheta+vzIon3*ctheta;
celback= cos(-pi/2+elionr); selback = sin(-pi/2+elionr);
vxIon5=vxIon4*celback-vyIon4*selback;
vyIon5=vxIon4*selback+vyIon4*celback;
vzIon5=vzIon4;
cazback= cos(azionr); sazback = sin(azionr);
vxIon6=vxIon5*cazback+vzIon5*sazback;
vyIon6=vyIon5;
vzIon6=-vxIon5*sazback+vzIon5*cazback;
vxIon7=vxIon6+vxGas;
vyIon7=vyIon6+vyGas;
vzIon7=vzIon6+vzGas;



[x,y,z] = sphere;
x2 = 2*x; y2 = 2*y; z2 = 2*z;
figure
hold on
surf(x,y,z,'FaceAlpha',0.2) % centered at (3,-2,0)
%surf(x2+XIon,y2+YIon,z2+ZIon,'FaceAlpha',0.2)
nvIon=sqrt(vxIon^2+vyIon^2+vzIon^2);
nvIon2=sqrt(vxIon2^2+vyIon2^2+vzIon2^2);
nvGas=sqrt(vxGas^2+vyGas^2+vzGas^2);
nvGas2=sqrt(vxGas2^2+vyGas2^2+vzGas2^2);
%quiver3(128*XIon/156,128*YIon/156,128*ZIon/156,2,0,0,0,'k','LineWidth',2,'MaxHeadSize',0.4)
%quiver3(128*XIon/156,128*YIon/156,128*ZIon/156,0,2,0,0,'k','LineWidth',2,'MaxHeadSize',0.4)
%quiver3(128*XIon/156,128*YIon/156,128*ZIon/156,0,0,2,0,'k','LineWidth',2,'MaxHeadSize',0.4)
quiver3(0,0,0,vxIon/nvIon,vyIon/nvIon,vzIon/nvIon,'r-','LineWidth',2,'MaxHeadSize',0.8)
quiver3(0,0,0,vxGas/nvGas,vyGas/nvGas,vzGas/nvGas,'b-','LineWidth',2,'MaxHeadSize',0.8)
%quiver3(0,0,0,(vxIon-vxGas)/nvrel,(vyIon-vyGas)/nvrel,(vzIon-vzGas)/nvrel,'m--','LineWidth',2,'MaxHeadSize',0.8)
quiver3(128*XIon/156,128*YIon/156,128*ZIon/156,vxIon2/nvIon2,vyIon2/nvIon2,vzIon2/nvIon2,'-','Color',[255/255,128/255,0],'LineWidth',2,'MaxHeadSize',0.8)
%quiver3(3,0,0,XImpact/nImpact,YImpact/nImpact,ZImpact/nImpact,'c-','LineWidth',2,'MaxHeadSize',0.8)
%quiver3(0,0,0,-(vxIon-vxGas)/nvrel,-(vyIon-vyGas)/nvrel,-(vzIon-vzGas)/nvrel,'y--','LineWidth',2,'MaxHeadSize',0.8)
%quiver3(XImpact3,YImpact3,ZImpact3,(vxIon-vxGas)/nvrel,(vyIon-vyGas)/nvrel,(vzIon-vzGas)/nvrel,'g-.','LineWidth',2,'MaxHeadSize',0.8)
%quiver3(0,0,0,XImpact,YImpact,ZImpact,'m-','LineWidth',2,'MaxHeadSize',0.8)
quiver3(0,0,0,XIon,YIon,ZIon,'k:','LineWidth',2,'MaxHeadSize',0.8,'ShowArrowHead','off')
%quiver3(XImpact3,YImpact3,ZImpact3,XtransnormVector,YtransnormVector,ZtransnormVector,'y--','LineWidth',2,'MaxHeadSize',0.8)
quiver3(128*XIon/156,128*YIon/156,128*ZIon/156,VxIonCMaftercolnorm,VyIonCMaftercolnorm,VzIonCMaftercolnorm,'--','Color',[255/255,128/255,0],'LineWidth',2,'MaxHeadSize',0.8)
quiver3(128*XIon/156,128*YIon/156,128*ZIon/156,VxGasCMaftercolnorm,VyGasCMaftercolnorm,VzGasCMaftercolnorm,'--','Color',[0/255,128/255,255/255],'LineWidth',2,'MaxHeadSize',0.8)
quiver3(XIon,YIon,ZIon,VxIonLABaftercolnorm,VyIonLABaftercolnorm,VzIonLABaftercolnorm,'r--','LineWidth',2,'MaxHeadSize',0.9)
quiver3(0,0,0,VxGasLABaftercolnorm,VyGasLABaftercolnorm,VzGasLABaftercolnorm,'b--','LineWidth',2,'MaxHeadSize',0.9)
xlabel('x')
ylabel('y')
zlabel('z')
view(105,28)
if save
filetowrite=['C:\Users\coupier\Documents\Ongoing Work\my paper\kinematics\sphere_theta_' num2str(impacttheta*360/(2*pi))...
    '_phi_' num2str(impactphi*360/(2*pi)) '_cor_' num2str(cor) '.png'];
print('-dpng',filetowrite)
end
%second figure with real velocities (not normalized)
[x,y,z] = sphere;
x2 = 2*x; y2 = 2*y; z2 = 2*z;
figure
hold on
surf(x,y,z,'FaceAlpha',0.2) % centered at (3,-2,0)
surf(x2+XIon,y2+YIon,z2+ZIon,'FaceAlpha',0.2)
quiver3(128*XIon/156,128*YIon/156,128*ZIon/156,2,0,0,0,'k','LineWidth',2,'MaxHeadSize',0.4)
quiver3(128*XIon/156,128*YIon/156,128*ZIon/156,0,2,0,0,'k','LineWidth',2,'MaxHeadSize',0.4)
quiver3(128*XIon/156,128*YIon/156,128*ZIon/156,0,0,2,0,'k','LineWidth',2,'MaxHeadSize',0.4)
quiver3(XIon,YIon,ZIon,vxIon,vyIon,vzIon,'r-','LineWidth',3,'MaxHeadSize',0.8)
quiver3(0,0,0,vxGas,vyGas,vzGas,'b-','LineWidth',3,'MaxHeadSize',0.8)
%quiver3(0,0,0,(vxIon-vxGas),(vyIon-vyGas),(vzIon-vzGas),'m--','LineWidth',2,'MaxHeadSize',0.8)
quiver3(128*XIon/156,128*YIon/156,128*ZIon/156,vxIon2,vyIon2,vzIon2,'-','Color',[255/255,128/255,0],'LineWidth',2,'MaxHeadSize',0.8)
quiver3(128*XIon/156,128*YIon/156,128*ZIon/156,vxGas2,vyGas2,vzGas2,'-','Color',[0/255,128/255,255/255],'LineWidth',2,'MaxHeadSize',0.8)
%quiver3(3,0,0,XImpact/nImpact,YImpact/nImpact,ZImpact/nImpact,'c-','LineWidth',2,'MaxHeadSize',0.8)
%quiver3(0,0,0,-(vxIon-vxGas)/nvrel,-(vyIon-vyGas)/nvrel,-(vzIon-vzGas)/nvrel,'y--','LineWidth',2,'MaxHeadSize',0.8)
quiver3(XImpact3,YImpact3,ZImpact3,(vxIon-vxGas),(vyIon-vyGas),(vzIon-vzGas),'g-.','LineWidth',2,'MaxHeadSize',0.8)
%quiver3(0,0,0,XImpact,YImpact,ZImpact,'m-','LineWidth',2,'MaxHeadSize',0.8)
quiver3(0,0,0,XIon,YIon,ZIon,'k:','LineWidth',2,'MaxHeadSize',0.8,'ShowArrowHead','off')
%quiver3(XImpact3,YImpact3,ZImpact3,XtransnormVector,YtransnormVector,ZtransnormVector,'y--','LineWidth',2,'MaxHeadSize',0.8)
quiver3(128*XIon/156,128*YIon/156,128*ZIon/156,VxIonCMaftercol,VyIonCMaftercol,VzIonCMaftercol,'--','Color',[255/255,128/255,0],'LineWidth',2,'MaxHeadSize',0.8)
quiver3(128*XIon/156,128*YIon/156,128*ZIon/156,VxGasCMaftercol,VyGasCMaftercol,VzGasCMaftercol,'--','Color',[0/255,128/255,255/255],'LineWidth',2,'MaxHeadSize',0.8)
quiver3(XIon,YIon,ZIon,VxIonLABaftercol,VyIonLABaftercol,VzIonLABaftercol,'r--','LineWidth',3,'MaxHeadSize',0.8)
quiver3(0,0,0,VxGasLABaftercol,VyGasLABaftercol,VzGasLABaftercol,'b--','LineWidth',3,'MaxHeadSize',0.8)
xlabel('x')
ylabel('y')
zlabel('z')
view(105,28)
if save
filetowrite=['C:\Users\coupier\Documents\Ongoing Work\my paper\kinematics\realspeed_theta_' num2str(impacttheta*360/(2*pi))...
    '_phi_' num2str(impactphi*360/(2*pi)) '_cor_' num2str(cor) '.png'];
print('-dpng',filetowrite)
end

%% orientation tests
%to reproduce elevation_rotate and azimuth_rotate of SIMION

vx=1/sqrt(2)
vy=1
vz=1/sqrt(2)
[vx2,vy2,vz2]=elevation_rotate(45,vx,vy,vz)
[vx3,vy3,vz3]=elevation_rotate_manura(45,vx,vy,vz)
%[vx2,vy2,vz2]=azimuth_rotate(90,vx,vy,vz)
figure
hold on
quiver3(0,0,0,vx,vy,vz,'k-','LineWidth',1,'MaxHeadSize',0.3)
quiver3(0,0,0,vx2,vy2,vz2,'r-','LineWidth',2,'MaxHeadSize',0.3)
quiver3(0,0,0,vx3,vy3,vz3,'b-','LineWidth',2,'MaxHeadSize',0.3)
xlabel('x')
ylabel('y')
zlabel('z')