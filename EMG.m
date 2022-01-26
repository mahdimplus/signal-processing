clc
clear 
close all
%% load audio signal 
d1=load('E:\master\signal processing progect\EMG_Wheelchair_Rmax.mat');
d2= load('E:\master\signal processing progect\EMG_Wheelchair_Rmin.mat');
%% T
t=size(d1.Times);
t1=size(d1.Biceps);
fs=t*1000/t1;
                       % ITS THE SAME WITH FS
% t1=size(d2.Temps);
% t2=size(d2.Biceps);
% fs1=t1*1000/t2;  
%% plot
Biceps_max=(d1.Biceps);
Triceps_max=(d1.Triceps);
Biceps_min=(d2.Biceps);
Triceps_min=(d2.Triceps);
subplot(2,2,1);
plot(d1.Biceps);
title('Biceps Rmax');
subplot(2,2,3);
plot(d1.Triceps);
title('Triceps Rmax');
subplot(2,2,2);
plot(d2.Biceps);
title(' Biceps Rmin ');
subplot(2,2,4);
plot(d2.Triceps);
title('Triceps Rmin');
%% correlation
figure 
r_Rmax =xcorr(Biceps_max  , Triceps_max );
r_Rmin =xcorr(Triceps_min , Biceps_min);
subplot(2,1,1);
plot(r_Rmax);
title('correlation of Biceps max & Biceps min');
grid on;
subplot(2,1,2);
plot(r_Rmin);
title('correlation of Triceps max & Triceps min');
%% pwelch
 figure
 Nfft = 1024;
[Biceps_max_welch,f1] = pwelch(d1.Biceps,gausswin(Nfft),Nfft/2,Nfft,fs);
[Triceps_max_welch,f2] = pwelch(d1.Triceps,gausswin(Nfft),Nfft/2,Nfft,fs);
[Biceps_min_welch,f3] = pwelch(d2.Biceps,gausswin(Nfft),Nfft/2,Nfft,fs);
[Triceps_min_welch,f4] = pwelch(d2.Triceps,gausswin(Nfft),Nfft/2,Nfft,fs);

subplot(4,1,1);
plot(f1,Biceps_max_welch);
ylabel('PSD ');
xlabel('f (Hz)');
title('Biceps max ');
grid on;

subplot(4,1,2);
plot(f2,Triceps_max_welch);
ylabel('PSD ');
xlabel('f (Hz)');
title('Triceps max  ');
grid on;

subplot(4,1,3);
plot(f3,Biceps_min_welch);
ylabel('PSD ');
xlabel('f (Hz)');
title('Biceps min ');
grid on;

subplot(4,1,4);
plot(f4,Triceps_min_welch);
ylabel('PSD ');
xlabel('f (Hz)');
title('Triceps min  ');
grid on;

 %% low path filter
order=2;
fcuthigh=200;
[b,a]=butter(order,fcuthigh/(fs/2),'low');   
 Biceps_max1=filter(b,a,Biceps_max);
 Biceps_min1=filter(b,a,Biceps_min);
 Triceps_max1=filter(b,a,Triceps_max);
 Triceps_min1=filter(b,a,Triceps_min);
%% downsample
Biceps_max2=downsample(Biceps_max1,2);
Biceps_min2=downsample(Biceps_min1,2);
Triceps_max2=downsample(Triceps_max1,2);
Triceps_min2=downsample(Triceps_min1,2);
%% envelope
Biceps_max3=Biceps_max2.^2;
Biceps_min3=Biceps_min2.^2;
Triceps_max3=Triceps_max2.^2;
Triceps_min3=Triceps_min2.^2;

filtsig1=filter(b,a,Biceps_max3);
filtsig2=filter(b,a,Biceps_min3);
filtsig3=filter(b,a,Triceps_max3);
filtsig4=filter(b,a,Triceps_min3);

[BicepsMaxWelch1,f5] = pwelch(filtsig1,gausswin(Nfft),Nfft/2,Nfft,fs);
[BicepsMinWelch2,f6] = pwelch(filtsig2,gausswin(Nfft),Nfft/2,Nfft,fs);
[TricepsMaxWelch3,f7] = pwelch(filtsig3,gausswin(Nfft),Nfft/2,Nfft,fs);
[TricepsMixWelch4,f8] = pwelch(filtsig4,gausswin(Nfft),Nfft/2,Nfft,fs);
figure
subplot(2,2,1);
plot(f5,BicepsMaxWelch1);
title('BicepsMaxWelch1');
subplot(2,2,2);
plot(f6,BicepsMinWelch2)
title('BicepsMinWelch2');
subplot(2,2,3);
plot(f7,TricepsMaxWelch3)
title('TricepsMaxWelch3');
subplot(2,2,4);
plot(f8,TricepsMixWelch4)
title('TricepsMixWelch4');

order=2;        %order of filter
fcuthigh=200;   %high cut frequency

filtsig5=filter(b,a,Biceps_max2);
filtsig6=filter(b,a,Biceps_min2);
filtsig7=filter(b,a,Triceps_max2);
filtsig8=filter(b,a,Triceps_min2);

[BicepsMaxWelch5,f9] = pwelch(filtsig5,gausswin(Nfft),Nfft/2,Nfft,fs);
[BicepsMinWelch6,f10] = pwelch(filtsig6,gausswin(Nfft),Nfft/2,Nfft,fs);
[TricepsMaxWelch7,f11] = pwelch(filtsig7,gausswin(Nfft),Nfft/2,Nfft,fs);
[TricepsMixWelch8,f12] = pwelch(filtsig8,gausswin(Nfft),Nfft/2,Nfft,fs);

figure
subplot(2,2,1);
plot(f9,BicepsMaxWelch5);
title('BicepsMaxWelch5');
subplot(2,2,2);
plot(f10,BicepsMinWelch6)
title('BicepsMinWelch6');
subplot(2,2,3);
plot(f11,TricepsMaxWelch7)
title('TricepsMaxWelch7');
subplot(2,2,4);
plot(f12,TricepsMixWelch8)
title('TricepsMixWelch8');  

%% co_contraction
maxamplitudebimax = max(filtsig1)-min(filtsig1);
maxamplitudetrimax = max(filtsig3)-min(filtsig3);
maxamplitudebimaxmin= max(filtsig2)-min(filtsig2);
maxamplitudetrimaxmin= max(filtsig4)-min(filtsig4);

e1=find( filtsig1 >.5*max(filtsig1)); %% Biceps_max3
o1=size(e1);

e2=find( filtsig2 >.5*max(filtsig2));% Biceps_min3
o2=size(e2);

e3=find( filtsig3 >.5*max(filtsig3)); % triceps_max3
o3=size(e3);

e4=find( filtsig4 >.5*max(filtsig4));%triceps_min3
o4=size(e4);

c=0;
p=0;
 for i=1:o1(1)
    for j=1:o3(1)
         if  e1(i)==e3(j)
         c=c+1;
         end
     end
 end
 for i1=1:o2(1)
    for j1=1:o4(1)
         if  e2(i1)==e4(j1)
         p=p+1;
         end
     end
 end
w=c/29921;
w1=p/31632;   
%% pre need
cyclemax=fix(d1.begin_cycle*1000);
cyclemax(33)=29921;
cyclemin=fix(d2.begin_cycle*1000);
cyclemin(31)=31632;
c1=0 ;
cb1=0;
cc1=0;
cd1=0;
pb1=0;
pc1=0;
pd1=0;
p1=0 ;
wcycle_50 =zeros(1,32);
wcycle_40 =zeros(1,32);
wcycle_30 =zeros(1,32);
wcycle_20 =zeros(1,32);
w1cycle_50=zeros(1,30);
w1cycle_40=zeros(1,30);
w1cycle_30=zeros(1,30);
w1cycle_20=zeros(1,30);
energy1_50=zeros(1,32);
energy1_40=zeros(1,32);
energy1_30=zeros(1,32);
energy1_20=zeros(1,32);
energy3_50=zeros(1,32);
energy3_40=zeros(1,32);
energy3_30=zeros(1,32);
energy3_20=zeros(1,32);
energy2_50=zeros(1,30);
energy2_40=zeros(1,30);
energy2_30=zeros(1,30);
energy2_20=zeros(1,30);
energy4_50=zeros(1,30);
energy4_40=zeros(1,30);
energy4_30=zeros(1,30);
energy4_20=zeros(1,30);

Biceps_max3cycle=Biceps_max1.^2;
Biceps_min3cycle=Biceps_min1.^2;
Triceps_max3cycle=Triceps_max1.^2;
Triceps_min3cycle=Triceps_min1.^2;

filtsig1cycle=filter(b,a,Biceps_max3cycle);
filtsig2cycle=filter(b,a,Biceps_min3cycle);
filtsig3cycle=filter(b,a,Triceps_max3cycle);
filtsig4cycle=filter(b,a,Triceps_min3cycle);

for x1=1:32
D1=cyclemax(x1+1)-cyclemax(x1);
filt1_slice=zeros(D1,1);
filt3_slice=zeros(D1,1);

for m=cyclemax(x1):cyclemax(x1+1)

    filt1_slice(m)=filtsig1cycle(m);  %biceps max
    filt3_slice(m)=filtsig3cycle(m);  %triceps max
   
end 

maxamplitudebimax  = max(filt1_slice)-min(filt1_slice);
maxamplitudetrimax = max(filt3_slice)-min(filt3_slice);

e1cycle_50=find(filt1_slice >.5*max(filt1_slice)); % Biceps_max3
e1cycle_40=find(filt1_slice >.4*max(filt1_slice));
e1cycle_30=find(filt1_slice >.3*max(filt1_slice));
e1cycle_20=find(filt1_slice >.2*max(filt1_slice));

e3cycle_50=find(filt3_slice >.5*max(filt3_slice)); % triceps_max3
e3cycle_40=find(filt3_slice >.4*max(filt3_slice));
e3cycle_30=find(filt3_slice >.3*max(filt3_slice));
e3cycle_20=find(filt3_slice >.2*max(filt3_slice));

%% energy bic max
energy1_50(x1)=sum(e1cycle_50.^2)*1e-9;
energy1_40(x1)=sum(e1cycle_40.^2)*1e-9;
energy1_30(x1)=sum(e1cycle_30.^2)*1e-9;
energy1_20(x1)=sum(abs(e1cycle_20).^2)*1e-9;
%% energy tri max
energy3_50(x1)=sum(e3cycle_50.^2);
energy3_40(x1)=sum(e3cycle_40.^2);
energy3_30(x1)=sum(e3cycle_30.^2);
energy3_20(x1)=sum(e3cycle_20.^2);

o1max=size(e1cycle_50);
o3max=size(e3cycle_50);

o5max=size(e1cycle_40);
o8max=size(e3cycle_40);

o6max=size(e1cycle_30);
o9max=size(e3cycle_30);

o7max=size(e1cycle_20);
o0max=size(e3cycle_20);

for z=1:o1max(1)
    for z1=1:o3max(1)
         if  e1cycle_50(z)==e3cycle_50(z1)
         c1=c1+1;
         end
    end
end
for z=1:o5max(1)
    for z1=1:o8max(1)
         if  e1cycle_40(z)==e3cycle_40(z1)
         cb1=cb1+1;
         end
    end
end
for z=1:o6max(1)
    for z1=1:o9max(1)
         if  e1cycle_30(z)==e3cycle_30(z1)
         cc1=cc1+1;
         end
    end
end
for z=1:o7max(1)
    for z1=1:o0max(1)
         if  e1cycle_20(z)==e3cycle_20(z1)
         cd1=cd1+1;
         end
    end
end
wcycle_50(x1) =(c1/29921)*1e3;%max
wcycle_40(x1) =cb1/29921*1e3; %max
wcycle_30(x1) =cc1/29921*1e3; %max
wcycle_20(x1) =cd1/29921*1e3; %max
end 
%% min
for s=1:29
    D1=cyclemin(s+1)-cyclemin(s);
    filt2_slice=zeros(D1,1);
    filt4_slice=zeros(D1,1);
    for m1=cyclemin(s):cyclemin(s+1)
           filt2_slice(m1)=filtsig2cycle(m1);  %biceps min
           filt4_slice(m1)=filtsig4cycle(m1);  %triceps min 
    end 

 maxamplitudebimin  = max(filt2_slice)-min(filt2_slice);
 maxamplitudetrimin = max(filt4_slice)-min(filt4_slice);

e2cycle_50=find( filt2_slice >.5*max(filt2_slice)); % Biceps_min3
e2cycle_40=find( filt2_slice >.4*max(filt2_slice)); % Biceps_min3
e2cycle_30=find( filt2_slice >.3*max(filt2_slice)); % Biceps_min3
e2cycle_20=find( filt2_slice >.2*max(filt2_slice)); % Biceps_min3
e4cycle_50=find( filt4_slice >.5*max(filt4_slice)); %triceps_min3
e4cycle_40=find( filt4_slice >.4*max(filt4_slice)); %triceps_min3
e4cycle_30=find( filt4_slice >.3*max(filt4_slice)); %triceps_min3
e4cycle_20=find( filt4_slice >.2*max(filt4_slice)); %triceps_min3

%% energy bic min
energy2_50(s)=sum(e2cycle_50.^2)*1e-9;
energy2_40(s)=sum(e2cycle_40.^2)*1e-9;
energy2_30(s)=sum(e2cycle_30.^2)*1e-9;
energy2_20(s)=sum(abs(e2cycle_20).^2)*1e-9;
%% energy tri min
energy4_50(s)=sum(e4cycle_50.^2);
energy4_40(s)=sum(e4cycle_40.^2);
energy4_30(s)=sum(e4cycle_30.^2);
energy4_20(s)=sum(e4cycle_20.^2);

o2min=size(e2cycle_50);
o4min=size(e4cycle_50);

o1min=size(e2cycle_50);
o3min=size(e4cycle_50);

o5min=size(e2cycle_50);
o6min=size(e4cycle_50);

o7min=size(e2cycle_50);
o8min=size(e4cycle_50);
 for z2=1:o2min(1)
    for z3=1:o4min(1)
         if  e2cycle_50(z2)==e4cycle_50(z3)
         p1=p1+1;
         end
    end
  end
  for z2=1:o1min(1)
    for z3=1:o3min(1)
         if  e2cycle_40(z2)==e4cycle_40(z3)
         pb1=pb1+1;
         end
    end
  end
  
  for z2=1:o5min(1)
    for z3=1:o6min(1)
         if  e2cycle_30(z2)==e4cycle_30(z3)
         pc1=pc1+1;
         end
    end
  end
  for z2=1:o7min(1)
    for z3=1:o8min(1)
         if  e2cycle_20(z2)==e4cycle_20(z3)
         pd1=pd1+1;
         end
    end
  end
w1cycle_50(s)=(p1/31632)*1e3;%min
w1cycle_40(s)=pb1/31632*1e3; %min
w1cycle_30(s)=pc1/31632*1e3; %min
w1cycle_20(s)=pd1/31632*1e3; %min
end   
%% saving C1
save('resultc1_50','wcycle_50','energy1_50','energy3_50')
save('resultc1_40','wcycle_40','energy1_40','energy3_40')
save('resultc1_30','wcycle_30','energy1_30','energy3_30')
save('resultc1_20','wcycle_20','energy1_20','energy3_20')
%% saving C2
save('resultc2_50','wcycle_50','energy2_50','energy4_50')
save('resultc2_40','wcycle_40','energy2_40','energy4_40')
save('resultc2_30','wcycle_30','energy2_30','energy4_30')
save('resultc2_20','wcycle_20','energy2_20','energy4_20')
%% classification 50%
uy_max50 = mean(wcycle_50);
ux_max50 = mean(energy1_50);
uy_min50=mean(w1cycle_50);
ux_min50=mean(energy2_50);
y50=[uy_max50 uy_min50 ];
x50=[ux_max50 ux_min50 ]; 
m50=(uy_max50-uy_min50)/(ux_max50-ux_min50);
x_50middle=(ux_max50+ux_min50)/2;
y_50middle=(uy_max50+uy_min50)/2;
xq1_50=0:1:10;
yq1_50=y_50middle+m50*(xq1_50 - x_50middle);
m50new=-(1/m50);
x50new= 0:0.001:10;
y50new=y_50middle+m50new*(x50new-x_50middle);
%% classification 40%
uy_max40 = mean(wcycle_40);
ux_max40 = mean(energy1_40)-2; 
uy_min40=mean(w1cycle_40);
ux_min40=mean(energy2_40);

y40=[uy_max40 uy_min40 ];
x40=[ux_max40 ux_min40 ];

m40=(uy_max40-uy_min40)/(ux_max40-ux_min40);
x_40middle=(ux_max40+ux_min40)/2;
y_40middle=(uy_max40+uy_min40)/2;

xq1_40=0:1:10;
yq1_40=y_40middle+m40*(xq1_40-x_40middle);

m40new=-(1/m40);
x40new= (0:0.0001:20);
y40new=y_40middle+m40new*(x40new-x_40middle);

%% classification 30%
uy_max30 = mean(wcycle_30);
ux_max30 = mean(energy1_30)-5;
uy_min30=mean(w1cycle_30);
ux_min30=mean(energy2_30);

y30=[uy_max30 uy_min30 ];
x30=[ux_max30 ux_min30 ];

m30=(uy_max30-uy_min30)/(ux_max30-ux_min30);
x_30middle=(ux_max30+ux_min30)/2;
y_30middle=(uy_max30+uy_min30)/2;

xq1_30=0:1:10;
yq1_30=y_30middle+m30*(xq1_30-x_30middle);

m30new=-(1/m30);
x30new= (0:0.0001:40);
y30new=y_30middle+m30new*(x30new-x_30middle);
%% classification 20%
uy_max20=sum(wcycle_20)/32;
ux_max20=sum((energy1_20)/32)-10; 
uy_min20=mean(w1cycle_20);
ux_min20=mean(energy2_20);

y20=[uy_max20 uy_min20 ];
x20=[ux_max20 ux_min20 ];

m20=(uy_max20-uy_min20)/(ux_max20-ux_min20);
x_20middle=(ux_max20+ux_min20)/2;
y_20middle=(uy_max20+uy_min20)/2;

xq1_20=0:40;
yq1_20=y_20middle+m20*(xq1_20-x_20middle);

m20new=-(1/m20);
x20new= 0:0.1:80;
y20new=y_20middle+m20new*(x20new-x_20middle);

%% monitor 50%threshold
figure
subplot(2,2,1)
scatter(energy1_50,wcycle_50,'filled')
xlabel('energy')
ylabel('index')
title('50%')
hold on
scatter(energy2_50,w1cycle_50,'r','filled')
hold on
scatter(x50,y50,'g')
plot(xq1_50,yq1_50)
plot(x50new,y50new)
legend('max','min')
%% monitor 40%threshold
subplot(2,2,2)
scatter(energy1_40,wcycle_40,'filled')
xlabel('energy')
ylabel('index')
title('40%')
hold on
scatter(energy2_40,w1cycle_40,'r','filled')
scatter(x40,y40,'g')
plot(xq1_40,yq1_40)
hold on
plot(x40new,y40new)
%% monitor 30%threshold
subplot(2,2,3)
scatter(energy1_30,wcycle_30,'filled')
xlabel('energy')
ylabel('index')
title('30%')
hold on
scatter(energy2_30,w1cycle_30,'r','filled')
scatter(x30,y30,'g')
plot(xq1_30,yq1_30)
hold on
plot(x30new,y30new)
%% monitor 20%threshold
subplot(2,2,4)
scatter(energy1_20,wcycle_20,'filled')
xlabel('energy')
ylabel('index')
title('20%')
hold on
scatter(energy2_20,w1cycle_20,'r','filled')
scatter(x20,y20,'g')
plot(xq1_20,yq1_20)
hold on
plot(x20new,y20new)
%% sensivity_50 & specificity_50
TP_50=0;
FP_50=0;
TN_50=0;
FN_50=0;
for k=1:32 
ytest_50=y_50middle+m50new*(energy1_50(k)-x_50middle);
if (wcycle_50(k)- ytest_50 )>0   %% max class  blue 
   TP_50=TP_50+1;  
else
   FP_50=FP_50+1;
end
end
for k=1:30
y1test_50=y_50middle+m50new*(energy2_50(k)-x_50middle);
if (w1cycle_50(k)-y1test_50 )<0
    TN_50=TN_50+1;
else
    FN_50=FN_50+1;
end
end
sensivity_50=TP_50*100/(TP_50+FP_50);
specifity_50=TN_50*100/(TN_50+FP_50);

%% sensivity_40 specificity_40
TP_40=0;
FP_40=0;
TN_40=0;
FN_40=0;
for k=1:32 
ytest_40=y_40middle+m40new*(energy1_40(k)-x_40middle);
if (wcycle_40(k)- ytest_40 )>0   %% max class  blue 
   TP_40=TP_40+1;  
else
   FP_40=FP_40+1;
end
end
for k=1:30
y1test_40=y_40middle+m40new*(energy2_40(k)-x_40middle);
if (w1cycle_40(k)-y1test_40 )<0
    TN_40=TN_40+1;
else
    FN_40=FN_40+1;
end
end
sensivity_40=TP_40*100/(TP_40+FP_40);
specifity_40=TN_40*100/(TN_40+FP_40);
%% sensivity_30 & specificity_30
TP_30=0;
FP_30=0;
TN_30=0;
FN_30=0;
for k=1:32 
ytest_30=y_30middle+m30new*(energy1_30(k)-x_30middle);
if (wcycle_30(k)- ytest_30 )>0   %% max class  blue 
   TP_30=TP_30+1;  
else
   FP_30=FP_30+1;
end
end
for k=1:30
y1test_30=y_30middle+m30new*(energy2_30(k)-x_30middle);
if (w1cycle_30(k)-y1test_30 )<0
    TN_30=TN_30+1;
else
    FN_30=FN_30+1;
end
end
sensivity_30=TP_30*100/(TP_30+FP_30);
specifity_30=TN_30*100/(TN_30+FP_30);
%%  sensivity_20 specificity_20
TP_20=0;
FP_20=0;
TN_20=0;
FN_20=0;
for k=1:32 
ytest_20=y_20middle+m20new*(energy1_20(k)-x_20middle);
if (wcycle_20(k)- ytest_20 )>0   %% max class  blue 
   TP_20=TP_20+1;  
else
   FP_20=FP_20+1;
end
end
for k=1:30
y1test_20=y_20middle+m20new*(energy2_20(k)-x_20middle);
if (w1cycle_20(k)-y1test_20 )<0
    TN_20=TN_20+1;
else
    FN_20=FN_20+1;
end
end
sensivity_20=TP_20*100/(TP_20+FP_20);
specifity_20=TN_20*100/(TN_20+FP_20);



