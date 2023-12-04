clear
e4init

pkg load control

load gta.data

y = gta(:,14);

%plot --format png
figure(1)
grid off
hold on
title('Global Annual Temperature Anomalies')
plot(y, 'r')

#[VAR,P,TVP,oar,results,MCNN,NVR]=autodhr(y,1,[],[50:70],[],[],0);
[VAR,P,TVP,oar]=autodhr(y,1,[],[50:70],[],[],0)

[Trend,season,cycle,irreg]=dhrfilt(y,P,TVP,VAR,12,0,0);

plot([Trend,y,irreg-1] )

[ARSPT, S, LAGS, AR, ROOTS, NORM, P] = aresp(Trend,oar);
[s,i]=sort(NORM,'descend');
[P,NORM](i(1:18),:)
63./[1:9]

#[VAR,P,TVP,oar,results,MCNN,NVR]=autodhr(y,1,[],[50:80],[inf,63,21,9],[1, 1 ; 1, 0],0)
[VAR,P,TVP,oar]=autodhr(y,1,[],[50:70],[inf,63,21,9])

[Trend,season,cycle,irreg]=dhrfilt(y,P,TVP,VAR,12,0,1);

plot([Trend] )

plot([Trend(:,2),y,cycle+0.7,irreg-1] )

plot([Trend(:,2)+cycle,y] )


