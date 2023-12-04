
load GTA.data
plot(GTA(:,14))

e4init
pkg load control

[VAR,P,TVP,oar]=autodhr(GTA(:,14),1,[],[50:60],[Inf,60],[1,1;1,0],1)
[Trend,season,cycle,irreg]=dhrfilt(GTA(:,14),P,TVP,VAR,12,0,0);
figure(1)
plot([Trend(:,1:2),GTA(:,14),cycle+0.7] )


[VAR,P,TVP,oar]=autodhr(GTA(:,14),1,[],[50:70],[Inf,60],[1,.985;1,0],0)
[Trend,season,cycle,irreg]=dhrfilt(GTA(:,14),P,TVP,VAR,12,0,0);
figure(1)
plot([Trend(:,1:2),GTA(:,14),cycle+0.7] )
figure(2)
plot(irreg)

[VAR,P,TVP,oar]=autodhr(GTA(:,14),1,[],[50:70],[Inf,60,30],[1,1;1,0],1)
[Trend,season,cycle,irreg]=dhrfilt(GTA(:,14),P,TVP,VAR,12,0,0);
figure(1); plot([Trend(:,1:2),GTA(:,14),cycle+1,irreg-1,-ones(size(irreg))] )

[VAR,P,TVP,oar]=autodhr(GTA(:,14),1,[],[50:80],[Inf,55.5],[1,1;1,0],1)
[Trend,season,cycle,irreg]=dhrfilt(GTA(:,14),P,TVP,VAR,12,0,0);
figure(3); plot([Trend(:,1:2),GTA(:,14),cycle+0.7] )

  
[VARm,Pm,TVPm,oarm]=autodhr(K,12,[],[100:110],[Inf,720],[1,1;1,0],1)
[Trendm,seasonm,cyclem,irregm]=dhrfilt(K,Pm,TVPm,VARm,12,0,0);
figure(3); plot([Trendm(:,1:2),K,cyclem+0.7] )

[VAR,P,TVP,oar,results,MCNN,NVR]=autodhr(GTA(:,14),1,[],[55:83],[inf,60,20],[1, 1 ; 1, 0],0)
[Trend,season,cycle,irreg]=dhrfilt(GTA(:,14),P,TVP,VAR,12,2,0);
figure(10); plot([Trend(:),GTA(:,14),sum(cycle,2)+0.7,irreg-1] )
[Trend,season,cycle,irreg]=dhrfilt(GTA(:,14),P,TVP,VAR,12,0,0);
figure(10); plot([Trend(:,1:2),GTA(:,14),cycle+0.7,irreg-1] )

[VAR,P,TVP,oar,results,MCNN,NVR]=autodhr(GTA(:,14),1,[],[55:83],[inf,63,21,9],[1, 1 ; 1, 0],0)
