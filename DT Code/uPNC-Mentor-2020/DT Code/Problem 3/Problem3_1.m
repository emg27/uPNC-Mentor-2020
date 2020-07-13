clc; 
out = pmdDataSetup(Data);
angles = out(:,1);
ang = [0 45 90 135 180 225 270 315];
firingRates = out(:,2:end);
meanRate =mean(firingRates,1);
firingRatemeans = [];
covarRates = cov(firingRates);
[PC,Eigen]=eig(covarRates);
eigenVal = fliplr(sum(Eigen));
PC = fliplr(PC);
varTot = sum(eigenVal);
embed = [];
figure;
run=1:length(eigenVal);
plot(run,eigenVal)
xlabel('Eigenvalue Rank')
ylabel('Variance')
check = [3 10];
varEx = [];
%%
for i =check
   varEx = sum(eigenVal(1:i))/varTot;
end
for trial=1:size(firingRates,1)
    firingRatemeans = [firingRatemeans;meanRate];
end

for i=1:3
   scores = (firingRates-firingRatemeans)*PC(:,i);
   embed = [embed,scores];
end
embeddedAngles = [angles, embed];
colors = jet(length(ang));

for i = 1:length(ang)
   testAng = ang(i);
   index = embeddedAngles(:,1)==testAng;
   embedAng = embeddedAngles(index,:);
   if i==1
        plot3(embedAng(:,2),embedAng(:,3),embedAng(:,4), 'Color', colors(i,:),'Marker', '.', 'linestyle', 'none')
   else 
       hold on
       plot3(embedAng(:,2),embedAng(:,3),embedAng(:,4), 'Color', colors(i,:),'Marker', '.', 'linestyle', 'none')
   end
end