clc; 
out = pmdDataSetup(Data);
angles = out(:,1);
ang = [0 45 90 135 180 225 270 315];
firingRates = out(:,2:end);
meanRate =mean(firingRates,1);
numTrials = length(out);
dim = 3;
Ldim= dim;
Ndim = 232;
W = ones(Ndim,Ldim);
crisscross = diag(diag(ones(Ndim)));
zeroedObservations=zeros(size(firingRates));
for i=1:size(firingRates, 1)
   zeroedObservations(i, :)=firingRates(i,:)-meanRate; 
end
zeroedObservations = zeroedObservations';
for i =1:500
    C = W*W'+crisscross;
    Z1 = W'*C^(-1)*zeroedObservations;
    covar = eye(Ldim)-W'*(inv(C))*W;
    W2 = zeroedObservations*Z1'*(numTrials*covar+Z1*Z1')^(-1);
    crisscross2 = (1/numTrials)/Ndim*(diag(diag(zeroedObservations*zeroedObservations'-W2*Z1*zeroedObservations')));
    W = W2;
    crisscross = crisscross2;
    if i == 500
        C1 =W*W'+crisscross*eye(Ndim);
        embed=W'*C^(-1)*zeroedObservations;
    end
end
embeddedAngles = [angles, embed'];
colors = jet(length(ang));

for i = 1:length(ang)
   testAng = ang(i);
   index = embeddedAngles(:,1)==testAng;
   embedAng = embeddedAngles(index,:);
   if i==1
       figure;
        plot3(embedAng(:,2),embedAng(:,3),embedAng(:,4), 'Color', colors(i,:),'Marker', '.', ...
                        'linestyle', 'none')
   else 
       hold on
       plot3(embedAng(:,2),embedAng(:,3),embedAng(:,4), 'Color', colors(i,:),'Marker', '.', ...
                        'linestyle', 'none')
   end
end
comp = svd(W)