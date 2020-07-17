
%

[spikes, pos, vel] = Problem4_0_1(Data);

emptyVals = find(sum(vertcat(spikes{:})) == 0);

for n = 1:1222
    spikes{n}(:,36) = [];
    spikes{n}(:,37) = [];
end 

% firingRateMean=mean(spikes,1);
% firingRateCovariance=cov(spikes);
% 
% II = firingRateMean; Mx1
% 
% V = firingRateCovariance; MxM
% stateModel =

% Zt = A*Zt(i) + Noise # i refers to the previous trial
% input Zt is a velocity value that will be used to get an approximate
% postion value

% STEPS
% 1) Preprocessing -> initiating the A, C, Q, R variables with appropriate
% size MxM, use pmdDataSetup, make sure that there are 50 of each 
% 2) Initializing first time stamp, 

%% SORTING DATA BY REACH ANGLE
data = pmdDataSetup(Data);

angles = [0 45 90 135 180 225 270 315];
reachAngles = data(:,1);
count = 0;
trainingData.veloX = {808};
trainingData.veloY = {808};
trainingData.posX = {808};
trainingData.posY = {808};
trainingData.spikes= {808};

sortedIndex = [];
for i =1:length(angles)
    index = reachAngles==angles(i);
    for j=1:length(index)
        if index(j)==1&&count<101
            sortedIndex = [sortedIndex j];
            count = count +1;
        elseif index(j)==1&&count>=101
            count = 0;
            break;
        end
    end
end


%%
for i = 1:length(sortedIndex)
    trainingData.veloX{i} = vel.x{sortedIndex(i)};
    vel.x{sortedIndex(i)}(:) = [];
    trainingData.veloY{i} = vel.y{sortedIndex(i)};
    vel.y{sortedIndex(i)}(:) = [];
    trainingData.posX{i} = pos.x{sortedIndex(i)};
    pos.x{sortedIndex(i)}(:) = [];
    trainingData.posY{i} = pos.y{sortedIndex(i)};
    pos.y{sortedIndex(i)}(:) = [];
    trainingData.spikes{i} = spikes{sortedIndex(i)};
    spikes{sortedIndex(i)} = [];
end



%%

velX = vel.x{1}(1);
velY = vel.y{1}(1);
posX = pos.x{1}(1);
posY = pos.y{1}(1);

% 4 x 1
Zt = [velX, posX, posY, velY]'; % would it be a tuple of the x and y value? Or multidirectional movements is 4 =  Mx1


velX = vel.x{1}(2);
velY = vel.y{1}(2);
posX = pos.x{1}(2);
posY = pos.y{1}(2);

% 4 x 1
Zt2 = [velX, posX, posY, velY]'; % second Zt value for intialization

% [ vel.x ] 
% [ pos.x ] 
% [ pos.y ]    
% [ vel.y ] 


% Training Phase


% Initialization


  % sum from 2 to T (T being the number of bins)
  
  T = length(vel.x{1}); % num of bins used
  
  
  %% A UPDATE
  
  firstASum = zeros(4,4);
  secondASum = zeros(4,4);
  
   for trial=1:size(trainingData.veloX) % outer loop goes through all the trails 
       T=length(trainingData.veloX{trial}); % inner loop iterates through bin sizes of spike counts per trial
        
        % redefine Zt for each time bin 

        for binPos=2:T 
            
          if binPos == T + 1
       
              break;
          end 
          

        if(isnan(trainingData.veloX{trial}(binPos)) || isnan(trainingData.veloY{trial}(binPos)) || isnan(trainingData.posX{trial}(binPos)) || isnan(trainingData.posY{trial}(binPos)))
            
            continue;
            
        end 
        
         if(isnan(trainingData.veloX{trial}(binPos - 1)) || isnan(trainingData.veloY{trial}(binPos - 1)) || isnan(trainingData.posX{trial}(binPos - 1)) || isnan(trainingData.posY{trial}(binPos-1)))
            
            continue;
            
        end 
        
       velX1 = trainingData.veloX{trial}(binPos);
       velY1 = trainingData.veloY{trial}(binPos);
       posX1 = trainingData.posX{trial}(binPos);
       posY1 = trainingData.posY{trial}(binPos);
        
        % check isnan for vel/pos values, if so then continue 
           
       velX2 = trainingData.veloX{trial}(binPos - 1);
       velY2 = trainingData.veloY{trial}(binPos - 1);
       posX2 = trainingData.posX{trial}(binPos - 1);
       posY2 = trainingData.posY{trial}(binPos - 1);
        
        
        Zt = [posX1, posY1, velX1, velY1]'; %changed order to pos, vel
        
        Zt2 = [posX2, posY2, velX2, velY2]';
                
        firstASum =  Zt2*Zt' + firstASum;

        secondASum = Zt*Zt' + secondASum;
        
        end                
        
   end 
   %%
   
   %COMPUTE A
   
   A = firstASum * secondASum^(-1); % is MxM so 4x4, computing from only Z1 and Z2
 
   %%
  
   %% Q UPDATE
   
   QSum = zeros(4,4);
   firstTermQ = 0;
   
    for trial=1:size(trainingData.veloX) % outer loop goes through all the trails 
         
       T=length(trainingData.veloX{trial}); % inner loop iterates through bin sizes of spike counts per trial
        
        % redefine Zt for each time bin 

        for binPos=2:T 
            
             if binPos == T + 1
       
              break;
             end 
             
         if(isnan(trainingData.veloX{trial}(binPos)) || isnan(trainingData.veloY{trial}(binPos)) || isnan(trainingData.posX{trial}(binPos)) || isnan(trainingData.posY{trial}(binPos)))
            
            continue;
            
        end 
        
         if(isnan(trainingData.veloX{trial}(binPos - 1)) || isnan(trainingData.veloY{trial}(binPos - 1)) || isnan(trainingData.posX{trial}(binPos - 1)) || isnan(trainingData.posY{trial}(binPos-1)))
            
            continue;
            
        end 
          
       velX1 = trainingData.veloX{trial}(binPos);
       velY1 = trainingData.veloY{trial}(binPos);
       posX1 = trainingData.posX{trial}(binPos);
       posY1 = trainingData.posY{trial}(binPos);
        
       velX2 = trainingData.veloX{trial}(binPos - 1);
       velY2 = trainingData.veloY{trial}(binPos - 1);
       posX2 = trainingData.posX{trial}(binPos - 1);
       posY2 = trainingData.posY{trial}(binPos - 1);
        
        
        
        Zt = [posX1, posY1, velX1, velY1]'; %changed order to pos, vel
        
        Zt2 = [posX2, posY2, velX2, velY2]';
            
        
        QSum = (Zt2-A*Zt)*(Zt2-A*Zt)' + QSum;


        end
        
        firstTermQ = 1/((T-1)+ firstTermQ);

        
    end 
      
   %% COMPUTE Q
   
  Q = firstTermQ * QSum; % M x M
  

%% C Update Loop

firstSumC = zeros(231,4);
secondSumC = zeros(4,1);
Xt = zeros(231,1);

for trial=1:size(trainingData.veloX) % outer loop goes through all the trails 
         
       T=length(trainingData.veloX{trial}); % inner loop iterates through bin sizes of spike counts per trial
        
       % COMPUTE Xt
 
  % gets the first column for spike counts, 233x1
       
         Xt = trainingData.spikes{trial}(2,:)'; % D = numNeurons needs to be 233 x 1

       
        
        % redefine Zt for each time bin 

        for binPos=1:T 
            
            
          if(isnan(trainingData.veloX{trial}(binPos)) || isnan(trainingData.veloY{trial}(binPos)) || isnan(trainingData.posX{trial}(binPos)) || isnan(trainingData.posY{trial}(binPos)))
            
            continue;
            
          end 
           
       velX1 = trainingData.veloX{trial}(binPos);
       velY1 = trainingData.veloY{trial}(binPos);
       posX1 = trainingData.posX{trial}(binPos);
       posY1 = trainingData.posY{trial}(binPos);
        
         Zt = [posX1, posY1, velX1, velY1]'; %changed order to pos, vel
                         
        firstSumC = (Xt*Zt') + firstSumC;
        
        secondSumC = (Zt*Zt')^(-1) + secondSumC; 
        
        end
        
end
%% COMPUTE C

        C = firstSumC * secondSumC; % D x M multiplying a 233x1 by a 1x4 results in 233 x 4


%% R Update Loop

 RSum = zeros(231,231);
 Xt = zeros(231,1);
 firstTermR = 0;
   
    for trial=1:size(trainingData.veloX) % outer loop goes through all the trails 
         
       T=length(trainingData.veloX{trial}); % inner loop iterates through bin sizes of spike counts per trial
        
       Xt = trainingData.spikes{trial}(2,:)'; % D = numNeurons needs to be 233 x 1

        % redefine Zt for each time bin 

        for binPos=1:T 
             
             if(isnan(trainingData.veloX{trial}(binPos)) || isnan(trainingData.veloY{trial}(binPos)) || isnan(trainingData.posX{trial}(binPos)) || isnan(trainingData.posY{trial}(binPos)))
            
            continue;
            
             end 
              
       velX1 = trainingData.veloX{trial}(binPos);
       velY1 = trainingData.veloY{trial}(binPos);
       posX1 = trainingData.posX{trial}(binPos);
       posY1 = trainingData.posY{trial}(binPos);
                
        Zt = [posX1, posY1, velX1, velY1]'; %changed order to pos, vel
        
        RSum = (Xt - C*Zt)*(Xt-C*Zt)'; + RSum;


        end
        
        firstTermR = 1/(T+ firstTermR);

        
    end 

%% COMPUTE R

R = firstTermR * RSum; % D x D (233 x 233)

  

%% Testing Phase

% predicting next xPos based on the kalfilter

PI = 0;
V = 0;

%% PREPROCESSING THE TESTING DATA, removing empty cells, put in new structs

for i=1:size(pos.x,2)
    if isempty(pos.x{i})
        pos.x{i} = 0;
    end 
     if isempty(pos.y{i})
        pos.y{i} = 0;
     end  
     if isempty(vel.x{i})
        vel.x{i} = 0;
     end 
     if isempty(vel.y{i})
        vel.y{i} = 0;
     end 
     if isempty(spikes{i})
        spikes{i} = 0;
    end 
end

posXRemovedZeros = pos.x(cellfun(@(x) ~isequal(x, 0), pos.x));
posYRemovedZeros = pos.y(cellfun(@(x) ~isequal(x, 0), pos.y));
velXRemovedZeros = vel.x(cellfun(@(x) ~isequal(x, 0), vel.x));
velYRemovedZeros = vel.y(cellfun(@(x) ~isequal(x, 0), vel.y));
spikesRemovedZeros = spikes(cellfun(@(x) ~isequal(x, 0), spikes));



%% LOOPING KALMAN FUNCTION

predictions.predictedmu = {808};
predictions.predictedsigma = {808};

for trial=1:size(pos.x) % outer loop goes through all the trails 
         
       T=length(posXRemovedZeros{trial}); % inner loop iterates through bin sizes of spike counts per trial
        
       Xt = spikesRemovedZeros{trial}'; % D = numNeurons needs to be 231 x 1
                                              
       xPos = posXRemovedZeros{trial};

       [mu,sig] = kalmanFilter(xPos, A, Q, R, C, Xt, V);

       predictions.predictedmu{trial} = mu(1,1:T);
%        predictions.predictedsigma{trial} = sig(;
       
        
end 

function [mu, sig] = kalmanFilter(xPos, A, Q, R, C, Xt, V)


T = size(Xt,2); % needs to be the bin size
% 
mu = zeros(1,T); % 2 by bin size
sig = zeros(1,T);

for n=1:T
    if n==1
        mu1(:,n) = xPos(n);
        sig1{n} = V;
    else 
        mu1(:,n) = A*mu(:,n-1);
        sig1{n} = A*sig{n-1}*A' + Q;
    end 
    
    Kt = sig1{n}*C'*inv(C*sig1{n}*C'+R);
    mu(:,n) = mu1(:,n) + Kt*(Xt(:,n) - C*mu1(:,n));
    sig{n} = sig1{n} - Kt*C*sig1{n}; 
    
end 


end 





%%










