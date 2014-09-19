%%	draw the allan variance plot and analyze noise of inertial sensor output

 %% Import data

GRAVITY = 9.80665;
samplePeriod = 10/1000;

dataNum = 100000;

acc = zeros(dataNum,3);
gyro = zeros(dataNum,3);
accTime = zeros(dataNum,1);
gyroTime = zeros(dataNum,1);
logTime = zeros(dataNum,1);

acc(:,1) = pinknoise(dataNum);
acc(:,2) = pinknoise(dataNum);
acc(:,3) = pinknoise(dataNum);

% gyro = originalData(:,6:8);

% accTime = (originalData(:,5) - originalData(1,5)) / 1000000  + 500;
% gyroTime = (originalData(:,9) - originalData(1,9)) / 1000000 + 500;
% logTime = originalData(:,10) / 1000;

% angularRate = zeros(size(gyro));
velocity = zeros(size(acc));
position = zeros(size(velocity));

% for i = 2 : length(angularRate)
% 	angularRate(i,:) = angularRate(i-1,:) + (gyro(i-1,:) + gyro(i,:)) * samplePeriod / 2;
% end

for i = 2 : length(velocity)
	velocity(i,:) = velocity(i-1,:) + (acc(i-1,:) + acc(i,:)) * samplePeriod / 2;
end

for i = 2 : length(velocity)
	position(i,:) = position(i-1,:) + (velocity(i-1,:) + velocity(i,:)) * samplePeriod / 2;	
end

%--------------------------------------------

maxJ = floor( log2(dataNum) - 1 ) ;
J = 0 : maxJ;
n = zeros(size(J));

for i = 1 : length(J)
	n(i) = power(2,J(i));
end

% t0 = 1 / samplePeriod;

% % Allan Variance plot for Gyro

% for i = 1 : length(J)

% 	sumTheta = [0; 0; 0];

% 	for gyroAxis = 1 : 3

% 		for k = 1 : dataNum - 2 * n(i)

% 			sumTheta(gyroAxis,1) = sumTheta(gyroAxis,1) + (angularRate(k+2*n(i),gyroAxis) - 2 * angularRate(k+n(i),gyroAxis) + angularRate(k,gyroAxis))^2;

% 		end

% 		sigmaGyro(i,gyroAxis) = sumTheta(gyroAxis,1) / (2 * (n(i) * samplePeriod) * (n(i) * samplePeriod) * (dataNum - 2 * n(i)));
% 		sigmaGyro(i,gyroAxis) = sqrt(sigmaGyro(i,gyroAxis));

% 	end
	
% end

% figure('Number','off','Name','Allan Deviation plot for Gyroscope')
% hold on;
% set(gca, 'XScale', 'log', 'YScale', 'log');
% grid on;
% loglog(n*samplePeriod,sigmaGyro(:,1), 'r-');
% loglog(n*samplePeriod,sigmaGyro(:,2), 'g-');
% loglog(n*samplePeriod,sigmaGyro(:,3), 'b-');
% title('Allan Variance plot for Gyroscope');
% legend('Gyro-X','Gyro-Y','Gyro-Z');


% Allan Variance plot for Acc

for i = 1 : length(J)

	sumTheta = [0; 0; 0];

	for accAxis = 1 : 3

		for k = 1 : dataNum - 2 * n(i)

			sumTheta(accAxis,1) = sumTheta(accAxis,1) + (velocity(k+2*n(i),accAxis) - 2 * velocity(k+n(i),accAxis) + velocity(k,accAxis))^2;

		end

		sigmaAcc(i,accAxis) = sumTheta(accAxis,1) / (2 * (n(i) * samplePeriod) * (n(i) * samplePeriod) * (dataNum - 2 * n(i)));
		sigmaAcc(i,accAxis) = sqrt(sigmaAcc(i,accAxis));

	end
	
end

figure('Number','off','Name','Allan Deviation plot for Accelerometer')
hold on;
set(gca, 'XScale', 'log', 'YScale', 'log');
grid on;
loglog(n*samplePeriod,sigmaAcc(:,1), 'r-');
loglog(n*samplePeriod,sigmaAcc(:,2), 'g-');
loglog(n*samplePeriod,sigmaAcc(:,3), 'b-');
title('Allan Variance plot for Accelerometer');
legend('Acc-X','Acc-Y','Acc-Z');


