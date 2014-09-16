%%	draw the allan variance plot and analyze noise of inertial sensor output

 %% Import data

GRAVITY = 9.80665;
samplePeriod = 10/1000;

originalData = load('Sensor.txt');	

dataNum = length(originalData);

acc = zeros(dataNum,3);
gyro = zeros(dataNum,3);
accTime = zeros(dataNum,1);
gyroTime = zeros(dataNum,1);
logTime = zeros(dataNum,1);

acc = originalData(:,2:4);
gyro = originalData(:,6:8);

accTime = (originalData(:,5) - originalData(1,5)) / 1000000  + 500;
gyroTime = (originalData(:,9) - originalData(1,9)) / 1000000 + 500;
logTime = originalData(:,10) / 1000;

angularRate = zeros(size(gyro));
velocity = zeros(size(acc));

for i = 2 : length(angularRate)
	angularRate(i,:) = angularRate(i-1,:) + (gyro(i-1,:) + gyro(i,:)) * samplePeriod / 2;
end

for i = 2 : length(velocity)
	velocity(i,:) = velocity(i-1,:) + (acc(i-1,:) + acc(i,:)) * samplePeriod / 2;
end

%--------------------------------------------

maxJ = floor( log2(dataNum) - 1 ) ;
J = 0 : maxJ;
n = zeros(size(J));

for i = 1 : length(J)
	n(i) = power(2,J(i));
end

% t0 = 1 / samplePeriod;

for i = 1 : length(J)

	disp(['Computing i=',num2str(i)]);

	sumTheta = 0;

	for k = 1 : dataNum - 2 * n(i)

		sumTheta = sumTheta + (angularRate(k+2*n(i),1) - 2 * angularRate(k+n(i),1) + angularRate(k,1))^2;

	end

	sigma(i) = sumTheta / (2 * (n(i) * samplePeriod) * (n(i) * samplePeriod) * (dataNum - 2 * n(i)));

	sigma(i) = sqrt(sigma(i));


end

grid on;
loglog(n*samplePeriod,sigma);

