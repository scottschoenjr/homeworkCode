
%Sep 13 Lecture Guide
%Brendan Nichols

%% NARROWBAND BEAMFORMING EXAMPLE
%repeated from Sep 8 lecture

%define an example line-array of sensors along the X-axis
N = 10; %number of sensors
d = 15; %spacing between sensors (m)
pos = [(0:N-1)*d; zeros(1, N); zeros(1, N)]; %positions of each sensor (3xN)

%generate some fake data collected by the array
%don't worry about this step, I am making up the magnitudes/phases of the 50 Hz wave that each sensor records
%normally, this information is obtained by taking actual time measurements on an array, and taking an FFT
f = 50; %frequency that was recorded (Hz)
c = 1500; %sound speed of recording medium (m/s)
data = exp(-1i*0.1047*(0:N-1)*d); %the complex amplitudes of the 50Hz components at each sensor

%% create a steering/weight/replica vector to look in a given direction
theta = 0; %look angle (degrees)
l = [cos(theta*pi/180); sin(theta*pi/180); 0]; %look direction (unit vector)
tau = sum(bsxfun(@times, l, pos), 1)/c; %time delays which would be experienced by each sensor if this is the true theta (sec)
w = exp(1i*2*pi*f*tau); %weight vector

%% beamform the data
b = w*data.';

%% do the above steps in a loop, searching many possible angles
theta = linspace(0, 360, 721); %search angles (deg)
b = zeros(size(theta)); %initialize beamformer output at each theta

for i = 1:length(theta)
    %create weight vector
    l = [cos(theta(i)*pi/180); sin(theta(i)*pi/180); 0]; %look direction (unit vector)
    tau = sum(bsxfun(@times, l, pos), 1)/c; %time delays which would be experienced by each sensor if this is the true theta (sec)
    w = exp(1i*2*pi*f*tau); %weight vector
    
    %beamform the data
    b(i) = w*data.';
end

%% plot the magnitude of the delayed and summed beamformer output vs. search direction
plot(theta, abs(b))

%% WIDEBAND BEAMFORMING EXAMPLE
%{
Notes on dimensions of variables used below
Dimension 1 = Time or frequency
Dimension 2 = Sensor number or search space
Dimension 3 = Spatial dimension (X/Y/Z)
%}

%define an example line-array of sensors along the X-axis
N = 10; %number of sensors
d = 15; %spacing between sensors (m)
pos = permute([(0:N-1)*d; zeros(1, N); zeros(1, N)], [3 2 1]); %positions of each sensor (1xNx3) (m)

%% generate some fake data collected by the array
fs = 500; %sample frequency (Hz)
T = 10; %time-duration of the data (sec)
Nt = T*fs; %number of time points
t = (0:Nt-1).'/fs; %time axis (s)

%create a function which describes the source signal
f = @(t) sin(2*pi*50*t).*exp(-((t-5)/0.1).^2);

%% create a data vector which contains delayed versions of the source signal
%NOTE: this is only to make an example dataset, do not use this code for HW2
delays = (0.009*(0:N-1)); %true angle of arrival = 25.8 deg
tdelayed = bsxfun(@minus, t, delays);
data = f(tdelayed); %the data vector (Nt x N)

%% now that we have some data, take it into the frequency domain
freq = (0:Nt-1).'/Nt*fs; %frequency axis (Nt x 1)
Data = fft(data, [], 1); %take the FFT of the data along the time dimension (#1)

%set up all the look directions to scan
theta = linspace(0, 180, 181); %angles to scan (degrees)
Ntheta = length(theta); %number of look angles

%initialize the output of the beamformer for each look direction
b = zeros(Nt, Ntheta); %beamformer output (Nt x Ntheta)

%% beamform the data, select which method to use with a flag
BFNAIVE = 0; %for reference: profiler self time = 5.2  s (x1)
BFFAST = 1;  %for reference: profiler self time = 0.75 s (x7)
BFFASTER = 2;%for reference: profiler self time = 0.37 s (x14)
BFMETHOD = BFNAIVE;

%% NAIVE BEAMFORMER
if(BFMETHOD == BFNAIVE)
%for each look direction (dimension = 2)
for j = 1:Ntheta
    %generate the look direction and time delays
    l = permute([cos(theta(j)*pi/180); sin(theta(j)*pi/180); 0], [3 2 1]); %look direction (unit vector) (1x1x3)
    tau = sum(bsxfun(@times, l, pos), 3)/c; %time delays which would be experienced by each sensor if this is the true theta (sec)
    
    %for each frequency up to Nyquist (dimension = 1)
    for i = 1:Nt/2
        %create the weight vector for this look direction and frequency
        w = exp(1i*2*pi*freq(i)*tau); %weight vector
        
        %beamform this frequency in the data with the weight vector
        %store the result for the same frequency and this look direction
        b(i, j) = abs(w*Data(i, :).');
    end
end
end

%% FAST BEAMFORMER
if(BFMETHOD == BFFAST)
%only work on frequencies within the actual band of the signal
band = [40 60]; %bandwidth limits (Hz)

%for each look direction (dimension = 2)
for j = 1:Ntheta
    %generate the look direction and time delays
    l = permute([cos(theta(j)*pi/180); sin(theta(j)*pi/180); 0], [3 2 1]); %look direction (unit vector) (1x1x3)
    tau = sum(bsxfun(@times, l, pos), 3)/c; %time delays which would be experienced by each sensor if this is the true theta (sec)
    
    %for each frequency in the band (dimension = 1)
    for i = find((freq >= band(1) & freq <= band(2))).'; %indices which are in the band
        %create the weight vector for this look direction and frequency
        w = exp(1i*2*pi*freq(i)*tau); %weight vector
        
        %beamform this frequency in the data with the weight vector
        %store the result for the same frequency and this look direction
        b(i, j) = abs(w*Data(i, :).');
    end
end
end

%% FASTER BEAMFORMER
if(BFMETHOD == BFFASTER)
%eliminate for loops, work with matrices when possible
%AND only work on frequencies within the actual band of the signal
band = [40 60]; %bandwidth limits (Hz)
ii = (freq >= band(1) & freq <= band(2)); %indices which are in the band

%for each look direction (dimension = 2)
for j = 1:Ntheta
    %generate the look direction and time delays
    l = permute([cos(theta(j)*pi/180); sin(theta(j)*pi/180); 0], [3 2 1]); %look direction (unit vector) (1x1x3)
    tau = sum(bsxfun(@times, l, pos), 3)/c; %time delays which would be experienced by each sensor if this is the true theta (sec)
    
    %generate a weight vector for all frequencies in the band at once
    w = zeros(Nt, N); %initialize the weight vector with zeros
    w(ii, :) = exp(1i*2*pi*bsxfun(@times, freq(ii), tau)); %fill the weight vector for the frequencies in the band
    
    %multiply all the weight vectors by the data vectors at once
    b(:, j) = abs(sum(w.*Data, 2));
end
end

%% take the ifft of the beamformed outputs to get them back into the time domain
b = 2*real(ifft(b, [], 1))/N;

%% plot the results
imagesc(t, theta, b.')
set(gca, 'ydir', 'normal')
colorbar
xlim([4.5 5.5])
xlabel('Time (s)')
ylabel('Look Angle (\circ)')

%% plot more meaningful results (example: the maximum value of the beamformed output)
plot(theta, max(b))