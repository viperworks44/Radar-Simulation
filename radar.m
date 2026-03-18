%4o set Radar Stampoulidis Georgios 03094
PRI = 1e-3; %Pulse Repetition Interval
DC = 0.1; %Duty Cycle
fs = 2e4; %Sampling Frequency
K = 101; %Number of Samples of FFT
c = 3e8; %Speed of Light
Pt = 1; %Transmit Power
G = 1; %Antenna Gain
Lo = 1; %Loss Factor
sigma = 1; %RCS
lamda = 1; %Wave Length
W = 0; %Noise
R0 = 5; %Distance
fc= 1e9; %Frequence of x
M = 20; %Number of Pulses
v = 0; %Velocity
SNR = 5; %SNR


t = 0:1/fs:PRI-1/fs;
pulse_width = DC*PRI;
x = rectpuls(t - pulse_width/2,pulse_width);
k = -K/2:(K/2-1);
PULSE = fft(x,K);
figure('Name','DFT With FFT help','NumberTitle','off');
stem(k/K,fftshift(abs(PULSE)));
title('Magnitude Spectrum');
xlabel('Signal Frequency/Sampling Frequency (f/fs)','FontSize',14);
ylabel('DFT Amplitude','FontSize',14);

%Range bin
DR = c*pulse_width/2 %range resolution
Rbin = R0/DR; %bin resolution
N = PRI*fs; %number of range bins
Ntarget = ceil(fs*2*R0/c);

fd = 2*v*fc/c;
la = c/fc;
Dfd = 1/(M*PRI);
DV = c*Dfd/(2*fc);
maxfd= 1/(2*PRI);
maxV = maxfd*c/(2*fc);
maxR = c*(PRI - DC*PRI)/2;

Pr = (Pt*G^2*lamda^2*sigma)/((4*pi)^3*R0^4*Lo);

y2d = zeros(N, M);
for m = 1:M
    for n = 1:N
        y2d(n,m) = sqrt(Pr)*rectpuls(((n-1)/fs)-(2*R0/c)-(pulse_width/2),pulse_width)*exp(-1i*4*pi*fc*R0/la)*exp(1i*2*pi*(v/la)*(((n-1)/fs)+(m-1)*PRI)) + sqrt(Pr/(10.^(SNR*0.6)*randn(1,1))); 
    end
end

Fbin = fd/Dfd
Dbin = round(Fbin)

for n = 1:N
    RDP(n,:) = fftshift(abs(fft(y2d(n,:),M)));
end

range = linspace(0,maxR,N-1);
velocity = linspace(-maxV+DV/2,maxV-DV/2,M+1);
figure('Name','Range-Doppler Profile (RDP)','NumberTitle','off');
imagesc(velocity,range,RDP);
title('Range-Doppler Profile');
xlabel('Velocity','FontSize',14);
ylabel('Range','FontSize',14);  
set(gca,'YDir','normal');

%---------------------#set4#---------------------
d = 0.25;
L = 10;
SNR = 5;
fi = pi/3;
y3d = zeros(N, M, L);

for l = 1:L
    for m = 1:M
        for n = 1:N
            y3d(n,m,l) = sqrt(Pr)*rectpuls(((n-1)/fs)-(2*R0/c)-(pulse_width/2),pulse_width)*exp(-1i*4*pi*fc*R0/la)*exp(1i*2*pi*(v/la)*(((n-1)/fs)+(m-1)*PRI))*exp(-1i*2*pi*d*(l-1)*cos(fi)/la) + sqrt(Pr/(10.^(SNR)*randn(1,1))) ;
        end
    end
end

[n_grid, m_grid, l_grid] = meshgrid(1:N, 1:M, 1:L);
figure('Name','3D Datacube','NumberTitle','off');
slice(n_grid,m_grid,l_grid, abs(y3d), 1:N, 1:M, 1:L);
xlabel('n');
ylabel('m');
zlabel('l');
title('3D imagesc Plot');
colorbar;

for theta1 = 0:179
    MagAz(theta1+1) = sin(pi*d*L*cos(deg2rad(theta1))/la)/sin(pi*d*cos(deg2rad(theta1))/la);
end

DbAz =mag2db(abs(MagAz)/max(abs(MagAz)));
figure('Name','DB vs Angle','NumberTitle','off');
plot(DbAz);
title('Azimuth Cut (elevation angle = 0.0)');
xlabel('Azimuth Angle (degrees)');
ylabel('Normalised Power (dB)');


for l=1:L
     ADP(l)=fftshift(abs(fft(y3d(1,1,l))));
end

figure('Name','Angle-Doppler Profile (RDP)','NumberTitle','off');
imagesc(ADP)
title('Angle-Doppler Profile');
xlabel('angle','FontSize',14);
ylabel('velocity','FontSize',14);  
set(gca,'YDir','normal');