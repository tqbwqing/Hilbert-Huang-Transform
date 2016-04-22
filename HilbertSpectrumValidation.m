% --------------------------------------------------------------------------------------------------
% Hilbert Spectrum Validation 
% a few examples derived by Huang's paper
% --------------------------------------------------------------------------------------------------

%% example 1: 2nd Stokes wave in deep water
duration = 512;
dt = 0.01; fs = 1/dt;
t = 0:dt:duration-dt;
a = 1;
k = 0.2;
f = 1/32;
omega = 2 * pi * f;
x = (1 / 2) * a^2 * k +...
    a * cos(omega * t) + ...
    (1 / 2) * a^2 * k * cos(2 * omega * t);
figure; plot(t, x)

imfs = eemd2(x, 0, 1);
input = imfs(:, 2).';

HilbertSpectrumCalcAnalog(input, fs, length(t)/1, [0 0.1], 1);
%[HS, frequency_axis, time_axis] = HilbertSpectrumCalc(input, fs, 200);
%% example 2
duration = 512;
dt = 1; fs = 1/dt;
t = 0:dt:duration-dt;
x = exp(-0.01 * t) .* cos(2 * pi * t / 32);
figure; plot(t, x);

%HilbertSpectrumCalcAnalog(x.', fs, length(t)/4, [0.02 0.04], 1);

[HS, frequency_axis, time_axis] = HilbertSpectrum(x.', length(t)/1);

E = 20 * log ( HS.^2 + 0.00001);
figure;
imagesc('XData', time_axis, 'YData', frequency_axis, 'CData', E);
axis([time_axis(1) time_axis(end) frequency_axis(1) frequency_axis(end)]);
xlabel('Time (sec)');
ylabel('Frequency (Hz)');

[S,freq] = hspec(x.',100);
%% example 3
duration = 1000;
dt = 0.1; fs = 1/dt;
t = 0:dt:duration-dt;
t2 = 0:dt:duration/2-dt;
x1 = cos(2 * pi * 0.03 * t2);
x2 = cos(2 * pi * 0.01 * t2);
x = [x1, x2];
figure; plot(t,x);

%HilbertSpectrumCalcAnalog(x, fs, length(t)/32, [0.00 0.05], 1);

%% example 4
duration = 512;
dt = 0.01; fs = 1/dt;
t = 0:dt:duration-dt;
x = cos(2 * pi * t / 30) + cos(2 * pi * t / 34);
figure; plot(t, x);

imfs = eemd2(x, 0, 1);
N = 5;
figure; 
for i = 1:1:N
    subplot(N, 1, i); plot(imfs(:, i));
end
input = imfs(:, 1:2).';
HilbertSpectrum(x, fs, length(t)/4, [0.01 0.05], 1);
%% example 4 - intrawave modulation
duration = 512;
dt = 1; fs = 1/dt;
t = 0:dt:duration-dt;
e = 0.3;
omega = 2 * pi /64;
x = cos(omega * t + e * sin(2 * omega * t));
figure; plot(t, x);
imfs = eemd2(x, 0.1, 100);
N = 5;
figure; 
for i = 1:1:N
    subplot(N, 1, i); plot(imfs(:, i));
end


