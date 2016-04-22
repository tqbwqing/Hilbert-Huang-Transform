%% example 1 -  2nd order Stokes wave
duration = 512; dt = 0.1; fs = 1/dt;
t = 0:dt:duration-dt;

a = 1;
k = 0.2;
f = 1/8;
omega = 2 * pi * f;
x = (1 / 2) * a^2 * k + a * cos(omega * t) + (1 / 2) * a^2 * k * cos(2 * omega * t);

imfs = eemd(x, 0, 1);
imfs = imfs(:, 2:end);

frequency_bin = length(t);
[HS, frequency_axis, time_axis] = HilbertSpectrum(imfs, frequency_bins);
frequency_axis = frequency_axis * fs / 2;
time_axis = time_axis * dt;

HSdB = 10 * log10(HS.^2 + 0.0001);  % Hilbert spectrum in dB
h = sum(HS, 2);                     % marginal spectrum
meanh = sum(HS, 2) / duration;      % mean marginal spectrum
IE = sum(HS.^2, 1);                 % instantaneous energy density level

DS = zeros(length(frequency_axis), 1);  % Degree of stationarity
for w = 1:1:length(frequency_axis)
    DS(w) = sum(1 - (HS(w, :) / meanh(w)), 2);
end
DS = DS / duration;

figure(1); 
plot(t, x)

figure(2);
thresh = sum(x.^2);
A = sum(imfs.^2, 1) > 0.1 * thresh ;
idx = find(A > 0);
for i = 1:1:length(idx)
    subplot(length(idx)+1, 1, i); plot(t, imfs(:, i));
end
subplot(length(idx)+1, 1, length(idx)+1); plot(t, sum(imfs(:, idx+1:end), 2))

figure(3);
imagesc('XData', time_axis, 'YData', frequency_axis, 'CData', HS);
axis([time_axis(1) time_axis(end) frequency_axis(1) frequency_axis(end)]);

figure(4);
imagesc('XData', time_axis, 'YData', frequency_axis, 'CData', HSdB);
axis([time_axis(1) time_axis(end) frequency_axis(1) frequency_axis(end)]);

figure(5); hold on;
X = fft(x); 
absX = abs(X)/length(X);
absXdB = 10 * log10(absX.^2 + eps);
F = (fs / 2) * linspace(0, 1, length(absX)/2);
plot(frequency_axis, 20*log10(h+eps)); plot(F, absXdB(1:end/2), 'r');

figure(6);
plot(time_axis, IE);

figure(7); 
plot(frequency_axis, DS);






