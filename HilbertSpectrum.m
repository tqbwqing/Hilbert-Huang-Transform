function [HS, frequency_axis, time_axis] = HilbertSpectrum(input, fs, frequency_bins, frequency_range, to_plot)

% --------------------------------------------------------------------------------------------------
%   HilbertSpectrumCalc: Hilbert spectrum calculator.
%   args: 
%       - input: matrix with the set of inputs signals, every signal has to be in every row.
%       - fs: sampling frequency, in the case of time series the fs=1
%       - frequency_bins: quantization in frequency axis
%       - frequency_range: two values that define the frequency range used in Hz
%
%   returns:
%       - HS: 2D matrix which contains the Hilbert spectrum
%       - frequency_axis
%       - time_axis
% --------------------------------------------------------------------------------------------------

[number_of_signals, number_of_samples] = size(input);
if number_of_samples == 1      % if there is only one singal which is row vector
    input = input.';
    number_of_samples = number_of_signals;
    number_of_signals = 1;
end

if nargin == 1
    fs = 1;
    frequency_bins = number_of_samples;
    frequency_range = [0 fs/2];
    to_plot = 1;
elseif nargin == 2
    frequency_bins = number_of_samples;
    frequency_range = [0 fs/2];
    to_plot = 1;
elseif nargin == 3
    frequency_range = [0 fs/2];
    to_plot = 1;
elseif nargin == 4
    to_plot = 1;
end

frequency_axis = linspace(frequency_range(1), frequency_range(2), frequency_bins);
time_axis = 0:1:(number_of_samples - 1);
HS = zeros(length(frequency_axis), length(time_axis));

if frequency_range(2) > fs/2
    disp('frequency range invalid');
    return
end

for i = 1:1:number_of_signals
    
    xh      = hilbert(input(i, :));
    xh_amp  = abs(xh);
    xh_if   = [diff(unwrap(angle(xh))), 0] * fs / 2 / pi;
    
    HS_temp  = zeros(length(frequency_axis), length(time_axis));
    [~, bin] = histc(xh_if, frequency_axis);
    idx = bin == 0;
    bin(idx) = 1;
    
    for n = 1:1:length(time_axis)        
        HS_temp(bin(n), n) = xh_amp(n);
    end
    
    HS = HS + HS_temp;
    
end

time_axis = time_axis / fs;

if to_plot == 1
    figure;
    imagesc('XData', time_axis, 'YData', frequency_axis, 'CData', HS);
    axis([time_axis(1) time_axis(end) frequency_axis(1) frequency_axis(end)]);
    xlabel('Time (sec)');
    ylabel('Frequency (Hz)');
end

end

