function [components, pool_prev_omegas, swarm_status] = SwarmDecomposition(signal, param_struct)

% -------------------------------------------------------------------------
%   swarmDecomposition6(): swarm decomposition v6
%   input:
%    1) signal: the input signal
%    2) param_struct: method's parameters
%       param_struct = {min_peak_thresh, component_std, sgolay_degree, sgolay_length}
%   output:
%    1) components: output of method
%    2) recognized omega
%
%   All ad-hoc subroutines needed are included
% -------------------------------------------------------------------------

swarm_status = 1;   % flag of algorithm status

x = signal(:);

% parameters
pks_thresh                  = param_struct.min_peak_thresh;
component_termination_crit  = param_struct.component_std;
smooth_degree               = param_struct.welch_window;
smooth_length               = param_struct.welch_noverlap;

% recognition of oscillatory behaviour - find oscillRecon_ref_max_ampl (ref)
[~, ~, oscillRecon_ref_max_ampl] = ...
    RecognizeOscillations(x, [smooth_degree, smooth_length, pks_thresh]);
%disp(size(recognized_omegas));
pool_prev_omegas = [];

% components + residue initialization 
components = [];            %zeros(length(x), length(recognized_omegas)+1);

% miscellaneous initializations of swarm filtering 
output_type  = 2;        % output_type = 2 -> sum, output_type = 1 -> median
memI         = 0.9; 
memII        = 0.9;            
output_coeff = 0.004;
R            = 1 * rms(x); 

% decomposition starts
s = x;
total_counter = 1;
while true
   
    % recognition of oscillatory content
    [inter_recognized_omegas, amps_rec_omega, ~, NoOscillRecon] = ...
        RecognizeOscillations(s, [smooth_degree, smooth_length, pks_thresh], oscillRecon_ref_max_ampl);
    
    %disp('----')
    % Termination Criterion 
    if NoOscillRecon == 1 || isempty(inter_recognized_omegas) || chi2gof(s) == 0
        break;
    end
   
    current_recognized_omega = inter_recognized_omegas(amps_rec_omega == max(amps_rec_omega));
    disp(['all=', num2str(inter_recognized_omegas), 'current=', num2str(current_recognized_omega)])
    %disp('inter_rec= '); disp(inter_recognized_omegas);
    %disp('curr_rec= '); disp(current_recognized_omega);
    [NoM, delta] = SwarmParameterEstimator(current_recognized_omega(1));
    disp(['NoM= ', num2str(NoM)]); %disp(['d= ', num2str(delta)]);
    intermediate_components = zeros(length(x), 1);

    % recursively finding a subcomponent
    counter = 1;
    swarm_input = s;
    while true
        
        swarm_params = [NoM, delta, memI, memII, output_coeff, R];
        if length(swarm_params)<6
            swarm_status = 0;
            break;
        end
        
        %disp(swarm_params)
        swarm_output = SwarmFiltering(swarm_input, swarm_params, output_type, 0);

        [swarm_output_alligned, ~] = AlignSignals(swarm_input, swarm_output, 2);
        SD = sum((swarm_input - swarm_output_alligned).^2)/sum(swarm_input.^2);
           
        if SD < component_termination_crit || counter == 5
            intermediate_components(:, 1) = swarm_output;
            break;
        else
            swarm_input = swarm_output;
            counter = counter + 1;
        end
    end
       
    if swarm_status == 0
        break;
    end
    
    % delay estimation
    [aligned_intermediate_components, ~] = AlignSignals(s, intermediate_components(:, 1), 2);
        
    % peeling
    s = s - aligned_intermediate_components;
    %figure; plot(s);
    
    % attribute result to component
    %disp(pool_prev_omegas)
    [binCount, binIdx] = histc(pool_prev_omegas, [0.90*current_recognized_omega 1.1*current_recognized_omega]);
%     disp(['binCount=', num2str(binCount)]);
%     disp(['binIdx=', num2str(binIdx)]);
    [~, idx_dist] = min(abs(pool_prev_omegas - current_recognized_omega));

    if  isempty(pool_prev_omegas) 
        
        pool_prev_omegas = [pool_prev_omegas current_recognized_omega(1)];
        components = [components aligned_intermediate_components];
        [pool_prev_omegas, sortIdx] = sort(pool_prev_omegas);
        components = components(:, sortIdx);
        
    elseif binCount == 0
        
        pool_prev_omegas = [pool_prev_omegas current_recognized_omega(1)];
        components = [components aligned_intermediate_components];
        [pool_prev_omegas, sortIdx] = sort(pool_prev_omegas);
        components = components(:, sortIdx);
        
    else
        
        Idx = find(binIdx ~= 0);
        if length(Idx) > 1
            Idx = idx_dist;
        end
        % disp(['IDX= ', num2str(Idx)])
        % disp(['Idx=', num2str(Idx)]);
        % disp(['comp= ', num2str(size(components))]);
        % disp(['allig= ', num2str(size(aligned_intermediate_components))]);
        components(:, Idx) = components(:, Idx) + aligned_intermediate_components;
        
    end
    
%     omega_distances = abs(recognized_omegas - current_recognized_omega(1));
%     [~, s_omega_idx] = min(omega_distances);
%     components(:, s_omega_idx) = components(:, s_omega_idx) + aligned_intermediate_components;

    total_counter = total_counter + 1;
    if total_counter == Inf
        break;
    end
    
end

if swarm_status ~= 0
    components(:, end+1) = s;     % residue
else
    components = 0;
    pool_prev_omegas = 0;
end

% -------------------------------------------------------------------------
% -----------------      Auxillary Functions    ---------------------------
% -------------------------------------------------------------------------

function [aligned_signal, delay] = AlignSignals(initial_signal, delayed_signal, type)

% -------------------------------------------------------------------------
%   alignSignals: signal alignment and delay estimation
%   *** not sure if it works properly ***
% -------------------------------------------------------------------------
if nargin == 2
    type = 1;
end
switch type
    case 1
        [Cyx, lags] = xcorr(delayed_signal, initial_signal);
        Cyx_a = hilbert(Cyx);
        [~, idx] = max(abs(Cyx_a));
        tg = abs(lags(idx));
        tc = atan2(imag(Cyx_a(idx)), real(Cyx_a(idx))) / (2*pi);
        disp(tc)
        delay = tg + tc;
        aligned_signal = [delayed_signal(delay+1:end); zeros(delay, 1)];
    case 2
        [C, lags] = xcorr(initial_signal, delayed_signal);
        [~, idx] = max(C);
        delay = abs(lags(idx));
        aligned_signal = [delayed_signal(delay+1:end); zeros(delay, 1)];
    otherwise
        disp('error in arguments')
        aligned_signal = delayed_signal;
        delay = 0;
end

function [recognized_omega, amplitudes_omega, ref_max_ampl, IsEmpty] = RecognizeOscillations(input, params, max_ampl)

% -------------------------------------------------------------------------
%   Welch based spectrum
% -------------------------------------------------------------------------

window = params(1);
noverlap = params(2);
pks_threshold = params(3);

x = input;
% X = fft(x); Sx = abs(X).^2;
% smoothedSx = (sgolayfilt(Sx, smooth_degree, smooth_length));
nfft = 2*2048;
Sx = pwelch(x, window, noverlap, nfft);
Sx = Sx(1:end-1);
smoothedSx = Sx;

%figure; plot(Sx)
F = linspace(0, 1 - 1 / length(Sx), length(Sx));

if nargin == 2
    ref_max_ampl = max(smoothedSx);
else 
    ref_max_ampl = max_ampl;
end
%figure;plot(F, smoothedSx(1:end/2)/ref_max_ampl);

% special case
if max(smoothedSx / ref_max_ampl) < pks_threshold
    
    IsEmpty          = 1;
    recognized_omega = 0;
    amplitudes_omega = 0;
    ref_max_ampl     = 0;
    
else
    
    [amplitudes_omega, locs] = ...
        findpeaks(smoothedSx/ ref_max_ampl, 'MINPEAKHEIGHT', pks_threshold, 'MINPEAKDISTANCE', 2);

    IsEmpty                 = 0;
    temp_recognized_omega   = F(locs);
    recognized_omega        = sort(temp_recognized_omega, 'ascend');
    
end

function [NoM, delta] = SwarmParameterEstimator(omega)

% -------------------------------------------------------------------------
%   paramEstimator: Summary of this function goes here
%       args: omega ->[0,1]
% -------------------------------------------------------------------------

% NoM = (3.908*omega + 0.68) ./ (omega.^3 - 0.4066*omega.^2 + 0.1662*omega - 0.004681);
% delta = -1.098*omega.^2 + 3.11*omega +0.03999;

NoM = round(33.46 * omega.^(-0.7347) - 29.1);
if omega > 0.82
    NoM = NoM - round(22.22*omega-16.22);
end
if NoM <= 1
    NoM = 2;
end
%NoM = round(38.19 * omega.^(-0.69) - 35.39);

%delta = -1.474*omega.^2 + 3.391*omega - 0.005;
delta = -1.5*omega.^2 + 3.454*omega - 0.01;
%delta = -0.912 * omega.^3 - 0.3117*omega.^2 + 3.219*omega - 0.004852;

function [output, swarmPos] = SwarmFiltering(input, params, type, waitBar)

% -------------------------------------------------------------------------
%   SWARMFILTERING Summary of this function goes here
%   Detailed explanation goes here
% -------------------------------------------------------------------------

parCounter = 1; NoM = params(parCounter);                       % Number of Members
k1 = 1; %parCounter = parCounter + 1; k1 = params(parCounter);           % External force's coeff: Fext = k1*diff + k0
k0 = 0; %parCounter = parCounter + 1; k0 = params(parCounter);           % External force's bias
parCounter = parCounter + 1; delta = params(parCounter);        % Delta parameter
parCounter = parCounter + 1; memTypeI = params(parCounter);     % Memory type I: memory in velocity
parCounter = parCounter + 1; memTypeII = params(parCounter);    % Memory type II: memory in position 
parCounter = parCounter + 1; outCoeff = params(parCounter);     % Coefficient multiplied to output
vmax = +Inf; %parCounter = parCounter + 1; vmax = params(parCounter);         % Maximum velocity of the Members of Swarm
parCounter = parCounter + 1; DoC = params(parCounter);          % Distance of Covergence
eta = 1; %parCounter = parCounter + 1; eta = params(parCounter);          % Coeff of Internal Force

% prey initialization
prey = input;% - mean(input);

% swarm initialization
swarmPos = zeros(length(prey), NoM);
A = abs(prey(1) - prey(2));
if A == 0
    A = 1;
end
if rem(NoM, 2) == 0
    dif = A / NoM;
    swarmPos(1, :) = [-A/2:dif:-dif, dif:dif:A/2];
else
    dif = (A / 2) / floor(NoM / 2);
    swarmPos(1, :) = -A/2:dif:A/2;
end
swarmPos(1, :) = swarmPos(1, :) + prey(1);
swarmVel = zeros(length(prey), NoM);

% swarm chasing
duration = length(prey);

if waitBar == 1
    wb = waitbar(0, 'wait');
end

for t = 2:1:duration
    
    tempDist = prey(t) - swarmPos(t-1, :);
    Fext = k1 * tempDist - sign(tempDist) * k0;
    Fint = zeros(1, NoM);
    for i = 1:1:NoM
        for j = 1:1:NoM
            y = swarmPos(t-1, i) - swarmPos(t-1, j); 
            if y == 0
                continue;
            else
                Fint(i) = Fint(i) + sign(y)* eta * log(abs(y) / DoC);
            end
        end
    end
    Fint = Fint / (NoM - 1);
    
    F = Fext + Fint;
    
    swarmVel(t, :) = memTypeI * swarmVel(t-1, :) + F * delta;
    for i = 1:NoM
        if swarmVel(t, i) > vmax
            swarmVel(t, i) = vmax;
        end
    end
    
    swarmPos(t, :) = memTypeII * swarmPos(t-1, :) + swarmVel(t, :) * delta;
    if waitBar == 1
        waitbar(t/duration, wb, 'wait...');
    end
end

if waitBar == 1
    close(wb);
end

if type == 1
    output = outCoeff * median(swarmPos, 2);
elseif type == 2
    output = outCoeff * sum(swarmPos, 2);
end
