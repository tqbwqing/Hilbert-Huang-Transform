function allmode = eemd(Y, Nstd, NE)

% -------------------------------------------------------------------------
% 
% Ensemble Empirical Mode Decomposition 
%
% INPUT:
%       Y: Inputted data;
%       Nstd: ratio of the standard deviation of the added noise and that of Y;
%       NE: Ensemble number for the EEMD
%
% OUTPUT:
%       A matrix of N*(m+1) matrix, where N is the length of the input
%       data Y, and m=fix(log2(N))-1. Column 1 is the original data, columns 2, 3, ...
%       m are the IMFs from high to low frequency, and comlumn (m+1) is the
%       residual (over all trend).
%
% NOTE:
%       It should be noted that when Nstd is set to 0 and NE is set to 1, the
%       program degenerates to a EMD program
%       
% -------------------------------------------------------------------------

xsize = length(Y);
dd    = 1:1:xsize;
Ystd  = std(Y);
Y     = Y / Ystd;

TNM   = fix(log2(xsize)) - 1;
TNM2  = TNM + 2;

allmode = zeros(xsize, TNM2);
for kk = 1:1:TNM2 
    for ii = 1:1:xsize
        allmode(ii, kk) = 0.0;
    end
end

for iii = 1:1:NE
    
    X1 = zeros(1, xsize);
    for i = 1:1:xsize,
        temp  = randn * Nstd;
        X1(i) = Y(i)+temp;
    end

    mode = zeros(xsize, 1);
    for jj = 1:1:xsize,
        mode(jj, 1) = Y(jj);
    end
    
    xorigin = X1;
    xend = xorigin;
    
    nmode = 1;
    while nmode <= TNM
        xstart = xend;
        iter = 1;
   
        while iter <= 10,
            [spmax, spmin, ~] = extrema(xstart);
            upper   = spline(spmax(:,1), spmax(:,2), dd);
            lower   = spline(spmin(:,1), spmin(:,2), dd);
            mean_ul = (upper + lower) / 2;
            xstart  = xstart - mean_ul;
            iter    = iter +1;
        end
        
        xend = xend - xstart;
   	    nmode = nmode+1;
        for jj = 1:1:xsize,
            mode(jj, nmode) = xstart(jj);
        end
    end
   
    for jj = 1:1:xsize,
        mode(jj, nmode+1) = xend(jj);
    end
   
    allmode = allmode + mode;
    
end

allmode = allmode / NE;
allmode = allmode * Ystd;

end

function [spmax, spmin, flag] = extrema(in_data)

flag  = 1;
dsize = length(in_data);

spmax(1, 1) = 1;
spmax(1, 2) = in_data(1);

jj = 2;
kk = 2;
while jj < dsize
    if (in_data(jj-1) <= in_data(jj)) && (in_data(jj) >= in_data(jj+1))
        spmax(kk,1) = jj;
        spmax(kk,2) = in_data (jj);
        kk = kk + 1;
    end
    jj = jj + 1;
end
spmax(kk, 1) = dsize;
spmax(kk, 2) = in_data(dsize);

if kk >= 4
    slope1 = (spmax(2, 2) - spmax(3, 2)) / (spmax(2, 1) - spmax(3, 1));
    tmp1   = slope1 * (spmax(1, 1) - spmax(2, 1)) + spmax(2, 2);
    
    if tmp1 > spmax(1, 2)
        spmax(1, 2) = tmp1;
    end

    slope2 = (spmax(kk-1, 2) - spmax(kk-2, 2)) / (spmax(kk-1, 1) - spmax(kk-2, 1));
    tmp2   = slope2 * (spmax(kk, 1) - spmax(kk-1, 1)) + spmax(kk-1, 2);
    if tmp2 > spmax(kk, 2)
        spmax(kk, 2) = tmp2;
    end
else
    flag=-1;
end

msize = size(in_data);
dsize = max(msize);

spmin(1, 1) = 1;
spmin(1, 2) = in_data(1);
jj = 2;
kk = 2;

while jj < dsize,
    if (in_data(jj-1) >= in_data(jj)) && (in_data(jj) <= in_data(jj+1))
        spmin(kk, 1) = jj;
        spmin(kk, 2) = in_data (jj);
        kk = kk + 1;
    end
    jj = jj + 1;
end
spmin(kk, 1) = dsize;
spmin(kk, 2) = in_data(dsize);

if kk >= 4
    slope1 = (spmin(2, 2) - spmin(3, 2)) / (spmin(2, 1) - spmin(3, 1));
    tmp1   = slope1 * (spmin(1, 1) - spmin(2, 1)) + spmin(2, 2);
    if tmp1 < spmin(1, 2)
        spmin(1, 2) = tmp1;
    end

    slope2 = (spmin(kk-1, 2) - spmin(kk-2, 2)) / (spmin(kk-1, 1) - spmin(kk-2, 1));
    tmp2   = slope2 * (spmin(kk, 1) - spmin(kk-1, 1)) + spmin(kk-1, 2);
    if tmp2 < spmin(kk, 2)
        spmin(kk, 2) = tmp2;
    end
else
    flag = -1;
end

end

