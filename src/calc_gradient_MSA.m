function [MSA,freq] = calc_gradient_MSA(dsvfile)
%[MSA,freq] = calc_gradient_MSA(dsvfile)
%   ...

GRD = dsvread(dsvfile); 

min_time = max([GRD.time(1),0]);
max_time = min([GRD.time(end),Inf]);
min_sample = round(1+(min_time-GRD.time(1))/GRD.time_dwell);
max_sample = round(1+(max_time-GRD.time(1))/GRD.time_dwell);
GRDtime = GRD.time(min_sample:max_sample); %[s]
GRDval  = GRD.values(min_sample:max_sample); %[mT/m]

npts = numel(GRDval);
dt = GRDtime(2)-GRDtime(1);
df = 1/(npts*dt);

%calculate mean squared amplitude
MSA = abs(fftshift(fft(GRDval))/npts);        %[mT/m]
freq = (floor(-npts/2):floor(npts/2-1))*df;   %[Hz]

end

