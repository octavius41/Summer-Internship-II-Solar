close all; clc;

fm = 50; dt = 0.002; fs = 1/dt; A = 10^-9; t=0:dt:12;r = 200*10^9;
vf = -1500:1:1500; vb = 1500:-1:-1499; v = [vf,vb];
N = length(t);
w = linspace(-fs/2,fs/2,N);

i_f = v./r;
nois = A*sin(2*pi*fm.*t);

d = i_f + nois;

[f1num,f1denum] = butter( 4, [40/(fs/2),60/(fs/2)], 'bandpass' );
%[f1num,f1denum] = cheby1( 1, 0.01, [0.196,0.204] );
fgain = 0.910836493150565;
f1num = f1num*fgain;
f1denum = f1denum*fgain;
h = freqz(f1num, f1denum);

stepsize = (4*10^4*dt);
N = stepsize;
frequencies = (0:N-1) * (fs / N);
j = 1;
tempt = 0:dt:0.160;
filtered_out = [];

for i = 1:stepsize:(length(t)-1)
 
        temp = d(i:i+stepsize);
        tempf = filter(f1num, f1denum, d(i:i+stepsize));
        %p = rad2deg(atan2(imag(tempf((i+79)/80),real(tempf((i+79)/80)))));  %rad2deg(atan2(imag(tempf((i+69)/70),real(tempf((i+69)/70)))))
        %signal_values = d(i);
        fft_result = fft(tempf,N);
        [~, max_index] = max(abs(fft_result));
        tempF = frequencies(max_index);
        tempA = tempf(74);

        temp_noise = tempA*sin(2*pi*tempF*tempt);
        filtered_out = [filtered_out, d(i:i+stepsize) - temp_noise];
        j = j+1;
 
end

figure(1)
plot(t,filtered_out(75:6075));

figure(2)
plot(v,filtered_out(75:6075));



