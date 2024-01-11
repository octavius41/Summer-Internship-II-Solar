close all; clc;

fm = 50; T = 0.002; fs = 1/T; A = 10^-9; t=0:T:12;r = 200*10^9; duration = 12;
vf = -1500:1:1500; vb = 1500:-1:-1499; vt = [vf,vb];
N = length(t);
w = linspace(-fs/2,fs/2,N);

it = vt./r;
nt = A*sin(2*pi*fm*t);
dt = it + nt;

figure(1)
plot(t,dt); title("time-current graph of source signal"); xlabel("time(s)");ylabel("current(A)");legend("source(t)");
grid on;

ffin = fftshift(fft(dt,N));

figure(2)
plot(w,abs(ffin)/N);title("frequency spectrum of source signal"); xlabel("frequency(f)");ylabel("Amplitude");legend("source(f)");
xlim([-100 100])

stepsize = round(N/60);
tsamp = 0:T:(T*(stepsize-1));

rcos = cos(2 * pi * fm * tsamp);
rsin = sin(2 * pi * fm * tsamp);

cutoff_frequency = 2;
[f1num,f1denum] = cheby1( 1, cutoff_frequency/(fs/2), 0.008 );
gain = 0.00023585659679385248;
f1num = f1num.*gain;
f1denum = f1denum.*gain;

phases = [];
amplitudes = [];

varamp = 1.68;%length(t)^(-2)/(A*1.5*duration);% %(2/max(rcos))*((mean(in_phase_filtered)^2+(mean(quadrature_filtered)^2)));

deriv = zeros(1,N);
phases = zeros(1,N);

for i = 1:1:(length(t)-stepsize)
    temp = dt(i : i + stepsize - 1);
    in_phase_signal = temp .* rcos;
    quadrature_signal = temp .* rsin;

    in_phase_filtered = filter(f1num, f1denum, in_phase_signal);
    quadrature_filtered = filter(f1num, f1denum, quadrature_signal);

    phases(i) =  atan2(mean(in_phase_filtered), mean(quadrature_filtered));

    tempp = atan2(mean(in_phase_filtered), mean(quadrature_filtered));
    tempa = (varamp*(mean(in_phase_filtered)^2 + mean(quadrature_filtered)^2)^(1/2))/max(rsin);

    deriv(i) = tempa*sin(2 * pi * fm * tsamp(1) + tempp);
end

    temp = dt(1700 : 1700 + stepsize - 1);
    in_phase_signal = temp .* rcos;
    quadrature_signal = temp .* rsin;

    in_phase_filtered = filter(f1num, f1denum, in_phase_signal);
    quadrature_filtered = filter(f1num, f1denum, quadrature_signal);

figure(3)
plot(t,phases);title("time-phase relation of calculated noise data point");xlabel("time(s)");ylabel("radians(n2*pi)");legend("phase");


figure(4)
plot(tsamp,in_phase_filtered);
hold on;
plot(tsamp,quadrature_filtered);
hold on;
plot(tsamp,atan2(mean(in_phase_filtered),mean(quadrature_filtered))); legend("inphase","quadrature","phase");
hold off;

rest = dt - deriv;

ffout = fftshift(fft(rest,N));

figure(5);
plot(t,nt,"green")
hold on;
plot(t,deriv,"red"); title("implemented noise and derived noise comparison");xlabel("time(s)");ylabel("amplitude");legend("implemented","derived");%xlim([5.9 6.1])
hold off;

figure(6);
plot(t,rest);title("time dependent graph of output signal"); xlabel("time(s)");ylabel("current(A)");legend("output(t)");
xlim([5.5 6.5]);grid on;
figure(7)
plot(vt,rest);title("voltage dependent graph of output signal"); xlabel("voltage(mV)");ylabel("current(A)");legend("output(mV)");
grid on;
figure(8)
plot(w,abs(ffout)/N);title("frequency spectrum of output signal"); xlabel("frequency(f)");ylabel("Amplitude");legend("source(f)");
xlim([-100 100])
figure(9)
plot(w,abs(ffin)/N);title("frequency spectrum of source signal"); xlabel("frequency(f)");ylabel("Amplitude");legend("source(f)");




figure(10)
plot(abs(freqz(f1num,f1denum,fs/2)));title("frequency spectrum of moving average filter");xlabel("frequency(f)");ylabel("Amplitude");legend("H(f)");

figure(11)
plot(t,dt)
hold on;
plot(t,rest,"green");title("time dependent graph of source and output signal"); xlabel("time(s)");ylabel("current(A)");legend("source(t)","output(t)");
hold off;

figure(12)
plot(vt,dt)
hold on;
plot(vt,rest,"green");title("voltage dependent graph of source and output signal"); xlabel("voltage(mV)");ylabel("current(A)");legend("source(mV)","output(mV)");
hold off;

figure(13)
plot(t,dt);
hold on;
plot(t,deriv);title("time dependent source signal and derived noise to be substracted");xlabel("time(s)");ylabel("current(A)");legend("source(t)","derived noise(t)");
hold off;


%{
figure(14)
subplot(221)
plot(t,dt)
hold on;
plot(t,rest,"green");title("time dependent graph of source and output signal"); xlabel("time(s)");ylabel("current(A)");legend("source(t)","output(t)");
hold off;
subplot(222)
plot(vt,dt)
hold on;
plot(vt,rest,"green");title("voltage dependent graph of source and output signal"); xlabel("voltage(mV)");ylabel("current(A)");legend("source(mV)","output(mV)");
hold off;
subplot(223)
plot(w,abs(ffin)/N,"red");title("frequency spectrum of source signal"); xlabel("frequency(f)");ylabel("Amplitude");legend("source(f)");
xlim([-100 100])
subplot(224)
plot(w,abs(ffout)/N,"red");title("frequency spectrum of output signal"); xlabel("frequency(f)");ylabel("Amplitude");legend("output(f)");
xlim([-100 100])
%}




