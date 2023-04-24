%Propagation of a group of waves
m    = 1;                             %mass
hbar = 6.626e-34/(2*pi);              %constant
t    = linspace(0,10e-31, 200);       %time interval
x    = linspace(0,2*pi*200e-33,200);  %x interval
E    = normrnd(1, 0.01, 1, 1000);     %Energy measurement with accuracy
w    = E/hbar;                        %Frequency
k    = sqrt(2*m*E)/hbar;              %k


[T,K,X] = meshgrid(t, k, x);
[T,W,X] = meshgrid(t, w, x);
PSI  = squeeze(sum(exp(1i*(K.*X - W.*T))));
PROB = PSI.*conj(PSI);
[X,T] = meshgrid(x, t);
surf(X, T, PROB);
xlabel('x');
ylabel('t');
zlabel('non normalized probability');


figure;
t    = 2e-31;                          %constant time
x    = linspace(0,2*pi*200e-33,2000);  %x interval
psi  = sum(exp(1i*(k'*x - w'*t)));
plot(x, psi.*conj(psi));




