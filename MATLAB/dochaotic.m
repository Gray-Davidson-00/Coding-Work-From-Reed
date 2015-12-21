function [t,omega,domega,th,xh,Esigma,Erho,Ek]=dochaotic()
%
%  Read the m-file
%
%  Some constants that need to be adjusted 
%  when dealing with the data
%
 transientcut  = 15000; % how many initial points to cut
 fftcut        = 2000; % how many high freq. points to keep
 edgecut       = 10000;  % the fft for nonperiodic functions has
                        % edge effects - we just cut of the edges 
                        % first few and last few points.
 endindex      = 180000;% there is timing jitter in the data so
                        % I just keep the first part of the data
                        % [For the part after the jitter use
                        % [ transientcut>185000 , endindex=-1]
 %
 % I think the data is in counts, this needs to be converted to angles 
 % Alison -- find out what the conversion is. I vaguely recall 4000 counts
 % per turn ??
 %
 anglefactor   = 2*pi/4000;


%
% Read data
%
if nargin<1
   %x=load('../Alison/PeriodicRun1/periodic.txt');
   x=load('MWWheel19Nov11_2249_1520');
else
    x=load(flnm);
end
endindex  = length(x(:,1));
%
% find and correct for the jumps due to integer 
% overflow of the decoder
%
% I) down jumps
%
diffx=x(2:end,2)-x(1:end-1,2);  
ijmp=find(diffx> 65000);
for i=1:length(ijmp)
    x(ijmp(i)+1:end,2)=x(ijmp(i)+1:end,2)-2^16;
end
%
% II) up jumps
%
diffx=x(2:end,2)-x(1:end-1,2);
clear ijmp;
ijmp=find(diffx< -65000);
for i=1:length(ijmp)
    x(ijmp(i)+1:end,2)=x(ijmp(i)+1:end,2)+2^16;
end
%
% cut off transient
%
if endindex>0 
    xx=x(transientcut:endindex,:);
else
 xx=x(transientcut:end,:);
end
 showall=1;
fig=0;
if showall==1
   fig=fig+1;figure(fig);clf;
   plot(x(:,1),x(:,2),xx(:,1),xx(:,2))
   xlabel 'Time'
   ylabel 'Counts ~ Angle'
   title 'Step 1 ; Transient Cut'
end
clear x
x=xx;
%
% Taking derivatives is tricky because it amplifies noise
%
%
% There is jitter in the timing
% 
  % Get some statistic about the timing jitter and recreate a proper time
  % vector
    dt=x(2:end,1)-x(1:end-1,1);
    mindt=min(dt); meandt=mean(dt); maxdt=max(dt);
    disp(['Times: mean(dt)=' num2str(meandt) ' min(dt)=' num2str(mindt) ' max(dt)=' num2str(maxdt)])
    dtav=round(10*mean(dt))/10; % average rounded to 100 of ms
    if (abs(mindt-meandt)/meandt > .1)||(abs(maxdt-meandt)/meandt>.1)
        disp('Too much time jitter')
        plot(dt)
        %return
    end
   % For now we are not resampling and assume that the variance in dt is small
   % so that we can just take it at face value and create an evenly spaced
   % time vector and assign to each time the 'Closest' data point
   N=length(x(:,1));
   if (mod(N,2)~=0) % the number is odd but we need even number
       N=N-1;
   end
   tmax = x(N,1)-x(1,1);
   t = linspace(0, tmax, N+1); t = t(1:length(t)-1)';
   theta=x(1:N,2)*anglefactor;
   
 plot(t,theta)
   %
 % There is measurement noise in the angle
 %
 % Now we smooth the data and take derivatives by using fft
 % and cutting off the higher frequencies that are there because
 % we oversample
 %
 %
 % Where to set the cutoff is nontrivial and you have to plot the data and check
 % that the data is not changed but noise is supressed 
 %
 %
 xfft=fft(theta,N);

    fig=fig+1;figure(fig);clf;
  semilogy(2*abs(xfft),'b');
  xlabel 'index'
  ylabel 'fft of x'
  
  
  xfft(fftcut+1:end-fftcut)=0;
  hold on
  semilogy(2*abs(xfft),'r'); 

  % smoothed variable
  theta=real(ifft(xfft));
  
  % first and second derivative
    omega  = fourdifft(real(theta),1,t(end));
    domega = fourdifft(real(theta),2,t(end));
  
  
    
  % Cutoff the edge to get rid of edge effects
    t      = t(edgecut+1:end-edgecut);
    theta  = theta(edgecut+1:end-edgecut);
    omega  = omega(edgecut+1:end-edgecut);
    domega = domega(edgecut+1:end-edgecut);
   
 %
 %
 % restrict time
    i      = find((t>200)&(t<1200));
    %i      = find((t>920)&(t<980));
    %i      = find((t>1230)&(t<1280));
    t      = t(i)-t(i(1));
    theta  = theta(i);
    omega  = omega(i);
    domega = domega(i);
   
    
    fig=fig+1;figure(fig);clf
    subplot(2,1,1)
    plot(t,theta,'r','MarkerSize',3)
    %  plot(x(1:end-1,1)-x(1,1),diff(x(:,2))*anglefactor/(t(2)-t(1)),'+r',t,omega,'b','MarkerSize',3)
    subplot(2,1,2)
    plot(t,omega,'b','MarkerSize',3)
 
    
      fig=fig+1;figure(fig);clf
    %subplot(2,2,3)       % delay embedding - the delay here is 4 seconds
    %n=round(2/dtav);imax=max(length(t),10000);
    plot(theta,omega)
    str=['\omega'];
    ylabel(str)
    xlabel '\theta'
    grid
    
  
  return


   function Dmf = fourdifft(f,m,tmax)

% Dmf = fourdifft(f,m) computes the m-th derivative of the function
% f(x) using the Fourier differentiation process.   The function 
% is assumed to be 2pi-periodic and the input data values f should 
% correspond to samples of the function at N equispaced points on [0, 2pi).
% The Fast Fourier Transform is used.
% 
%  Input:
%  f:      Vector of samples of f(x) at x = 0, 2pi/N, 4pi/N, ... , (N-1)2pi/N
%  m:      Derivative required (non-negative integer)
%
%  Output:
%  Dmf:     m-th derivative of f
%
%  S.C. Reddy, J.A.C. Weideman 2000. Corrected for MATLAB R13 
%  by JACW, April 2003. 
%
%  Scaled with respect to interval tmax by Lucas Illing, 2010


     f = f(:);                       % Make sure f is a column vector
     N = length(f);

     N1 =  floor((N-1)/2);           % Set up wavenumbers           
     N2 = (-N/2)*rem(m+1,2)*ones(rem(N+1,2));
   wave = [(0:N1)  N2 (-N1:-1)]';

%Dmf = ifft(((i*wave).^m).*fft(f));   % Transform to Fourier space, take deriv,
                                     % and return to physical space.
Dmf = ifft(((i*wave*2*pi/tmax).^m).*fft(f)); 
if max(abs(imag(f))) == 0; Dmf = real(Dmf); end  % Real data in, real derivative out
                                       
