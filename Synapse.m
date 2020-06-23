%------------------------------Complete-------------
function Synapse()    
%---PARAMETERS TO VARY----------------
numTrials=5; %number of stimuli
numCells=5;
numVesicles=10; %nymber of vesicles
prob= 0.2 ;   %probability of release
gEmax=1e-4;  %conductance produce when a single vesicle released; default=1.0e-4
Eexc=0;  % equilibrium potentials for synapses
NaOn=1%turns off sodium channels at 0
KOn=1;%turns off Potassium channels 

%-----------HH & Passive constants---------- ------------
delT=0.05;totalTime=60;Rec_length=round(totalTime/delT); T=delT*(1:Rec_length); %discreitzaiton parameters
gNamax=0.2;gKmax=3;R=260;tau=10;C=tau/R;gl=1/R;  %conductances
Erest=-70;EK=-90;ENa=70;
 
%-----------vectors and matrices-------------------
%Voltage and injected current
    V=zeros(1,Rec_length);
    avgV=zeros(1,Rec_length);   
    Iinj=zeros(1,Rec_length); %injected current
 
%Hodgkin/Huxley & passive parameters
    n=zeros(1,Rec_length); m=zeros(1,Rec_length);h=ones(1,Rec_length);
    IK=zeros(1,Rec_length);INa=zeros(1,Rec_length); Ileak=zeros(1,Rec_length);Iionic=zeros(1,Rec_length);  
    
%synaptic parameters
    gEsingle=zeros(1,Rec_length); %conductance time course from single vesicle released
    gE=zeros(1,Rec_length); %total conductance (=nVesicles*gEmax*gEsingle)
    
%----constants and vectors for analyses and plotting
pkPSP=zeros(1,numTrials);           %will contain list of amplitudes during the experiment
rec_point=20; %time point on trace where voltage is measured 
 
%-------main program-----------
figure(1);clf;%hold on;
gEsingle=gen_PSG(); %calculates the conductance from a single vesicle (see below for function)
for ii=1:numTrials
    gE=gE*0;
    for jj=1:numCells
    V=V*0+Erest;   %for each trial, V is set to resting potential 
    nse=normrnd(0,1e-4,[1,Rec_length]); %add noise for realism
    nVesicles=vesicles(); %calls function to calculate the number of vesicles released in trial (see below)
    gE = gE+nVesicles*gEmax*gEsingle+nse; %scales the alpha function by gEmax and nVesicles
    end
    HH();
 figure(1);plot(T,V);
        set(gca,'fontsize',14);title('membrane potential','FontSize',14);xlabel('time (ms)');ylabel('membrane potential (mV)','fontsize',14);
     %uses hodgkin huxley to generate voltage
  %  avgV=avgV+V; %compiles average voltage
   % pkPSP(ii)=V(round(rec_point/delT))-Erest; %measures the amplitude at t=rec_point
   
    end

% HH();
%  figure(1);plot(T,V);
%         set(gca,'fontsize',14);title('membrane potential','FontSize',14);xlabel('time (ms)');ylabel('membrane potential (mV)','fontsize',14);
% avgV=avgV/numTrials; %calculates 
% figure(1);plot(T,avgV,'k','LineWidth',5);
% figure(2);clf;histogram(pkPSP,10*numVesicles);% makes histogram of PSP amplitudes
%  
 
 
%------------Auxiliary functions------------------------------
    function HH()   %calculates voltage with HH equations
        for t=2:Rec_length
            [ntau,ninf,mtau,minf,htau,hinf]=calc_params(V(t-1));
            n(t)=n(t-1)+(ninf-n(t-1))*(delT/ntau);
            IK(t)=KOn*gKmax*n(t)^3*(V(t)-EK);
 
            m(t)=m(t-1)+(minf-m(t-1))*(delT/mtau);
            h(t)=h(t-1)+(hinf-h(t-1))*(delT/htau);
            INa(t)=NaOn*gNamax*m(t)^3*h(t)*(V(t-1)-ENa);
            Iexc(t)=gE(t)*(V(t-1)-Eexc);
            %Iinh(t)=gItrain(t)*(V(t-1)-Einh);
 
            Ileak(t)=gl*(V(t-1)-Erest);
            Iionic(t)=INa(t)+IK(t)+Ileak(t);
            Isyn(t)=Iexc(t);%+Iinh(t);
            V(t)=V(t-1)+(Iinj(t)-Iionic(t)-Isyn(t))*(delT/C);
        end
        function [ntau,ninf,mtau,minf,htau,hinf]=calc_params(v)
            ntaumin=0.05;
            alph = 0.03*(-v-65.5)/(exp((-v-65.5)/2)-1);
            beta = 12*exp(-(v+70.5)/5);     
            ntau=1/(alph+beta);
            ninf=alph/(alph+beta);
            ntau=max(ntau,ntaumin); %need to set lower bound or divide by near zero occurs
 
            mtaumin=0.05; htaumin=0.5; %fudge to prevent divide by zero
            alph = 0.4*(-v-60.5)/(exp((-v-60.5)/2)-1);
            beta = 64*exp(-(v+75.5)/5);     
            minf = alph/(alph+beta);
            mtau = 1/(alph+beta);   mtau=max(mtau,mtaumin);             
            alph = 0.28*exp(-(v+59)/2);
            beta = 4/(1+exp(-(v+10)/10));
            hinf = alph/(alph+beta);
            htau =1/(alph+beta);    htau=max(htau,htaumin);                 
        end
    end
 
%calculates # vesicles released per trial
    function num=vesicles()
     num=0;
     for j=1:numVesicles
      if (unifrnd(0,1)<prob)        %generates uniform random number between 0 and 1
          num=num+1;                %keep track of number of released vesicles
      end
     end
    end
function PSG=gen_PSG() %for calculating the time course of the synaptic conductance
    alph=0.2;
    PSG=zeros(1,Rec_length);pTau=4;
    s_delay=round(5/delT); %adds delay so can look at baseline voltage
    for i=s_delay+1:Rec_length
        PSG(i)=(((i-s_delay)*delT)*alph)*exp(-(i-s_delay)*delT*alph);   %alpha function
    end
    PSG=PSG/max(PSG); %n    ormalizes PSG to a peak of 1
end
        
end

