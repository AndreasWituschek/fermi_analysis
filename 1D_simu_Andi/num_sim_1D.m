clear
tic;

BaseName = '2020-04-03'; % date
%FEL harmonic
h = 5;
amp = 0.85;
satur = 1;

%integration stepsize
Delta_t = 10;

% E-field amplitude (peak)
E_0 = 1*10^6; % in V/m

%%%%Wavelenght of laser%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%WDp = 787; % Pulse carrier freq., in nm for Rb D Lines
WDp = 52.235; % Pulse carrier freq., in nm for He 1s4p line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%pulse parameters (chirp, pulseduration)

ch2 = 0.0; % quadratic chirp in fs^2, here: is loaded from input files
ch3 = 0.0; % qubic chirp, fs^(-3), 
TDp = 99.0; % seed pulse duration (FWHM of intensity) in fs, for transform limited pulse

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discretization steps and time propagatios 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_tau12=200; % number of points along tau_12 (has to be even)
step_tau = 5; % step size along tau in fs


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some conversion units & parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CmeV = 0.1239842/1000.;     % cm-1 ----> eV
wlam=1239.8424121;      % nm --- eV 
dh =  1.5193;           % 10/hbar
Deb = 3.33564*10^(-30); % Debye in SI
Ev = 1.6021766e-19;         % eV in SI
speed_of_light = 299792458.0; % in m/s
epsilon_0 = 8.854187817620389*10^(-12); % in F/m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%for Runge-Kutta integration%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rk=4; %% rk = 1,2,3,4 ----- order of the Runge-Kutta scheme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Pulse parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pulse is detemined in the frequency domain later in the script via the following definiton: 
% E(w) = exp(-TDp^2/(8*log(2))*w^2 - 1i/2*ch2*w^2- 1i/6*ch3*w^3) 

TDe = TDp * sqrt(2);  % Pulse duration (FWHW of amplitude), fs 

coeff= 5;
% A number, which determines the boundary when the pulse amplitude is numerically zero: E(t=TDe*coeff) = 0;  
% to save time, this could be adusted to t_FWHM_chirped (for 1510fs^2
% and TDp = 23, 15 is fine, not yet checke how small one can choose
% that for example for 500fs^2 and so on, for unchirped 5 is fine (TDp
% =23fs
WWp  = sqrt(8*log(2))/TDp * 5; % domain of the integration for the conversion of the pulse envelope into the time domain.  
stWW = WWp/100;                % integration step for the conversion of the pulse envelop from the frequency domain to the time domain.  
NWW = 2*round(WWp/stWW)+1; % number of steps used for the calculation of the pulse (conversion from freq to time domain?)
% IMPORTANT
% In the code, the pulse enevelop in the time domain, E(t), is numerically calculated from E(w) and normalized to 1, E(t=0)=1. 
% For purely Gaussian pulses, coeff=3.5 is enough; for chirped pulses, it should be higher. For Schmales  Spektrum, e.g., it should be around 15. 
% The code plots pulse envelops in the time & frequency domain. If they look strange, increase coeff and/or decrease stWW.  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hamiltonian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N0=1; % # of levels in the ground state
N1 = 1; % # levels in singly excited state manifold
NN = N0 + N1; % total # of levels

%Energies of states in (cm^-1)
%E1=0; %Rb 5S1/2
%E2=12816.545; % Rb 5P3/2
E1=0; %Rb 5S1/2
E2=191492.711; % Rb 5P3/2

%create Hamiltionian for system 

H0 = zeros(N0,N0); % initialise part of System Hamiltonian belonging to ground state
H0(1,1) = E1;
H1 = zeros(N1,N1); % initialise part of System Hamiltonian belonging to singly excited manifold
H1(1,1) = E2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transitio dipole moments and up/down transition operator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X = zeros(NN,NN); 

% give the values of the transition dipole moments here in Debye
% X is the up-transition operator 
% the values are given in units of the tranistion dipole moment from the ground 
% state to the lowest excited state

%%%%%%%%%%%%%%%% My definition %%%%%%%%%%%%%%%%%%%%%%

%Mu21R =   10.6586; % from Rb ground S1/2 ---> P3/2
Mu21R = 0.3328; % for He 1s --> 4p
X(1,2) = Mu21R/Mu21R; % up transition operator 

Hdeb = Mu21R; % is the real value of the transition dipole moment X(1,2) in Debye; 
wp=wlam/WDp; % Pulse carrier frequency, eV

Z0 = eye(N0,N0);
Z1 = eye(N1,N1);

H0 = H0*CmeV; % Hamiltonian in eV
H1 = H1*CmeV; % Hamiltonian in eV

% Hamiltonian in the rotating frame
H1 = H1-Z1*wp;
H = blkdiag(H0,H1);

Xd=X'; % down transition operator?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The main body of the code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tau1=0.; tau2=0.; 

%F1_true = F1; F2_true = F2; F3_true = F3; F4_true = F4;
fi1=0.; fi2=0.;

%tau_initial=500*0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Rho_fin_loop = zeros(NN,size(E0_vec,1));
%Popu_fin_loop = zeros(NN,size(E0_vec,1));

%RP_signal = zeros(size(E0_vec,1),max_tau12,max_tau34,maxT); % initialize matrix for RP signal
%NRP_signal = zeros(size(E0_vec,1),max_tau12,max_tau34,maxT); % initialize matrix for NRP signal
%FiAm_vec = zeros(size(E0_vec,1),1); % initialize array for FiAm values (to save at the end to know which parameters have been used)
%Delta_t_vec = zeros(size(E0_vec,1),1);% initialize array for Delta_t values (to save at the end to know which parameters have been used)

% calculate FiAm out of E_0
FiAm = 0.5*E_0*Mu21R*Deb/Ev; % in eV 
%The meaning of FiAm:
% It gives the energy of interaction of the laser field with the transition dipole moment X(1,2); 
% FiAm = 0.0001 corresponds to weak pulses; Strong-field effects start from FiAm = 0.01; 
% To recalculate all this stuff into the REAL values of the tranisition dipole moments you should do the following: 
% Let 
%HDeb= Mu21R; % is the real value of the transition dipole moment X(1,2) in Debye; 
%Then 
%E_field = FiAm * Ev / (HDeb * Deb)
% is the real amplitude of your lase field in V/m which corresponds to FiAm
disp('Is FiAm = '+string(FiAm) + ' << 0.01 ?')
F1=FiAm; F2=FiAm;

Delta_tR = Delta_t*dh; %??

%initial condition (population in ground state before laser pulses)
RhoMuNu = zeros(NN,1);
RhoMuNuIn = zeros(NN,1);
RhoMuNuIn_0 = zeros(NN,1);
RhoMuNuIn_0(1,1) = 1.;

ttau12 = zeros(max_tau12,1); %initialize tau vector


Popi = zeros(max_tau12,1); %initialize vector for population after interaction with the two pulses
ac = zeros(max_tau12,1); %initialize vector for AC intensity

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for Runge-Kutta-Integration
Ark = zeros(4,4); Brk = zeros(4); Crk = zeros(4);
KK1 = zeros(NN,1); KK2 = zeros(NN,1); 
KK3 = zeros(NN,1); KK4 = zeros(NN,1); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(rk == 1) 
    Brk(1)=1.; 
end

if(rk == 2)
    Ark(2,1)=0.5;  Brk(2)=1.; Crk(2)=0.5; 
end

if(rk == 3) 
    Ark(2,1)=0.5; Ark(3,1)=-1.; Ark(3,2)=2.;
    Brk(1)=1./6.; Brk(2)=4./6.; Brk(3)=1./6.;
    Crk(2)=0.5; Crk(3)=1.; 
end

if(rk == 4) 
    Ark(2,1)=0.5; Ark(3,2)=0.5; Ark(4,3)=1.;
    Brk(1)=1./6.; Brk(2)=2./6.; Brk(3)=2./6.; Brk(4)=1./6.;
    Crk(2)=0.5; Crk(3)=0.5; Crk(4)=1.;
end


norm = sat_harm(1,h,satur,amp);
iterator = 0:19;
for value = iterator %defining 10 phase steps to cycle phase from 0 to to 2pi
    phase = 1/length(iterator)*2*pi*value;
    for itau12 = 1:max_tau12+1 %start of tau loop
        tau12 = ((itau12-1)-max_tau12/2)*step_tau;
        ttau12(itau12,1) = tau12;



        tau1 = -tau12; %center of pulse 1
        tau2 = 0; %center of pulse 2

        ttt = 2*coeff*TDe+tau12; %? range for integration of schr?dinger equation in fs?
        tau0 = -coeff*TDe-tau12; % start point of integration?

        maxJ = round(ttt/Delta_t) + 300*0; %? maximum number of steps for integration of Schroedinger equation?

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Pulses in the time domain
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        maxJ2 = 2*maxJ+3;
        Pulse1 = zeros(maxJ2,1);
        Pulse2 = zeros(maxJ2,1);
        PulseT1 = zeros(maxJ2,1);
        PulseT2 = zeros(maxJ2,1);
        TRK = zeros(maxJ2,1);
        %PulseT = zeros(maxJ,1); % delete?
        PulseW = zeros(NWW,1);
        WWW = zeros(NWW,1);

        for jt2 = 1:maxJ2
            tt2 = tau0 + (jt2-1.)*Delta_t/2.0;
            Pulse1(jt2,1) = exp(-2*log(2)*((tt2-tau1)/TDp)^2)*exp(1i*phase); %definition of pulse directly in TD
            Pulse2(jt2,1) = exp(-2*log(2)*((tt2-tau2)/TDp)^2);

    % Definition of pulses in Frequency domain:
             TRK(jt2,1) = tt2;
    %                     
    %         A1 = 1; A2=1;
    %         if (abs(tt2-tau1) > coeff*TDe) 
    %             A1=0; 
    %         end
    %         if (abs(tt2-tau2) > coeff*TDe) 
    %             A2=0; 
    %         end
    % 
    %                     
    %         for iw=1:NWW
    %             WW = (iw-1)*stWW - WWp;
    %             WWW(iw,1)=WW;
    %             EW = exp(-TDp^2/(8*log(2))*WW^2 - 1i/2*ch2*WW^2 - 1i/6*ch3*WW^3); 
    %             PulseW(iw,1)=EW;
    % 
    %             Pulse1(jt2,1)=Pulse1(jt2,1)+A1*EW*exp(1i*WW*(tt2-tau1)); 
    %             Pulse2(jt2,1)=Pulse2(jt2,1)+A2*EW*exp(1i*WW*(tt2-tau2));
    %         end % end iw loop
        end % end jt2 loop

        %Normalization
        Pu1 = max(abs(Pulse1));
        Pu2 = max(abs(Pulse2));

        Pulse1=Pulse1/Pu1; %*(1-tau1/500*0.2);
        Pulse2=Pulse2/Pu2;
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % saving of pulse in time and freq domain (normalized)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if ((itau12==1))
            fn = char(BaseName + "_amp" + string(amp*100) + "_pulse");
            FileName_pulse = ['//mpmnsh01.public.ads.uni-freiburg.de/mpmnsh01/user/Projekte/BMBF/FERMI/DataAnalysis/Beamtime2/combined/scan_031/simulations/phase_cycled/',fn,'.mat'];
            save(FileName_pulse, 'TDp', 'WDp', 'ch2', 'ch3', 'WWW', 'PulseW', 'TRK', 'Pulse1')
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Initial conditions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        RhoMuNuIn = RhoMuNuIn_0;

        tm = zeros(maxJ,1); % new preallocating

        pulses = zeros(maxJ,1);

        for jt = 1:maxJ % start of jt loop
            tm(jt) = tau0 + (jt-1.)*Delta_t;
            tc=tm(jt);

            FX = F1*sat_harm(Pulse1(2*jt-1,1)+Pulse2(2*jt-1,1),h,satur,amp)/norm;
            pulses(jt) = FX(1);
            FXd = FX';
            SFt = -X*FX -Xd*FXd + H;

            RhoMuNuIn_K = RhoMuNuIn;

            KK1 = -1i*SFt*RhoMuNuIn_K;

            if (rk > 1)
                tc=tm(jt)+Delta_t*Crk(2);

                FX = F1*sat_harm(Pulse1(2*jt-1,1)+Pulse2(2*jt-1,1),h,satur,amp)/norm;
                FXd = FX';
                SFt = - X*FX - Xd*FXd + H;
                RhoMuNuIn_K = RhoMuNuIn + Delta_tR*Ark(2,1)*KK1;
                KK2 = -1i*SFt*RhoMuNuIn_K;
            end

            if (rk > 2)
                tc=tm(jt)+Delta_t*Crk(3);

                FX = F1*sat_harm(Pulse1(2*jt-1,1)+Pulse2(2*jt-1,1),h,satur,amp)/norm;
                FXd = FX';
                SFt = - X*FX - Xd*FXd + H;
                RhoMuNuIn_K = RhoMuNuIn + Delta_tR*(Ark(3,1)*KK1+Ark(3,2)*KK2);
                KK3 = -1i*SFt*RhoMuNuIn_K;
            end

            if (rk > 3)
                tc=tm(jt)+Delta_t*Crk(4);

                FX = F1*sat_harm(Pulse1(2*jt-1,1)+Pulse2(2*jt-1,1),h,satur,amp)/norm;
                FXd = FX';
                SFt = - X*FX - Xd*FXd + H;
                RhoMuNuIn_K = RhoMuNuIn + Delta_tR*(Ark(4,1)*KK1+Ark(4,2)*KK2+Ark(4,3)*KK3);
                KK4 = -1i*SFt*RhoMuNuIn_K;
            end

            if (rk == 1) 
                Rho_fin = RhoMuNuIn +  Delta_tR*Brk(1)*KK1;
            end 
            if (rk == 2) 
                Rho_fin = RhoMuNuIn +  Delta_tR*(Brk(1)*KK1+Brk(2)*KK2);
            end 
            if (rk == 3) 
                Rho_fin = RhoMuNuIn +  Delta_tR*(Brk(1)*KK1+Brk(2)*KK2+Brk(3)*KK3);
            end
            if (rk == 4) 
                Rho_fin = RhoMuNuIn +  Delta_tR*(Brk(1)*KK1+Brk(2)*KK2+Brk(3)*KK3+Brk(4)*KK4);
            end
            RhoMuNuIn = Rho_fin;
        end % end of jt loop

        %Popu excited state
        Pop_flu = 0;
        for jj =1:N1
            Pop_flu = Pop_flu + Rho_fin(jj+1,1)*Rho_fin(jj+1,1)';
        end
        Popi(itau12,1) = Pop_flu;
        ac(itau12,1) = sum(abs(pulses).^2); 
    end % end of tau loop

%      % direct data plotting
%      figure('Name', 'Population vs tau12')
%      plot(ttau12,Popi, 'b')
%      title(phase)
%     figure('Name','Pulses');
%     plot(1:length(pulses),real(pulses)/F1)
    %data saving
    fn = char(BaseName + "_amp" + string(amp*100) + "_" + string(value));
    FileName = ['//mpmnsh01.public.ads.uni-freiburg.de/mpmnsh01/user/Projekte/BMBF/FERMI/DataAnalysis/Beamtime2/combined/scan_031/simulations/phase_cycled/',fn,'.mat'];
    save(FileName, 'Popi', 'ttau12', 'step_tau', 'max_tau12', 'wp','Delta_t', 'FiAm','Mu21R')

    fn = char(BaseName + "_amp" + string(amp*100) + "_" + string(value) + "_ac");
    FileName = ['//mpmnsh01.public.ads.uni-freiburg.de/mpmnsh01/user/Projekte/BMBF/FERMI/DataAnalysis/Beamtime2/combined/scan_031/simulations/phase_cycled/',fn,'.mat'];
    save(FileName, 'ac', 'ttau12', 'step_tau', 'max_tau12', 'wp','Delta_t', 'FiAm','Mu21R')
end

% figure('Name','Pulses');
% plot(1:length(pulses),real(pulses)/F1)


toc;

% %saturated harmonic generation function:
% function hg = sat_harm(x,h,satur,amp)
%     if satur == 1
%         hg = (amp*x)^h/(1+abs((amp*x)^h));
%     else
%         hg = x^h;
%     end
% end

% % saturated harmonic generation function:
function hg = sat_harm(x,h,satur,amp)
    if satur == 1
%         max_x_jv= h*(1+sqrt(2/3)*h^(-2/3));
%         max_jv = besselj(h,max_x_jv);
%         one = ones(length(x));
%         one(find(x<0)) = -1;
        seedenvelope=abs(x); %seed envelope
        seedphase=unwrap(angle(x)); %seed phase

        hg = besselj(h,h*amp*seedenvelope)*exp(1i*h*seedphase);
    else
        hg = x^h;
    end
end
% z = linspace(-2,2,100)+1i*linspace(-2,2,100);
% figure()
% plot(real(z),real(besselj(h,h*z,1)));
% figure()
% plot(real(z),imag(besselj(h,h*z,1)));
% figure()
% plot(real(z),besselj(h,h*real(z)));
% plot(1:length(Pulse1),Pulse1)
% hold on
% plot(1:923,Pulse2)
% plot(1:923,PulseT)
% hold off