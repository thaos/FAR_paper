// Based on Geoffroy et al. (2013)

lines(0);
clear();
stacksize(100000000);

models_name = ["BCC";"CCCMA";"CNRM";"CSIRO";"GFDL";"INM";"IPSL";"MIROC";"MPIM";"MRI";"NCC"];

F    = [6.70; 7.60; 7.10; 5.10; 6.60; 6.20; 6.40; 8.50; 8.20; 6.60; 6.20]; 
lamb = [1.21; 1.03; 1.08; 0.61; 1.34; 1.51; 0.78; 1.58; 1.14; 1.26; 1.11];
gamm = [0.71; 0.69; 0.47; 0.92; 0.87; 0.66; 0.60; 0.64; 0.76; 0.73; 0.95];
c    = [7.80; 7.30; 10.0; 6.20; 8.10; 8.90; 9.90; 8.20; 7.30; 8.90; 7.80];
c0   = [52.0; 64.0; 130.; 68.0; 110.; 314.; 98.0; 134.; 70.0; 62.0; 100.];


// Two-box model for the climate system (ocean-atmosphere).
// c  dT/dt   = F - lamb T - gamm (T - T_0)
// c0 dT_0/dt = gamm (T - T_0)
// Based on Held et al. (2010)

//-- Arguments :
// forcing : step - linear - stabilization
// F_infty : amplitude of the forcing
// c       : atm/upper ocean heat capacity
// c0      : deep ocean heat capacity
// lamb    : feedback parameter
// gamm    : ocean heat uptake coefficient
// N       : number of points simulated   

function [T, To, H] = held_model(forcing, F_infty, c, c0, lamb, gamm, N)


//-0- Forcing functions
dt=1; //-- timestep (year)
n=1:N+1;
t=0:N;

if (forcing == "abrupt")
  FF = F_infty*bool2s(t>=0);
elseif (forcing == "linear")
  te = N;
  FF = F_infty*t/te;
elseif (forcing == "stabilization")
  te = N/2;
  ts = 0:te;
  FF = zeros(1,N+1);
  FF(1,1:te+1) = F_infty*ts/te;
  FF(1,te+2:$) = F_infty;
elseif (forcing == "periodic")
  omega = 2*%pi / 10 ;
  FF = F_infty*sin(omega*t);
else
  FF = zeros(1,N+1);
end

disp("etape 0");

//-1- Numerical solutions (explicit scheme)
T=zeros(1,N+1);
To=zeros(1,N+1);
H=zeros(1,N+1);
T(1)=0;
To(1)=0;
for n=2:N+1
  T(n) = (T(n-1) + dt/c*(FF(n-1) - lamb*T(n-1) - gamm*(T(n-1)-To(n-1))));
  To(n)= (To(n-1) + dt/c0*(gamm*(T(n-1)-To(n-1))));
  H(n) = gamm*(T(n) - To(n));
end
// T_num = To(2:$)';
T_num = To(2:$)';



//-2- Parameters for exact analytical solution
//-- Miscellaneaous
b     = (lamb + gamm)/c + gamm/c0;
bstar = (lamb + gamm)/c - gamm/c0;
delta = b*b - 4*(lamb*gamm)/(c*c0);

//-- Constants
phi_f = c/(2*gamm)*(bstar - sqrt(delta));
phi_s = c/(2*gamm)*(bstar + sqrt(delta));

//-- Relaxation times
tau_f = (c*c0)/(2*lamb*gamm)*(b - sqrt(delta));
disp(tau_f);
tau_s = (c*c0)/(2*lamb*gamm)*(b + sqrt(delta));
disp(tau_s);

//-- ECS contributions
a_f   =  (phi_s*tau_f*lamb)/(c*(phi_s - phi_f));
a_s   = -(phi_f*tau_s*lamb)/(c*(phi_s - phi_f));

//-- Relations verification
// disp([a_f + a_s; phi_f*a_f + phi_s*a_s; tau_f*a_f + tau_s*a_s - (c + c0)/lamb; tau_f*(1 - a_f) + tau_s*(1 - a_s) - c0/gamm; c + phi_f*c0 - lamb*tau_f; c + phi_s*c0 - lamb*tau_s;]);



//-3- Exact analytical solutions
//-- Exact analytical solution for step forcing
if (forcing == "abrupt")
  T_h_ana = F_infty/lamb * ( a_f*(1-exp(-t/tau_f)) + a_s*(1-exp(-t/tau_s)) );
  T_ana = T_h_ana(2:$)';
//-- Exact analytical solution for linear forcing
elseif (forcing == "linear")
  T_l_ana = F_infty/lamb*(t/te) - F_infty/lamb*(1/te) * ( a_f*tau_f*(1-exp(-t/tau_f)) + a_s*tau_s*(1-exp(-t/tau_s)) ) ;
  T_ana = T_l_ana(2:$)';
//-- Exact analytical solution for linear+stabilization forcing
elseif (forcing == "stabilization")
  tbef = 0:te;
  taft = te+1:N;
  T_s_ana = zeros(1,N+1);
  T_s_ana(1:te+1)   = F_infty/lamb*(tbef/te) - F_infty/lamb*(1/te) * ( a_f*tau_f*(1-exp(-tbef/tau_f)) + a_s*tau_s*(1-exp(-tbef/tau_s)) ) ;
  T_s_ana(te+2:N+1) = F_infty/lamb - F_infty/lamb*(1/te) * ( a_f*tau_f*(1-exp(-te/tau_f))*exp(-(taft-te)/tau_f) + a_s*tau_s*(1-exp(-te/tau_s))*exp(-(taft-te)/tau_s) ) ;
  T_ana = T_s_ana(2:$)';
elseif (forcing == "periodic")
  T_p_ana = F_infty * ( a_f / (lamb*(1 + (omega*tau_f)^2)) * ( sin(omega*t) - omega*tau_f*cos(omega*t) ) + a_s / (lamb*(1 + (omega*tau_s)^2)) * ( sin(omega*t) - omega*tau_s*cos(omega*t) ) );
  T_ana = T_p_ana(2:$)';
else
  T_ana = 0;
end

if (forcing == "bode")

omega_=0.0001:0.01:1000;
omega=omega_';

u_f = a_f ./ ( 1 + (omega*tau_f).^2);
u_s = a_s ./ ( 1 + (omega*tau_s).^2);

// disp(size(u_f));
// disp(size(u_s));

disp(size((u_f + u_s).^2));
disp(size(1 + omega.^2));

modH = (u_f + u_s)^2 .* ( 1 + omega.^2 );

disp(size(modH));

G = 20*log(modH);

xset("window",0);
plot(omega,G,"black");
f=gcf();
a=gca();
// a.thickness=6;
a.log_flags="lnn";
// disp(a);

plot(omega,-40*log(omega) + 20*log(a_f/tau_f^4),"red");

a.data_bounds(2,2) = 10;

disp(2*%pi*tau_f^2/a_f);

// xs2eps(0,"bode","p");


end

endfunction
