function [u_out,s_out] = Control(params,ctrl_params,states,input_store,desired)
%% Desired
phi_d       = desired(1);
theta_d     = desired(2);
psi_dot_d   = desired(3);
h_d         = desired(4);
phi_dot_d   = desired(5);
theta_dot_d = desired(6);
psi_ddot_d  = desired(7);
h_dot_d     = desired(8);

% phi_d     = 0;
% theta_d   = deg2rad(17.12);
% psi_dot_d = 0;
% h_d       = 0;

% phi_dot_d   = 0;
% theta_dot_d = 0;
% psi_ddot_d  = 0;
% h_dot_d     = 0;



%% Parameters
m      = params.m     ;
Ixx    = params.Ixx   ;
Iyy    = params.Iyy   ;
Izz    = params.Izz   ;
Ixz    = params.Ixz   ;
S      = params.S     ;
b      = params.b     ;
c      = params.c     ;
Tmax   = params.Tmax  ;
g      = params.g     ;
rho    = params.rho   ;
Gamma  = params.Gamma ;
Gamma1 = params.Gamma1;
Gamma2 = params.Gamma2;
Gamma3 = params.Gamma3;
Gamma4 = params.Gamma4;
Gamma5 = params.Gamma5;
Gamma6 = params.Gamma6;
Gamma7 = params.Gamma7;
Gamma8 = params.Gamma8;
%% Variables allocation
u     = states(1);
v     = states(2);
w     = states(3);
p     = states(4);
q     = states(5);
r     = states(6);
phi   = states(7);
theta = states(8);
psi   = states(9);
h     = states(12);

V = sqrt(u.^2 +v.^2 +w.^2);

alpha = atan(w/u);
beta  = asin(v/V);
%%
R = [1,  sin(phi)*tan(theta),  cos(phi)*tan(theta);
     0,  cos(phi),            -sin(phi);
     0,  sin(phi)/cos(theta),  cos(phi)/cos(theta)];


phi_dot   = R(1,:)*[p; q; r];
theta_dot = R(2,:)*[p; q; r];
psi_dot   = R(3,:)*[p; q; r];
h_dot     = -[-sin(theta), sin(phi)*cos(theta), cos(phi)*cos(theta)]*[u; v; w];

%% Aerodynamic coefficients
aero_out = aero_lookup(0);
CY0     = 0;
CYbeta  = aero_out(1);
CYp     = aero_out(2);
CYr     = aero_out(3);
CY_da    = aero_out(4);
% CY_de    = aero_out(5) + aero_out(6);
CY_dr    = aero_out(7);

Cl0     = 0;
Clbeta  = aero_out(8);
Clp     = aero_out(9);
Clr     = aero_out(10);
Cl_da   = aero_out(11);
Cl_de   = aero_out(12) + aero_out(13);
Cl_dr   = aero_out(14);

Cn0     = 0;
Cnbeta  = aero_out(15);
Cnp     = aero_out(16);
Cnr     = aero_out(17);
Cn_da   = aero_out(18);
Cn_de   = aero_out(19) + aero_out(20);
Cn_dr   = aero_out(21);

CD0     = aero_out(22);
CDq     = aero_out(23);
CD_de   = aero_out(24) + aero_out(25);
CDalpha = 0.12;

CL0     = aero_out(26);
CLq     = aero_out(27);
CL_de   = aero_out(28) + aero_out(29);
CLalpha = 4.0107;

Cm0     = aero_out(30);
Cmq     = aero_out(31);
Cm_de   = aero_out(32) + aero_out(33);
Cmalpha = -0.3292;

qbar = 1/2*rho*V^2;

CL = CL0 + CLalpha*alpha;
CD = CD0 + CDalpha*alpha;


CX    = -CD    * cos(alpha) + CL    * sin(alpha);
CXq   = -CDq   * cos(alpha) + CLq   * sin(alpha);
CX_de = -CD_de * cos(alpha) + CL_de * sin(alpha);

CZ    = -CD    * sin(alpha) - CL    * cos(alpha);
CZq   = -CDq   * sin(alpha) - CLq   * cos(alpha);
CZ_de = -CD_de * sin(alpha) - CL_de * cos(alpha);

Cp0    = Gamma3 * Cl0    + Gamma4 * Cn0;
Cpbeta = Gamma3 * Clbeta + Gamma4 * Cnbeta;
Cpp    = Gamma3 * Clp    + Gamma4 * Cnp;
Cpr    = Gamma3 * Clr    + Gamma4 * Cnr;
Cp_de  = Gamma3 * Cl_de  + Gamma4 * Cn_de;
Cp_da  = Gamma3 * Cl_da  + Gamma4 * Cn_da;
Cp_dr  = Gamma3 * Cl_dr  + Gamma4 * Cn_dr;

Cr0    = Gamma4 * Cl0    + Gamma8 * Cn0;
Crbeta = Gamma4 * Clbeta + Gamma8 * Cnbeta;
Crp    = Gamma4 * Clp    + Gamma8 * Cnp;
Crr    = Gamma4 * Clr    + Gamma8 * Cnr;
Cr_da  = Gamma4 * Cl_da  + Gamma8 * Cn_da;
Cr_dr  = Gamma4 * Cl_dr  + Gamma8 * Cn_dr;
Cr_de  = Gamma4 * Cl_de  + Gamma8 * Cn_de;







% Angular accelerations
f_phi = Gamma1 * p * q - Gamma2 * q * r + (qbar*S*b) * ...
       (Cp0 + Cpbeta * beta + Cpp * (b*p/2/V) + Cpr * (b*r/2/V)) ...
       + sin(phi)*tan(theta) * (Gamma5 * p * r - Gamma6 * (p^2 - r^2) + (qbar*S*c/Iyy) * (Cm0 + Cmalpha * alpha + Cmq * (c*q/2/V))) ...
       + cos(phi)*tan(theta) * (Gamma7 * p * q - Gamma1 * q * r + (qbar*S*b) * (Cr0 + Crbeta * beta + Crp * (b*p/2/V) + Crr * (b*r/2/V))) ...
       + q * (sin(phi) * sec(theta)^2 * (q * cos(phi) - r * sin(phi)) + cos(phi) * tan(theta) * (p + q * sin(phi) * tan(theta) + r * cos(phi) * tan(theta))) ...
       + r * (-sin(phi) * tan(theta) * (p + q * sin(phi) * tan(theta) + r * cos(phi) * tan(theta)) + cos(phi) * sec(theta)^2 * (q*cos(phi) - r*sin(phi)));



f_theta = cos(phi) * (Gamma5*p*r - Gamma6*(p^2 - r^2) + (qbar*S*c/Iyy) * (Cm0 + Cmalpha * alpha + Cmq * (c*q/2/V))) ...
        - sin(phi) * (Gamma7*p*q - Gamma1*q*r + (qbar*S*b) * (Cr0 + Crbeta * beta + Crp * (b*p/2/V) + Crr * (b*r/2/V))) ...
        + q * (-sin(phi) * (p + q * sin(phi) * tan(theta) + r * cos(phi) * tan(theta))) ...
        + r * (-cos(phi) * (p + q * sin(phi) * tan(theta) + r * cos(phi) * tan(theta)));


f_psi = 0;


psi_ddot = + cos(phi)*sec(theta) * (Gamma7*p*q - Gamma1*q*r + (qbar*S*b) * (Cr0 + Crbeta * beta + Crp * (b*p/2/V) + Crr * (b*r/2/V))) ...
           + sin(phi)*sec(theta) * (Gamma5*p*r - Gamma6*(p^2 - r^2) + (qbar*S*c/Iyy) * (Cm0 + Cmalpha * alpha + Cmq * (c*q/2/V))) ...
           - r * sin(phi)*sec(theta) * (p + q * sin(phi) * tan(theta) + r * cos(phi) * tan(theta)) ...
           + q * cos(phi)*sec(theta) * (p + q * sin(phi) * tan(theta) + r * cos(phi) * tan(theta)) ...
           + r * cos(phi)*tan(theta)*sec(theta) * (p + q * sin(phi) * tan(theta) + r * cos(phi) * tan(theta)) ...
           + q * sin(phi)*tan(theta)*sec(theta) * (p + q * sin(phi) * tan(theta) + r * cos(phi) * tan(theta));


% Translational acceleration
f_h = sin(theta) * (r*v - q*w - g*sin(theta) + (qbar*S/m) * (CX + CXq * (c*q/2/V))) ...
    - cos(theta) * sin(phi) * (p*w - r*u + g * cos(theta) * sin(phi) + (qbar*S/m) * (CY0 + CYbeta*beta + CYp*(b*p/2/V) + CYr*(b*r/2/V))) ...
    - cos(phi) * cos(theta) * (q*u - p*v + g * cos(theta)*cos(phi) + (qbar*S/m) * (CZ + CZq * (c*q/2/V))) ...
    + u * cos(theta) * (q*cos(phi)-r*sin(phi)) ...
    - v * ( cos(phi) * cos(theta) * (p + q * sin(phi) * tan(theta) + r * cos(phi) * tan(theta)) - sin(phi) * sin(theta) * (q * cos(phi) - r * sin(phi)) ) ...
    + w * ( sin(phi) * cos(theta) * (p + q * sin(phi) * tan(theta) + r * cos(phi) * tan(theta)) + cos(phi) * sin(theta) * (q * cos(phi) - r * sin(phi)) );





% Control derivatives
g_phi_de = sin(phi) * tan(theta) * (qbar*S*c/Iyy) * Cm_de;
g_phi_da = (qbar*S*b) * Cp_da + cos(phi) * tan(theta) * (qbar*S*b) * Cr_da;
g_phi_dr = (qbar*S*b) * Cp_dr + cos(phi) * tan(theta) * (qbar*S*b) * Cr_dr;
g_phi_dt = 0;

g_theta_de = cos(phi) * (qbar*S*c/Iyy) * Cm_de;
g_theta_da = -sin(phi) * (qbar*S*b) * Cr_da;
g_theta_dr = -sin(phi) * (qbar*S*b) * Cr_da;
g_theta_dt = 0;

g_psi_de = sin(phi) * sec(theta) * (qbar*S*c/Iyy) * Cm_de;
g_psi_da = cos(phi) * sec(theta) * (qbar*S*b) * Cr_da;
g_psi_dr = cos(phi) * sec(theta) * (qbar*S*b) * Cr_dr;
g_psi_dt = 0;

g_h_de = sin(theta) * (qbar*S/m) * CX_de - cos(phi) * cos(theta) * (qbar*S/m) * CZ_de;
g_h_da = -cos(theta) * sin(phi) * (qbar*S/m) * CY_da;
g_h_dr = -cos(theta) * sin(phi) * (qbar*S/m) * CY_dr;
g_h_dt = sin(theta) * Tmax / m;
%%
% f = [f_phi; f_theta; f_psi; f_h];

G = [g_phi_de    g_phi_da     g_phi_dr      g_phi_dt
     g_theta_de  g_theta_da   g_theta_dr    g_theta_dt
     g_psi_de    g_psi_da     g_psi_dr      g_psi_dt
     g_h_de      g_h_da       g_h_dr        g_h_dt];


s1 = phi_dot   - phi_dot_d   + ctrl_params.k1*(phi-phi_d);
s2 = theta_dot - theta_dot_d + ctrl_params.k2*(theta-theta_d);
s3 = psi_ddot  - psi_ddot_d  + ctrl_params.k3*(psi_dot-psi_dot_d);
s4 = h_dot     - h_dot_d     + ctrl_params.k4*(h-h_d);

temp = zeros(4,1);
temp(1) =    f_phi     + ctrl_params.k1*(phi_dot-phi_dot_d)         + ctrl_params.lambda1 *abs(s1)^ctrl_params.a1 * sign(s1); %*tanh(s1*ctrl_params.a1)   ;
temp(2) =    f_theta   + ctrl_params.k2*(theta_dot-theta_dot_d)     + ctrl_params.lambda2 *abs(s2)^ctrl_params.a2 * sign(s2); %*tanh(s2*ctrl_params.a2)   ;
temp(3) =    f_psi     + ctrl_params.k3*(psi_ddot-psi_ddot_d)       + ctrl_params.lambda3 *abs(s3)^ctrl_params.a3 * sign(s3); %*tanh(s3*ctrl_params.a3)   ;
temp(4) =    f_h       + ctrl_params.k4*(h_dot-h_dot_d)             + ctrl_params.lambda4 *abs(s4)^ctrl_params.a4 * sign(s4); %*tanh(s4*ctrl_params.a4)   ;

if abs(det(G))<1e-10
    u = zeros(4,1);
else
    u = -G^-1 * temp;
end

u_out = [  min(max(-deg2rad(25),u(1)),deg2rad(10));
           min(max(-deg2rad(35),u(2)),deg2rad(35));
           min(max(-deg2rad(30),u(3)),deg2rad(30))*0.1+0.9*input_store(3)
           min(max(0.15,u(4)),1)                    ];

s_out = [s1;s2;s3;s4];
end