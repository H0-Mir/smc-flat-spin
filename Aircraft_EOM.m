function states_dot = Aircraft_EOM(~,params,states,input)

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


V = sqrt(u.^2 +v.^2 +w.^2);

alpha = atan(w/u);
beta  = asin(v/V);

delta_e = input(1);
delta_a = input(2);
delta_r = input(3);
delta_t = input(4);

%% Aerodynamic Coefficients
aero_out = aero_lookup(0);
CY0     = 0;
CYbeta  = aero_out(1);
CYp     = aero_out(2);
CYr     = aero_out(3);
CY_da   = aero_out(4);
CY_de   = aero_out(5) + aero_out(6);
CY_dr   = aero_out(7);

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

CL = CL0 + CLalpha*alpha ;
CD = CD0 + CDalpha*alpha ;
% Cm = Cm0 + Cmq*q*(c/2/V) + Cm_de*delta_e;

% D = qbar*S*CD;
% L = qbar*S*CL;

CX    = -CD * cos(alpha) + CL * sin(alpha);
CXq   = -CDq * cos(alpha) + CLq * sin(alpha);
CX_de = -CD_de * cos(alpha) + CL_de * sin(alpha);

CZ    = -CD * sin(alpha) - CL * cos(alpha);
CZq   = -CDq * sin(alpha) - CLq * cos(alpha);
CZ_de = -CD_de * sin(alpha) - CL_de * cos(alpha);

Cp0    = Gamma3 * Cl0 + Gamma4 * Cn0;
Cpbeta = Gamma3 * Clbeta + Gamma4 * Cnbeta;
Cpp    = Gamma3 * Clp + Gamma4 * Cnp;
Cpr    = Gamma3 * Clr + Gamma4 * Cnr;
Cp_de  = Gamma3 * Cl_de + Gamma4 * Cn_de;
Cp_da  = Gamma3 * Cl_da + Gamma4 * Cn_da;
Cp_dr  = Gamma3 * Cl_dr + Gamma4 * Cn_dr;

Cr0    = Gamma4 * Cl0    + Gamma8 * Cn0;
Crbeta = Gamma4 * Clbeta + Gamma8 * Cnbeta;
Crp    = Gamma4 * Clp    + Gamma8 * Cnp;
Crr    = Gamma4 * Clr    + Gamma8 * Cnr;
Cr_da  = Gamma4 * Cl_da  + Gamma8 * Cn_da;
Cr_dr  = Gamma4 * Cl_dr  + Gamma8 * Cn_dr;
Cr_de  = Gamma4 * Cl_de  + Gamma8 * Cn_de;


%% EOM

% Rotation matrix from body frame to inertial frame
C = [cos(theta)*cos(psi), sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi), cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi);
     cos(theta)*sin(psi), sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi), cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi);
     -sin(theta),         sin(phi)*cos(theta),                            cos(phi)*cos(theta)];

% Rotation matrix from body angular rates to euler angles rates
R = [1,  sin(phi)*tan(theta),  cos(phi)*tan(theta);
     0,  cos(phi),            -sin(phi);
     0,  sin(phi)/cos(theta),  cos(phi)/cos(theta)];


phi_dot   = R(1,:)*[p; q; r];
theta_dot = R(2,:)*[p; q; r];
psi_dot   = R(3,:)*[p; q; r];

% Translational
u_dot = r*v - q*w - g*sin(theta) + (qbar*S/m) * ...
        (CX + CXq * (c*q/2/V) + CX_de * delta_e) + (Tmax/m) * delta_t;



v_dot = p*w - r*u + g*cos(theta)*sin(phi) + (qbar*S/m) * ...
        (CY0 + CYbeta * beta + CYp * (b*p/2/V) + CYr * (b*r/2/V) + ...
         CY_de * delta_e + CY_da * delta_a + CY_dr * delta_r);

w_dot = q*u - p*v + g*cos(theta)*cos(phi) + (qbar*S/m) * ...
        (CZ + CZq * (c*q/2/V) + CZ_de * delta_e);



% Compute position derivatives in inertial frame
X_dot = C(1,:) * [u; v; w];
Y_dot = C(2,:) * [u; v; w];
h_dot = -C(3,:) * [u; v; w];


% Rotational
p_dot = Gamma1 * p * q - Gamma2 * q * r + (qbar*S*b) * ...
        (Cp0 + Cpbeta*beta +  Cpp * (b*p/2/V) + Cpr * (b*r/2/V) + ...
         Cp_de * delta_e + Cp_da * delta_a + Cp_dr * delta_r);

q_dot = Gamma5 * p * r - Gamma6 * (p^2 - r^2) + (qbar*S*c/Iyy)*...
                  ( Cm0 + Cmalpha*alpha + Cmq*(c*q/2/V) + Cm_de*delta_e );
r_dot = Gamma7 * p * q - Gamma1 * q * r + (qbar*S*b) * ...
        (Cr0 + Crbeta * beta +  Crp * (b*p/2/V) + Crr * (b*r/2/V) + ...
         Cr_de * delta_e + Cr_da * delta_a + Cr_dr * delta_r);


states_dot = [ u_dot;     v_dot;      w_dot;
               p_dot;     q_dot;      r_dot;
               phi_dot;   theta_dot;  psi_dot;
               X_dot;     Y_dot;      h_dot    ];

end