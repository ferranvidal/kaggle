clc; clear all; close all;


outfolder = [pwd,'/Simulation_data/'];
resume = 0;
porder = 3;
nodetype = 1;
fac_beta2 = 1;
nodesPar = 0;
M = 12;
typ = 'Au';
suff = '_new';


%% Generate mesh
D0 = .25;
G0 = 0.002;
H3 = 0.15; %gold thickness

% non plot
n1 = 5; n2 = 8; n4 = 8;
n3 = 3; 

if strcmp(typ,'Au')
    d = 1e-4*[1 2.5 4.5 7];
    a2 = 2; a1 = -7; a3 = 7;
    
elseif strcmp(typ,'Ti')
    d = 1e-4*[1 2.5 4.5 7];
    d = 3.88*d;
    a2 = 2; a1 = -5; a3 = 5;
end


if G0 == 0.01
    
    
    lambda_vec = [2.78;2.91];
    
    G = 0.0095;
    R = 0.253/2;
elseif G0 == 0.007
    
    
    lambda_vec = [2.95;3.35];
    
    G = 0.0071;
    R = 0.251/2;
elseif G0 == 0.005
    
    
    lambda_vec = [3.55;3.66];
    G = 0.0049;
    R = 0.255/2;
    
elseif G0 == 0.004
    
    lambda_vec = [3.5;4.2];
    
    G = 0.004;
    R = 0.254/2;
elseif G0 == 0.003
    
    
    lambda_vec = [6;7.5];
    
    G = 0.0031;
    R = 0.252/2;
elseif G0 == 0.002
    
    lambda_vec = [4.5;5.3];
    G = 0.0019;
    R = 0.254/2;
    
elseif G0 == 0.001
    
    lambda_vec = [8.1;11];
    
    G = 0.0011;
    R = 0.252/2;
    
end

% lambda_vec = 4.35;

D  = D0;% Reference lengthscale is 1*Lo micron
H1 = .35; %glass
H2 = G;
H4 = .4; %air
H0 = .05;
H5 = .05;
m0 = 5;
m5 = 5;
m1=6; m2 = 3; m3 = 9; m4=6;

b0 = 0;
b2 = 0;
b1 = -4; b3 = 6; b4 = 4;

R0 = R/2; R1 = R+G+R0;
air_indx = 3; metal_indx = 2; glass_indx = 1; alH_indx = 4; alV_indx = 5;

params2d_space =[n1 n2 n3 n4 a1 a2 a3];
params2d_length =[R0 R R+G R1 D d];
dz = H3*[0.005 0.015];

params3d = {{[H1-H0,H0], H2,H3-H2, [H5,H4-H5]},{[m0,m1],m2,m3, [m4,m5]},{[b0,b1], 0,b3, [b4,b0]},dz};
mesh = mkmesh_annular3d_MIR(porder(1),params2d_length,params2d_space,params3d,[glass_indx metal_indx air_indx alH_indx alV_indx],1,nodetype);

z = unique(mesh.p(:,3));
mesh.zfilm = H1+[0,H3];
mesh.zgoldmid = H1 + H2 + (H3-H2)/2;
fieldEnh_vol = @(p) all(p(:,3)>H1-1e-4 & p(:,3)<H1+H3+1e-4 & sqrt(p(:,1).^2+p(:,2).^2)>=R-1e-6 & sqrt(p(:,1).^2+p(:,2).^2)<=R+G+1e-6);
mesh.fieldEnh_vol = fieldEnh_vol;
fieldEnh = @(p) all(p(:,3)>H1+H2-1e-4 & p(:,3)<H1+H3+1e-4 & p(:,1) >R-1e-6 & p(:,1)<R+G+1e-6 & abs(p(:,2)) < 1e-5);
mesh.fieldEnh = fieldEnh;

power1 = @(p) all(abs(p(:,3)-z(end-(m5-1)))<1e-4);
power2 = @(p) all(abs(p(:,3)-z(end-(m5-1)+1))<1e-4);
power3 = @(p) all(abs(p(:,3)-z(end-(m5-1)+2))<1e-4);
power0 = @(p) all(abs(p(:,3)-z(m0+1))<1e-4);
mesh.power0 = power0;
mesh.power1 = power1;
mesh.power2 = power2;
mesh.power3 = power3;

npv = mesh.npv;
nd = 3;
ne = mesh.ne;

mesh.permeability = ones(npv,ne); % permeability mu
mesh.permittivity = ones(npv,ne,3); % permittivity epsilon


% Units (micron scale)
if strcmp(typ,'Au')
%     omegapE = 8.45; %eV
%     gammaE = 0.047; %eV
%     einf = 9.84;
    omegapE = 8.84; %eV
    gammaE = 0.103; %eV
    einf = 9.84;
    vf = 1.39e6; %m/s
    
elseif strcmp(typ,'Ti')
    omegapE = 2.8; %eV
    gammaE = 0.082; %eV
    einf= 2.2;
    vf = 1.79e6; %m/s
end

c = 3e8; %m/s
Lc = 1e-6; % m
Tc = 1e12; % Hz or 1/s
hinv = 1/6.582119514*1e16; %1/(eV.s)
transform = hinv*Lc/c;
omegap = omegapE*hinv*Lc/c; % adimensional
gamma = gammaE*hinv*Lc/c; % adimensional
omega0 = 2*pi*Lc*Tc/c; % adimensional
beta2 = fac_beta2*3/5*(vf/c)^2;

data.source = 0;
data.bcm = [1;1;2;2;3;3];
data.fbou = @fbou_3d;
data.k = [0 0 1];
data.p = [1 0 0];
data.bcm_n = [1;1;2;2;2;2];
data.omega0 = omega0;
data.tau = 2*omega0*300/min(lambda_vec);
data.tau_n = 1*omegap/sqrt(beta2);

%% RB Parameters

if M == 1
    lambda = lambda_vec;
else
    lambda = linspace(lambda_vec(1),lambda_vec(2),M)';
end

f = 300./lambda;

%% Materials

layer_index = [1 2 3];
metal_index = 2;
gap_index = [4 5];
anisotropy = [1 1 1 1 1];

epsilon = zeros(M,sum(anisotropy));
epsilon(:,1) = sapphire(lambda);
epsilon(:,3) = 1;
eal = alumina2(lambda);
epsilon(:,4) = eal;
epsilon(:,5) = eal;

if fac_beta2 == 0
    epsilon(:,2) = drudemodel(omega0*f,omegap,gamma,einf);
else
    epsilon(:,2) = einf;
end
mat_index.layer_index = layer_index;
mat_index.gap_index = gap_index;
mat_index.metal_index = metal_index;
mat_index.anisotropy = anisotropy;
mesh.materialIndex = mat_index;
num_mat = 5;

%% Outputs
FE_v = zeros(M,1);
FE_f = zeros(M,1);

Pt = zeros(M,1);
Pr = zeros(M,1);
P0 = zeros(M,1);
fields = cell(M,1);
[~,A] = transmission(mesh,zeros(npv,nd,ne),zeros(npv,nd,ne),power0);

k = data.k;
p = data.p;
kxp = cross(k,p);
kX = mesh.dgnodes(:,1,:)*k(1) + ...
    mesh.dgnodes(:,2,:)*k(2) + ...
    mesh.dgnodes(:,3,:)*k(3);

if k(3) > 0
    id_incident = layer_index(1);
else
    id_incident = layer_index(end);
end
n_incident = sqrt(epsilon(:,id_incident));

decay_mid = @(p) all(abs(p(:,2))<1e-5 &  abs(p(:,3) - H1-H2-(H3-H2)/2)< 1e-5);
psurf = decay_plasmon3d(mesh,zeros(mesh.npv,3,mesh.ne),zeros(mesh.npv,3,mesh.ne),decay_mid);
% num_psurf = length(psurf1);

%% Generate model
id = ['G=',num2str(G0*1000),'_T=',num2str(1000*H3),'_b=',num2str(fac_beta2),'_',typ,suff];

fprintf('G %i / R %.1f / T %i / b2 %.2f \n',G0*1000,R*1000,1000*H3,fac_beta2)

if resume
    load([outfolder,num2str(resume),'_',id,'.mat']);
    mm = (resume+1):M;
else
    mm = 1:M;
end


if fac_beta2 > 0
    nodesPar = 0;
else
    parpool('local',nodesPar)
end

parfor (m = mm, nodesPar)
    
    m
    fm = f(m);
    omega = fm*omega0;
    epsilonm = epsilon(m,:);
    mesh_m = mesh;
    
    k = 1;
    for n = 1:num_mat
        tm = mesh_m.mate == n;
        if anisotropy(n) == 1
            mesh_m.permittivity(:,tm,:) = epsilonm(k);
            k = k + 1;
        elseif anisotropy(n) == 3
            mesh_m.permittivity(:,tm,1) = epsilonm(k);
            mesh_m.permittivity(:,tm,2) = epsilonm(k+1);
            mesh_m.permittivity(:,tm,3) = epsilonm(k+2);
            k = k + 3;
        else
            error('Specify 1 component (isotropic) or 3 components (3D anisotropic)' )
        end
    end
    
    
    if fac_beta2 == 0
        [EDG,VDG,eh] = hdg_maxwell_a(mesh_m,data,omega,0);
    else
        param = {omega,omegap,gamma,beta2,1};
        [EDG,VDG,JDG,UDG,eh] = hdg_maxwellhydro_a(mesh_m,data,param,0);
    
    end
    HDG = -1i/omega*VDG;
    
    omegaeff = 1i*omega0*f(m)*n_incident(m);
    
    Ei = zeros(npv,nd,ne);
    Hi = zeros(npv,nd,ne);
    Ei(:,logical(p),:) = exp(omegaeff*kX);
    Hi(:,logical(kxp),:) = n_incident(m)*exp(omegaeff*kX);
    
    [F,V] = field_enhancement_vol(mesh_m,EDG,fieldEnh_vol);
    FE_v(m) = F/V;
    [F,V] = field_enhancement(mesh_m,EDG,fieldEnh);
    FE_f(m) = F/V;
    
    Pt(m,:) = 0.5*transmission(mesh_m,EDG,HDG,power3)/A;
    Pr(m) = 0.5*transmission(mesh_m,EDG-Ei,HDG-Hi,power0)/A;
    P0(m) = 0.5*transmission(mesh_m,Ei,Hi,power0)/A;
    
    
    SDG = cross(EDG,conj(HDG),2);
    [~,tm1] = decay_plasmon3d(mesh_m,EDG,HDG,decay_mid);
    fields{m} = struct('EDGx',EDG(:,1,:),'SDGz',SDG(:,3,:),'decay',tm1);
    
%     if m < M
%         save([outfolder,num2str(m),'_',id,'.mat'],'data','mesh','M','D0','lambda_vec','f','lambda','FE','Pt','Pr','P0','fields','psurf');
%     end
%     if m > 1
%         delete([outfolder,num2str(m-1),'_',id,'.mat'])
%     end
%     
    fprintf('G %i / R %.1f / T %i / b2 %.2f \n',G0*1000,R*1000,1000*H3,fac_beta2)
    m
end
save([outfolder,id,'.mat'],'data','mesh','M','D0','lambda_vec','f','lambda','FE_v','FE_f','Pt','Pr','P0','fields','psurf','epsilon');
fprintf('G %i / R %.1f / T %i\n',G0*1000,R*1000,1000*H3)
