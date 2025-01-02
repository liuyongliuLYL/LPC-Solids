cells=1000;
den = fread(fopen("./Out/density_all.dat", 'rb'),'double');den = reshape(den,cells,[]);
p = fread(fopen("./Out/pressure_all.dat", 'rb'),'double');p = reshape(p,cells,[]);
u = fread(fopen("./Out/velocity_all.dat", 'rb'),'double');u = reshape(u,cells,[]);
sigma = fread(fopen("./Out/sigma_x_all.dat", 'rb'),'double');sigma = reshape(sigma,cells,[]);
%s_x = fread(fopen("./Out/s_x_all.dat", 'rb'),'double');s_x = reshape(s_x,cells,[]);
%s_y = fread(fopen("./Out/s_y_all.dat", 'rb'),'double');s_y = reshape(s_y,cells,[]);
e = fread(fopen("./Out/i_energy_all.dat", 'rb'),'double');e = reshape(e,cells,[]);
% E = fread(fopen("./Out/t_energy_all.dat", 'rb'),'double');E = reshape(E,cells,[]);
X = fread(fopen("./Out/xpoint_all.dat", 'rb'),'double');X = reshape(X,cells,[]);
X_La = X;   X_La(:, 2:end) = repmat(X_La(:, 1), 1, size(X_La, 2)-1);
T = fread(fopen("./Out/t_all.dat", 'rb'),'double');T = reshape(T,cells,[]);
%ce = fread(fopen("./Out/ce.dat", 'rb'),'double');ce = reshape(ce,cells,[]);

r_1 = readBinary("./Out/r_1.dat");r_1 = reshape(r_1,cells,[]);
p_1 = -readBinary("./Out/p_1.dat");p_1 = reshape(p_1,cells,[]);
sigma11_1 = -readBinary("./Out/sigma11_1.dat");sigma11_1 = reshape(sigma11_1,cells,[]);
sigma22_1 = -readBinary("./Out/sigma22_1.dat");sigma22_1 = reshape(sigma22_1,cells,[]);
big_omega_1 = readBinary("./Out/big_omega_1.dat");big_omega_1 = reshape(big_omega_1,cells,[]);
big_phi_1 = readBinary("./Out/big_phi_1.dat");big_phi_1 = reshape(big_phi_1,cells,[]);
u_1 = fread(fopen("./Out/velocity_1.dat", 'rb'),'double');u_1 = reshape(u_1,cells,[]);
v_1 = fread(fopen("./Out/volume_1.dat", 'rb'),'double');v_1 = reshape(v_1,cells,[]);


figure;contourf(real(X_La),real(T),real(den),50);c = colorbar;c.Label.String = 'density';

Ri = 500; Ri = 1;Ri = floor(cells/2);
figure; plot(T(Ri,:),(r_1(Ri,:)),'-','LineWidth', 1);
ax = gca;ax.YScale = 'log';

