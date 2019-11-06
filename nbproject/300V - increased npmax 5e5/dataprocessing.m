%Loading in the Data (Only once)

if exist('auxnp', 'var') == 0
    auxnp = load('auxnp.dat'); end
if exist('dist', 'var') == 0
    dist = load('dist.dat'); end
if exist('pos1', 'var') == 0
    pos1 = load('pos1.dat'); end
if exist('pos2', 'var') == 0
    pos2 = load('pos2.dat'); end
if exist('pos3', 'var') == 0
    pos3 = load('pos3.dat'); end
if exist('q_m_frag1', 'var') == 0
    q_m_frag1 = load('q_m_frag1.dat'); end
if exist('q_m_frag2', 'var') == 0
    q_m_frag2 = load('q_m_frag2.dat'); end
if exist('q_m_frag3', 'var') == 0
    q_m_frag3 = load('q_m_frag3.dat'); end
if exist('rho', 'var') == 0
    rho = load('rho.dat'); end
if exist('vel1', 'var') == 0
    vel1 = load('vel1.dat'); end
if exist('vel2', 'var') == 0
    vel2 = load('vel2.dat'); end 
if exist('vel3', 'var') == 0
    vel3 = load('vel3.dat'); end

%Definining Constants

e_0 = 8.854e-12; % Permittivity of Free Space [Farads / meter]
k_b = 1.3806e-23; % Boltzmann Constant [J / K]
q = -1.602e-19; % Charge of an Electron [Coulomb]

%Domain Size Specifications - Make sure this aligns with values in
%   stats.dat!

dgz =   2E-04; %[m]
dgxy =   1E-05; %[m]
X = 75;
Y = 75;
Z = 24;
max_ion_angle = 5; %[degree]

dsizez =   dgz * Z; %[m]
dsizexy =  dgxy * X; %[m]

%2D Plot of Beam Divergence

pt_size = 2;
figure
hold on
scatter(pos1(:,3)*1e3, pos1(:,1)*1e3, pt_size+2, 'filled')
scatter(pos2(:,3)*1e3, pos2(:,1)*1e3, pt_size, 'filled')
scatter(pos3(:,3)*1e3, pos3(:,1)*1e3, pt_size, 'filled')
xlabel('z (mm)'); ylabel('x (mm)'); title('Beam Divergence in the x-z plane')

legend('Monomers', 'Dimers', 'Neutrals')
hold off

%3D Plot of Divergence 

figure
hold on
scatter3(pos1(:,1),pos1(:,2),pos1(:,3), pt_size, 'filled')
scatter3(pos2(:,1),pos2(:,2),pos2(:,3), pt_size, 'filled')
scatter3(pos3(:,1),pos3(:,2),pos3(:,3), pt_size, 'filled')
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)')
title('Beam Divergence')
hold off

%Num

z_grid = linspace(0, dsizez, 1001);
grid_bin_count = zeros(length(z_grid)-1, 1);
bin_volume = zeros(length(grid_bin_count), 1);
data_bins = discretize(pos2(:,3), z_grid);

for n = 1 : length(data_bins)
    current_bin = data_bins(n);
    grid_bin_count(current_bin) = grid_bin_count(current_bin) + 1;
end

for n = 1 : length(bin_volume)
    bin_volume(n) = pi / 3 * z_grid(n+1)^3*tand(max_ion_angle)^2;
    
    if n > 1
        a = n-1;
        for a_ind = a:-1:1
            bin_volume(n) = bin_volume(n) - bin_volume(a_ind);
        end
    end
end

begin_truncator = 50;
bin_density = grid_bin_count./bin_volume;

figure
plot(z_grid(1, begin_truncator:length(z_grid)), bin_density(begin_truncator-1:length(bin_density), 1));

%Debye Length Calculations
eV_to_K = 11604.518;
electron_temp = .025 * eV_to_K;
debye_length = sqrt((e_0*k_b*electron_temp)./(bin_density .* q^2));

figure
plot(z_grid(1, begin_truncator:length(z_grid)), debye_length(begin_truncator-1:length(debye_length), 1));


