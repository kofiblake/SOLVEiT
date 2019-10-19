
% dist = load('dist.dat');
% pos1 = load('pos1.dat');
% pos2 = load('pos2.dat');
% pos3 = load('pos3.dat');
% q_m_frag1 = load('q_m_frag1.dat');
% q_m_frag2 = load('q_m_frag2.dat');
% q_m_frag3 = load('q_m_frag3.dat');
% rho = load('rho.dat');
% vel1 = load('vel1.dat');
% vel2 = load('vel2.dat');
% vel3 = load('vel3.dat');

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

pt_size = 2;
figure
hold on
scatter(pos1(:,3), pos1(:,1), pt_size, 'filled')
scatter(pos2(:,3), pos2(:,1), pt_size, 'filled')
scatter(pos3(:,3), pos3(:,1), pt_size, 'filled')
xlabel('z (m)'); ylabel('x (m)'); title('Beam Divergence in the x-z plane')
legend('Monomers', 'Dimers', 'Neutrals')
hold off

figure
hold on
scatter3(pos1(:,1),pos1(:,2),pos1(:,3), pt_size, 'filled')
scatter3(pos2(:,1),pos2(:,2),pos2(:,3), pt_size, 'filled')
scatter3(pos3(:,1),pos3(:,2),pos3(:,3), pt_size, 'filled')
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)')
title('Beam Divergence')
hold off


