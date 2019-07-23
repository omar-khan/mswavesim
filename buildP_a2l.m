function [P_a2l_x, P_a2l_z, P_l2a_x, P_l2a_z] = buildP_a2l(m0)
    [n ndel]=size(m0);
    aus(1:n)=1;
    aus=aus';

    P_a2l_x=cross(m0,kron((aus),[ 0 0 1]));
    P_a2l_x=P_a2l_x./kron(sqrt(sum(P_a2l_x'.^2)'),[1 1 1]);
    P_a2l_z=cross(P_a2l_x,m0);

    P_l2a_x = P_a2l_x';
    P_l2a_z = P_a2l_z';
end