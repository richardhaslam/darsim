function [J_c, Residual_c] = Restrict(J, Residual, R, Pp, Ps, N, maxLevel)
% Restrict system to DLGR Grid
Jop = J(1:N, 1:N);
Jos = J(1:N, N+1:2*N);
Jwp = J(N+1:2*N, 1:N);
Jws = J(N+1:2*N, N+1:2*N);
Jopc = R(1).matrix*Jop*Pp(1).matrix;
Josc = R(1).matrix*Jos*Ps(1).matrix;
Jwpc = R(1).matrix*Jwp*Pp(1).matrix;
Jwsc = R(1).matrix*Jws*Ps(1).matrix;
Residual_cp = R(1).matrix*Residual(1:N);
Residual_cs = R(1).matrix*Residual(N+1:2*N);
for i=2:maxLevel
    Jopc = R(i).matrix*Jopc*Pp(i).matrix;
    Josc = R(i).matrix*Josc*Ps(i).matrix;
    Jwpc = R(i).matrix*Jwpc*Pp(i).matrix;
    Jwsc = R(i).matrix*Jwsc*Ps(i).matrix;
    Residual_cp = R(i).matrix*Residual_cp;
    Residual_cs = R(i).matrix*Residual_cs;
end
J_c = [Jopc, Josc; Jwpc, Jwsc];
Residual_c = [Residual_cp; Residual_cs];
end