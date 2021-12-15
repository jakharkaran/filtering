function [S11_residual, S12_residual, S22_residual, ...
    S11_leonard, S12_leonard, S22_leonard, ...
    S11_cross, S12_cross, S22_cross,...
    S11_reynolds, S12_reynolds, S22_reynolds,...
    S11_reynolds_1, S12_reynolds_1, S22_reynolds_1,...
    S11_reynolds_2, S12_reynolds_2, S22_reynolds_2] = residualStressComponents2D( ...
    U_DNS,V_DNS, filterType,coarseGrainingType, Delta, N_LES)

%Components of Residual Stress Term for 2D Turbulence
% U_DNS = Uf (filtered/resolved) + Ud (residual)
% Uf_c filtered coarse grained
% Udf_c residual filtered coarse grained

Uf = filter2D(U_DNS,filterType,'none', Delta, N_LES);
Vf = filter2D(V_DNS,filterType,'none', Delta, N_LES);

Ud = U_DNS - Uf;
Vd = V_DNS - Vf;

Udf_c = filter2D(Ud,filterType,coarseGrainingType, Delta, N_LES);
Vdf_c = filter2D(Vd,filterType,coarseGrainingType, Delta, N_LES);

UdUdf_c = filter2D(Ud.*Ud,filterType,coarseGrainingType, Delta, N_LES);
VdVdf_c = filter2D(Vd.*Vd,filterType,coarseGrainingType, Delta, N_LES);
UdVdf_c = filter2D(Ud.*Vd,filterType,coarseGrainingType, Delta, N_LES);
         
Uf_c = filter2D(U_DNS,filterType,coarseGrainingType, Delta, N_LES);
Vf_c = filter2D(V_DNS,filterType,coarseGrainingType, Delta, N_LES);

Uff_c = filter2D(Uf,filterType,coarseGrainingType, Delta, N_LES);
Vff_c = filter2D(Vf,filterType,coarseGrainingType, Delta, N_LES);

UUf_c = filter2D(U_DNS.*U_DNS,filterType,coarseGrainingType, Delta, N_LES);
VVf_c = filter2D(V_DNS.*V_DNS,filterType,coarseGrainingType, Delta, N_LES);
UVf_c = filter2D(U_DNS.*V_DNS,filterType,coarseGrainingType, Delta, N_LES);

UfUf_f_c = filter2D(Uf.*Uf,filterType,coarseGrainingType, Delta, N_LES);
VfVf_f_c = filter2D(Vf.*Vf,filterType,coarseGrainingType, Delta, N_LES);
UfVf_f_c = filter2D(Uf.*Vf,filterType,coarseGrainingType, Delta, N_LES);

UfUd_f_c = filter2D(Uf.*Ud,filterType,coarseGrainingType, Delta, N_LES);
VfVd_f_c = filter2D(Vf.*Vd,filterType,coarseGrainingType, Delta, N_LES);
UfVd_f_c = filter2D(Uf.*Vd,filterType,coarseGrainingType, Delta, N_LES);
VfUd_f_c = filter2D(Vf.*Ud,filterType,coarseGrainingType, Delta, N_LES);

S11_residual = UUf_c - Uf_c.*Uf_c;
S12_residual = UVf_c - Uf_c.*Vf_c;
S22_residual = VVf_c - Vf_c.*Vf_c;

S11_leonard = UfUf_f_c - Uff_c.*Uff_c;
S12_leonard = UfVf_f_c - Uff_c.*Vff_c;
S22_leonard = VfVf_f_c - Vff_c.*Vff_c;

S11_cross = 2* (UfUd_f_c - Uff_c.*Udf_c);
S12_cross = UfVd_f_c + VfUd_f_c - Uff_c.*Vdf_c - Vff_c.*Udf_c;
S22_cross = 2* (VfVd_f_c - Vff_c.*Vdf_c);

S11_reynolds = UdUdf_c - Udf_c.*Udf_c;
S12_reynolds = UdVdf_c - Udf_c.*Vdf_c;
S22_reynolds = VdVdf_c - Vdf_c.*Vdf_c;

S11_reynolds_1 = UdUdf_c;
S12_reynolds_1 = UdVdf_c ;
S22_reynolds_1 = VdVdf_c ;

S11_reynolds_2 = Udf_c.*Udf_c;
S12_reynolds_2 = Udf_c.*Vdf_c;
S22_reynolds_2 = Vdf_c.*Vdf_c;

end