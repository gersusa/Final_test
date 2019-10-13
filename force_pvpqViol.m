function [mpcPF] = force_pvpqViol(mpcPF,mpcLim)
%force_PVPQVIOL 
%This function indetyfy 4 differents types of QV violations 
%Type1_viol : Generators with reactive power higher than a maxLim and voltage near to Setpoint, the solution is leave the same voltage and Qmax
%Type2_viol : Generators with reactive power near to max y voltage higher than setPoint, the solution is leave is the buses the prefault voltage 
%Type3_viol : Generators with reactive power lower than a minLim and voltage near to Setpoint, the solution is leave the same voltage and Qmin
%Type4_viol : Generators with reactive power near to min y voltage lower than setPoint, the solution is leave is the buses the prefault voltage 
%mpcPF is the analysis case, it must had been  edited with forceVoltage and checkRigidConstraints
%mpcLim has the reactive limits


%% MatPower Constants
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
%% PV/PQ violation evaluation

% Bus indices of busses with generation
[~,bgOn]= ismember(mpcPF.gen(:,GEN_BUS),mpcPF.bus(:,BUS_I));
TOL=1e-4;
%Save initial conditions
Vset=mpcPF.gen(:,VG);
Vnow=mpcPF.bus(bgOn,VM);
Qgnow=mpcPF.gen(:,QG);
Qmax=mpcLim.gen(:,QMAX);
Qmin=mpcLim.gen(:,QMIN);

%CHECK AND CORRECT VIOLATIONS
% Type1_viol= mpcPF.gen(:,GEN_STATUS)~=0 & (abs(Vnow-Vset)<TOL) & ((Qgnow-Qmax)>TOL);
% mpcPF.gen(Type1_viol,QG)=mpcPF.gen(Type1_viol,QMAX);


Type2_viol=mpcPF.gen(:,GEN_STATUS)~=0 & (abs(Qgnow-Qmax)<TOL) & (Vnow-Vset>TOL);
mpcPF.bus(ismember(mpcPF.bus(:,1),mpcPF.gen(Type2_viol,GEN_BUS)),VM)=unique(mpcPF.gen(Type2_viol,VG),'stable');


% Type3_viol=mpcPF.gen(:,GEN_STATUS)~=0 & (abs(Vset-Vnow)<TOL) & (Qgnow-Qmin<-TOL);
% mpcPF.gen(Type3_viol,QG)=mpcPF.gen(Type3_viol,QMIN);


Type4_viol=mpcPF.gen(:,GEN_STATUS)~=0 & (abs(Qgnow-Qmin)<TOL) & (Vnow-Vset<-TOL);
mpcPF.bus(ismember(mpcPF.bus(:,1),mpcPF.gen(Type4_viol,GEN_BUS)),VM)=unique(mpcPF.gen(Type4_viol,VG),'stable');



end

