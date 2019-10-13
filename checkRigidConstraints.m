function [ mpcPFout ] = checkRigidConstraints( mpcPF, forceV )
%CHECKRIGIDCONSTRAINTS Check and force rigid constraints
%   This functions verifies and forces rigid constraints in a MatPower
%   case. It forces generator reactive values within limits, and changes
%   scheduled voltages of generators if necessary to make sure that PV/PQ
%   definitions are complied with.
%   
%   Copyright (c) 2019, Gers USA
%   by Tomas Valencia tvalencia@gersusa.com
%% MatPower Constants
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
%% Initialize
mpcPFout = mpcPF;
%% Verify rigid constraints
if mpcPFout.success || forceV
    % Force Q_G within [QMIN, QMAX]
    mpcPFout.gen(:,QG) = min(max(mpcPFout.gen(:,QG),mpcPFout.gen(:,QMIN)),...
        mpcPFout.gen(:,QMAX));
    mpcPFout.gen(mpcPFout.gen(:,GEN_STATUS)==0,QG) = 0;
    % Make sure that PV/PQ limits are respected.
    % In other words, if Q @ limit, VG must be unreached.
    g_qMin = mpcPFout.gen(:,QG)==mpcPFout.gen(:,QMIN);
    [~,idb_g_qMin] = ismember(mpcPFout.gen(g_qMin,GEN_BUS),...
        mpcPFout.bus(:,BUS_I));
    mpcPFout.gen(g_qMin,VG)=min(mpcPFout.gen(g_qMin,VG),...
        mpcPFout.bus(idb_g_qMin,VM));
    g_qMax = mpcPFout.gen(:,QG)==mpcPFout.gen(:,QMAX);
    [~,idb_g_qMax] = ismember(mpcPFout.gen(g_qMax,GEN_BUS),...
        mpcPFout.bus(:,BUS_I));
    mpcPFout.gen(g_qMax,VG)=max(mpcPFout.gen(g_qMax,VG),...
        mpcPFout.bus(idb_g_qMax,VM));
    
end

end

