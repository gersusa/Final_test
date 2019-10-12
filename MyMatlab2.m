function MyMatlab2(InFile1, InFile2, InFile3, InFile4, TimeLimitInSeconds, ScoringMethod, NetworkModel)
%MYMATLAB2 Version1 Script. This script generates solution2.txt
%
% Input Arguments:
% TimeLimitInSeconds and ScoringMethod are not used yet, the algorithm acts
% as it would for ScoringMethod = 0.
% InFile1: CON file
% InFile2: INL file
% InFile3: RAW file
% InFile4: ROP file
% pfs is in external notation

[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

%[mpc,contingencies] = convert2mpc(InFile3,InFile4,InFile2,InFile1);
load('mpc.mat');
% Get switched shunts data

[~,pfs] = runAllCONS(mpcOPF, contingencies,mpcOPF_or, 'AC',1); %Force voltages activate


order_save = mpcOPF_or.order.bus.e2i; % busses w/ fixed shunts

%%
a=tic;
if ~isempty(pfs)
    mpcOPF2 = mpcOPF;
    mpcOPF_or2 = mpcOPF_or;
    [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
    parfor i=1:numel(pfs)
        pfs{i}.order.bus.e2i= order_save;
        pfs{i}.indexMap=mpcOPF2.indexMap;
        pfs{i}.comp=mpcOPF_or2.comp;
        pfs{i}.bus=[pfs{i}.bus mpcOPF_or2.bus(:,PD) mpcOPF_or2.bus(:,QD)];    %BE CAREFULL WITH THE NEW CONSTANTS INDEX, PD(7) QD(8)
        pfs{i}.gen=[pfs{i}.gen mpcOPF_or2.gen(:,PMIN) mpcOPF_or2.gen(:,PMAX) mpcOPF_or2.gen(:,QMIN) mpcOPF_or2.gen(:,QMAX)]; %BE CAREFULL WITH THE NEW CONSTANTS INDEX, pmin(5) pmax(6) qmin(7) qmax(8)
        pfs{i} = fixGen2Normal_edit(gen2shunts_sol2(pfs{i}));

    end
end
disp('Creating solution2.txt')
create_solution2(pfs,contingencies,fixGen2Normal(gen2shunts(mpcOPF)));
toc(a)
end
