function create_solution1( mpc )
%CREATE_SOLUTION1 This routine creates the 'solution1.txt' file.
%
%   The function reads the dispatch data from MPC, which must be a MatPower
%   case as returned by the convert2mpc function, and creates the 
%   solution1.txt file according to the rules of the GO Competition.
%
%   Current version: v5
%
%   About this version:
%   Deleted argument WRITEFILE which serves no purpose anymore
%
%   About version 4:
%   Adapted code to consider switched shunts at nonzero
%   values.
%
%   About version 3:
%   Switched shunts not considered. Needs to be updated.
%   Changed format to match Challenge1 
%
%   Copyright (c) 2019, Gers USA
%   by Tomas Valencia tvalencia@gersusa.com

%% Define named indices into gen, comp matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
%%- Constants for reactive compensation matrix
SW_BUS=1;% Bus nb
FS=2;% Fixed shunt (MVAR @ 1.0 p.u)
SSMIN=4;% Switched shunt min (MVAR @ 1.0 p.u)
SSMAX=3;% Switched shunt max (MVAR @ 1.0 p.u)
SW_STATUS=5; % Switched shunt status



%%
name_solution1 = 'solution1.txt';
%% Write data to file

% Arrange gen data into a cell, so that it can be passed to fprintf
% function.
genCell = cell(4,size(mpc.gen,1));
genCell(1,:) = num2cell(mpc.gen(:,1)'); % Bus number
genIDs = values(mpc.indexMap.gen.mpc2psse,num2cell(1:size(mpc.gen,1)));
genIDs = cellfun(@(x)x{1},...
    regexp(genIDs,'.-(.{1,2})','tokens')); % Gen ID (without bus numbr)
genCell(2,:) = genIDs;
genCell(3:4,:) = num2cell(mpc.gen(:,PG:QG)'); % P & Q

bcs = zeros(size(mpc.bus,1),1);
fsb = mpc.order.bus.e2i(mpc.comp(:,SW_BUS)); % busses w/ fixed shunts
bcs(fsb) = mpc.bus(fsb,BS) - mpc.comp(:,FS); % swShunts = BS - fixed shunts
fileID = fopen(name_solution1,'w');
fprintf(fileID,'--bus section\n');
fprintf(fileID,'bus id,v(p.u.),theta(deg),bcs(MVAR at v=1p.u)\n');
fprintf(fileID,'%d,%2.10f,%2.10f,%2.10f\n',[mpc.bus(:,[BUS_I,VM,VA]),...
     bcs]');
fprintf(fileID,'--generator section\n');
fprintf(fileID,'i,unit id,pg(MW),qg(MVar)\n');
fprintf(fileID,'%d,''%s'',%2.10f,%2.10f\n',genCell{:});
fclose(fileID);


end
