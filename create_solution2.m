function create_solution2( pfs, contingencies, mpcOPF)
%CREATE_SOLUTION2 This function must be called after solving a
%GOCompetition scenario with solve_opf. It creates the required
%solution2.txt file.
%   
%   RESULT = CREATE_SOLUTION2(PFS, CONTINGENCIES, WRITEFILE)
%
%   PFS must be an array of Matpower cases as returned by function
%   RUNALLCONS and converted to external numbering. Each struct in the 
%   array must correspond to a contingency in CONTINGENCIES, and they must
%   appear in the same order. This means both arrays must have the same
%   length. If any conversions (shunts, fixed generators) were done, they
%   must be undone before calling this function.
%
%   MPCOPF is the pre-contingency base case, which is used for
%   contingencies that take out elements that were already offline.

%   Current version: v6
%   BE CAREFULL: VM and VA index could change according the contingency case, if is a online contingence
%   the indexs VM and VA are definited according the cut fields in runALLcons, but if is a
%   offline or GenUniqueArea contingence the index are definited according MPC_OPF 



%   About this version:
%   Deleted argument WRITEFILE which serves no purpose anymore
%   Added argument MPCOPF and modified the code so that now cases with
%   offline contingencies are considered.
 
%
%   About version: v5
%
%   About this version:
%   Deleted argument WRITEFILE which serves no purpose anymore
%   Added argument MPCOPF and modified the code so that now cases with
%   offline contingencies are considered.
%
%   About version 4:
%   Adapted code to consider switched shunts at nonzero
%   values.
%
%   About version 3: changed input argument from single base case to
%   struct array containing all post-contingency cases. This function no
%   longer evaluates all contingencies.
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
name_solution2 = 'solution2.txt';
%% Contingencies

if ~isempty(contingencies.branch)
    conLabB = keys(contingencies.branch);
    lconKB = length(conLabB);
else
    lconKB = 0;
end
if ~isempty(contingencies.gen)
    conLabG = keys(contingencies.gen);
    lconKG = length(conLabG);
else
    lconKG = 0;
end

if isfield(contingencies,'offline')
    if ~isempty(contingencies.offline.branch)
        conOffBK = keys(contingencies.offline.branch);
        lOffB = length(conOffBK);
    else
        lOffB = 0;
    end
    if ~isempty(contingencies.offline.gen)
        conOffGK = keys(contingencies.offline.gen);
        lOffG = length(conOffGK);
    else
        lOffG = 0;
    end
    if isfield(contingencies.offline,'oneGenInArea')
        conUniqueGenInAreaK = keys(contingencies.offline.oneGenInArea);
        conUniqueGenInAreaV = values(contingencies.offline.oneGenInArea);
        lconKOneGen = length(conUniqueGenInAreaK);
    else
        lconKOneGen = 0;
    end
end

conKeys = lconKB + lconKG + lOffB + lOffG + lconKOneGen;

%% Validate data
if length(pfs)~=lconKB + lconKG + lOffB + lOffG + lconKOneGen
    error('Array length mismatch')
end

%% Calculate branch contingencies PF data and write to file
fileID = fopen(name_solution2,'w');

for k=1:conKeys
    % Get contingency label    
    
    if k<=lconKB   %branch contingency
        contingencyLabel = conLabB{k};
    elseif k<=lconKB+lconKG  %gen contingency
        contingencyLabel = conLabG{k-lconKB};
    elseif k<=lconKB+lconKG+lOffB  
        contingencyLabel = conOffBK{k-lconKB-lconKG};
    elseif k<=lconKB+lconKG+lOffB+lOffG
        contingencyLabel = conOffGK{k-lconKB-lconKG-lOffB};
    else
        contingencyLabel = conUniqueGenInAreaK{k-lconKB-lconKG-lOffB-lOffG};
    end
    
    if k<=lconKB+lconKG
        
        VM=2;
        VA=3;
        resultPF = pfs(k);
        	
    else
        
      if k> lconKB+lconKG+lOffB+lOffG && isfield(contingencies.offline,'oneGenInArea')
          VM=8;
          VA=9;
			% CON element was an uniqueGenInArea element. PF is the same OPF result
			% Set \Delta = 0
			mpcOPF.gen(mpcOPF.indexMap.gen.psse2mpc(conUniqueGenInAreaV{k-lconKB-lconKG}),[PG,QG,GEN_STATUS]) = 0;
		    resultPF = mpcOPF;
            resultPF.delta = 0;         
           
      else
        VM=8;
        VA=9;
		resultPF = mpcOPF;
		resultPF.delta = 0;
      end
    end
    
    
    % Arrange gen data into a cell, so that it can be passed to fprintf
    % function.
    genCell = cell(4,size(resultPF.gen,1));
    genCell(1,:) = num2cell(resultPF.gen(:,1)'); % Bus number
    genIDs = values(resultPF.indexMap.gen.mpc2psse,num2cell(1:size(resultPF.gen,1)));
    genIDs = cellfun(@(x)x{1},...
        regexp(genIDs,'.-(.{1,2})','tokens')); % Gen ID (without bus numbr)
    genCell(2,:) = genIDs;
    genCell(3:4,:) = num2cell(resultPF.gen(:,2:3)'); % P & Q
    
    % Get switched shunts data
    bcs = zeros(size(resultPF.bus,1),1);
    fsb = resultPF.order.bus.e2i(resultPF.comp(:,SW_BUS)); % busses w/ fixed shunts
    bcs(fsb) = resultPF.bus(fsb,BS) - resultPF.comp(:,FS); % swShunts = BS - fixed shunts
    
    % Write to file
    % Contingency label section
    fprintf(fileID,'--contingency\n');
    fprintf(fileID,'label\n');
    fprintf(fileID,'%s\n',contingencyLabel);
    % Bus section
    fprintf(fileID,'--bus section\n');
    fprintf(fileID,'bus id,v(p.u.),theta(deg),bcs(MVAR at v=1p.u)\n');
    fprintf(fileID,'%d,%2.10f,%2.10f,%2.10f\n',...
        [resultPF.bus(:,[BUS_I,VM,VA]),bcs]');
	% Generator section
    fprintf(fileID,'--generator section\n');
    fprintf(fileID,'i,unit id,pg(MW),qg(MVar)\n');
    fprintf(fileID,'%d,''%s'',%2.10f,%2.10f\n',genCell{:});

	% Delta section
    fprintf(fileID,'--delta section\n');
    fprintf(fileID,'delta (MW)\n');
    fprintf(fileID,'%2.10f\n',resultPF.delta);

end

fclose(fileID);
end