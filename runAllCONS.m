function [ result, pfs, voltPenalty] = runAllCONS( mpcOPF, contingencies,...
    mpcLim, pf_model, limitTime, forceV)
%RUNALLCONS Run all contingencies on a base case
%   This function applies all contingencies to a base case, one by one and
%   returns information about all the violations that occurred.
%
%   [RESULT, PFS] = RUNALLCONS(MPCOPF, CONTINGENCIES)
%   [RESULT, PFS] = RUNALLCONS(MPCOPF, CONTINGENCIES, MPCLIM)
%   [RESULT, PFS] = RUNALLCONS(MPCOPF, CONTINGENCIES, MPCLIM, PF_MODEL)
%   [RESULT, PFS] = RUNALLCONS(MPCOPF, CONTINGENCIES, MPCLIM, PF_MODEL,
%   FORCEV)
%
%   MPCOPF is the input MatPower case. It must be a pre-contingency case,
%   usually the result of an OPF. This is the case that all contingencies
%   will be run on.
%  
%   MPCLIM is an input MatPower case that will be used only to retrieve the
%   limits for evaluating constraint violations. This additional case can
%   be used if MPCOPF's limits have been modified during the optimization
%   cycle. This argument is optional. If absent, MPCOPF will be used.
%
%   CONTINGENCIES is a struct that contains all the contingencies that will
%   be evaluated. This struct must have the format returned by the
%   CONVERT2MPC function.
%
%	PF_MODEL is a string with the selected power flow model: 'DC' or 'AC'
%
%   FORCEV is a boolean that indicates whether voltages should be forced
%   within admissible values after evaluating contingencies when creating
%   output cases PFS. False by default.
%
%   RESULT is a struct where the aggregated results of the evaluation are
%   stored. The format of the struct is:
%   result
%         .sigma
%               .sigmaTh
%                       .val
%                       .sigmaTh_i
%                       .summary
%               .sigmaMis
%               .sigma_c
%                       .sigmaTh
%                       .sigmaTh_i
%                       .sigmaTh_c
%                       .H_c
%                       .H
%         .voltViol
%         .noConvs
%
%   sigmaTh.val contains the values of the slack variables of the branch
%   overload constraints. sigmaTh.val is a N-by-2 matrix, where N is the
%   number of contingencies evaluated and the two columns correspond to the
%   'from' and 'to' ends of each branch. Each row contains the sigma values
%   of the branch with highest violation for each contingency. Values are
%   in p.u. of MPCOPF.baseMVA base.
%
%   sigmaTh_i contains the index of the branch with highest violation for
%   each contingency.
%   
%   sigmaTh.summary contains the branch indices of all branches with
%   nonzero violations and its corresponding highest violation across all
%   contingencies. sigmaTh.summary is a N_b-by-3 matrix, where N_b is the
%   number of branches with overload violations in at least one
%   contingency. sigmaTh.summary(:,1) contains the branch indices;
%   sigmaTh.summary(:,2) contains the overload values in p.u. of line
%   rating and sigmaTh.summary(:,3) contains the index of the contingency
%   where the worst overload occurs.
%
%   sigmaMis contains the penalties associated with bus power unbalances
%   for each contingency. sigmaMis is a NC-by-2 matrix, where NC is the
%   number of contingencies evaluated. Columns one and two correspond to
%   real and reactive power penalties respectively. The i-th row contains 
%   the penalties caused by mismatches in power balance equations in all
%   busses for the i-th contingency. In other words:
%   sigmaMis(i,1) = sum(evalPenalty(real(POWER_MISMATCH)))
%   sigmaMis(i,2) = sum(evalPenalty(imag(POWER_MISMATCH))),
%   where POWER_MISMATCH is a NB-by-1 complex vector containing power
%   balance equation mismatches in all NB busses of the system for
%   contingency i. See EVALPENALTY for information about this function.
%
%   sigma.sigma_c is a struct that contains information about the
%   contingencies that cause branch overloads. See below for further
%   information. For all sigma.sigma_c struct data, NOL is the number of
%   contingencies that cause at least one branch overload.
%
%   sigma.sigma_c.sigmaTh and sigma.sigma_c.sigmaTh_i are 1-by-NOL cell 
%   arrays. Each element of these cell arrays is a
%   NBR-by-1 double array, where NBR is the number of overloaded elements
%   in the corresponding contingency. 
%   sigmaTh contains overloads in p.u. MVA in the base of the system 
%   (MPC.baseMVA).
%   sigmaTh_i contains the index of the overloaded branch (in internal
%   indexing).
%   Note that the index of these cell arrays does not correspond to the 
%   index of contingencies as referred in other struct data in RESULT.
%   See sigma.sigma_c.sigmaTh_c for more.
%
%   sigma.sigma_c.sigmaTh_c is the index of each contingency in the cell
%   arrays in the indexing system of all other data structures in RESULT,
%   i.e. an indexing system that goes from 1 to the total number of
%   contingencies, starting with branch contingencies and then gen
%   contingencies.
%
%   sigma.sigma_c.H_c is a 1-by-NOL cell array. Each element corresponds to
%   an overload contingency. H_c is a Linear Shift Factor matrix for the
%   overloaded lines of the corresponding contingencies. Each row
%   corresponds to one line, in the same order as in sigmaTh and sigmaTh_i.
%   Each column corresponds to a generation bus with online generation.
%
%   sigma.sigma_c.H is the vertical concatenation of all H_c matrices.
%
%   voltViol is a NC-by-1 boolean matrix depicting which contingencies
%   caused any voltage violations, where NC is the number of contingencies.
%   TRUE means a violation occurred, FALSE means no voltage violation
%   occurred.
%
%   voltPenalty is a stimation of penalization cost due to force voltage and 
%   penalize the mismtmatch power.
 
%   ----------------------------------------------------------------
%   V4 AND ABOVE NOTE: No reactive redispatch is enabled, so both columns
%   are currently identical.
%   ----------------------------------------------------------------
%
%   noConvs. Same as voltViol, for non-convergence issues.
%
%   PFS is an optional output argument. An array of N structs, where N is
%   the number of contingencies. pfs(i) is the resulting MatPower case of
%   MPCOPF after applying the i-th contingency, in internal numbering.
%
%   Contingency numbering
%   ---------------------
%   Contingencies are numbered sequentially, always
%   starting with branch contingencies, respecting the order contained in
%   struct contingencies.
%
%   About version: v10
%   Contingencies have been edited to internal notation before enter in the
%   parFor.
%   Aditional fields have been removed in the mpcOPF struct and the pfs
%   results in each contingences have been cutted to save just relevant
%   information
%   The new pfs struct is in external notation, with the fields:

%      pfs(i).bus     [Index_bus Vm VA  PD QD BS VMAX VMIN]
%      pfs(i).bus     [Index_bus Vm VA  VMAX VMIN  BS] %now
%      pfs(i).branch  [From_bus To_bus Pf Qf Pt Qt]
%      pfs(i).gen     [Gen_bus PG QG VG]

%   About version: v9
%   Modified RESULT struct, deleted old fields that are no longer used.
%   Deleted some fields from struct PFS to make it lighter.
%   Commented out the calculation of GSF matrix H to save time, since this
%   matrix is not used.
%   
%   About version 8:
%   Adapted the step that verifies QMIN and QMAX limits to new formulation
%   that uses reactive generators.
%
%   About version 7:
%   Adapted the code to the new version of SOLVESCOPF, which modifies QMAX
%   and QMIN limits of generators, which need to have their original values
%   before evaluating contingencies.
%
%   About version 6.1b:
%   Added optional arguments: DC/AC model input and forceV to force
%   voltages within admissible values.
%
%   About version 6:
%   Modifications to the output structure RESULT: added sigma_c and
%   sigmaMis
%
%   About version 5:
%   Fixed error in internal/external indexing that rendered incorrect.   
%
%   About version 4: 
%   RUNPF_DROOP now uses internal indexing, this function was adapted
%   apropriately.
%   About version 3:
%   Bugs when using pfs with parallel computing
%   corrected.
%   Reactive dispatch function commented out, as its output was being
%   ignored anyway and it takes a long time to run.
%
%   Copyright (c) 2019, Gers USA by Tomas Valencia tvalencia@gersusa.com
%   Edited by Dario Arango   dario.arango@gers.com.co

%% Validate inputs
tRACtot = tic;
if nargin<6
    forceV = false;
end
if nargin<5
    limitTime = 1e5;
end
if nargin<4
    pf_model = 'AC';
end
if ~isfield(mpcOPF,'reactgen')
    error('Missing field reactgen in mpc. Use newest version of code')
end
%% MatPower Constants
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
%% Memory pre-allocation
% Get number of contingencies
if ~isempty(contingencies.branch)
    conKeysB = values(contingencies.branch);
    lconKB = length(conKeysB);
else
    conKeysB = 0;
    lconKB = 0;
end
if ~isempty(contingencies.gen)
    conKeysG = values(contingencies.gen);
    lconKG = length(conKeysG);
else
    conKeysG = 0;
    lconKG = 0;
end

conKeys = lconKB + lconKG;

% Pre-allocate memory for all matrices
sigmaTh = zeros(conKeys,2);
sigmaTh_i = zeros(conKeys,2);
noConvs = false(conKeys,1);
voltViol = false(conKeys,1);
sigmaTh_d = cell(1,conKeys);
sigmaTh_d_i = cell(1,conKeys);
sigmaTh_c = cell(1,conKeys);
sigmaMis = zeros(conKeys,2);

% Allocate memory for pfs struct array if required
if nargout>=2
    pfs_flag = 1;
    f = fieldnames(mpcOPF)';
    f{end+1} = 'iterations';
    f{end+1} = 'delta';
    % Delete unnecessary fields to avoid data transmission overhead
    f(ismember(f,{'branch','gencost','indexMap','comp','order','reactgen','userfcn','softlims','qflow','om','x',...
        'mu','var','lin','nle','nli','qdc','raw'}))=[];
    c = cell(length(f),1);
    pfs(conKeys) = cell2struct(c,f);
else
    pfs_flag = 0;
end


%% Convert to internal numbering (necessary for RUNPF_DROOP)
mpcOPF.order.state = 'i';
mpcOPF.isxfmr = e2i_data(mpcOPF,mpcOPF.isxfmr,'branch'); % CAMBIO ALSDJFLKJDFLAJSDFLJSALDFJSDLJFLSJD
mpcOPF.order.state = 'e';
if isfield(mpcOPF,'userfcn')
    mpcOPF = rmfield(mpcOPF,'userfcn');
end
if isfield(mpcLim,'userfcn')
    mpcLim = rmfield(mpcLim,'userfcn');
end
if isfield(mpcOPF,'order')
    mpcOPF = rmfield(mpcOPF,'order');
end
if isfield(mpcLim,'order')
    mpcLim = rmfield(mpcLim,'order');
end
mpcOPF = ext2int(mpcOPF);
mpcLim = ext2int(mpcLim);
%%-- Copy generator active & reactive limits from MPCLIM to MPCOPF before
%%evaluating contingencies. This is necessary to evaluate with correct 
%%limits in case they have been modified during the optimization process.
%%This does not apply for reactive generators, as their limits are
%%voltage-dependent.
rg = ismember((1:size(mpcOPF.gen,1))',mpcOPF.reactgen);  %Indices de los compensadores reactivos, de estado activos en interno
mpcOPF.gen(~rg,[QMIN,QMAX]) = mpcLim.gen(~rg,[QMIN,QMAX]);
mpcOPF.gen(:,[PMIN,PMAX]) = mpcLim.gen(:,[PMIN,PMAX]);
%
contingency(conKeys)=struct;

%Branch contingency  
% conBranch_int=zeros(length(conKeysB),1);

%Branch Contingencies struct array in internal index - conBranch_int
 for k=1:length(conKeysB)%For because only one dimensional indexing is supported 
% conBranch_ext= mpcOPF.indexMap.branch.psse2mpc(conKeysB{k});
 %conBranch_int(k)= find(mpcOPF.order.branch.status.on==conBranch_ext); 
 contingency(k).branch =find(mpcOPF.order.branch.status.on==mpcOPF.indexMap.branch.psse2mpc(conKeysB{k}));  %Al parecer guarda la contingencua pertinente
 contingency(k).gen = [];
 end
 
%Gen contingenciess   
   for k=lconKB+1:conKeys%For because only one dimensional indexing is supported 
   contingency(k).branch =[];
   contingency(k).gen.e= mpcOPF.indexMap.gen.psse2mpc(conKeysG{k-lconKB});
   [~, contingency(k).gen.ns]= ismember(contingency(k).gen.e,mpcOPF.order.gen.status.on); %
   if contingency(k).gen.ns==0
   % Means generator is already OFF (probably fixed gen
   contingency(k).gen.bf= mpcOPF.bus(:,BUS_I)== mpcOPF.order.bus.e2i(mpcOPF.order.ext.gen(contingency(k).gen.e,GEN_BUS));
   contingency(k).gen.SMAX=[mpcOPF.gen(contingency(k).gen.e,PMAX) mpcOPF.gen(contingency(k).gen.e,QMAX)];
   else
   contingency(k).gen.o= mpcOPF.order.gen.e2i(contingency(k).gen.ns); %int (sorted)
   end
   end 

   %Save order field to posterior use
   order_save=mpcOPF.order;
   %Remove fields to avoid redundant and unused  information
   fields={'gencost','indexMap','comp','reactgen','order','softlims','om','x','mu','var','lin','nle','nli','qdc','raw'};  
   mpcOPF=rmfield(mpcOPF, fields);
   
%% Run branch contingencies
%nb, ng , nbranch
% t_rdb=tic;   %Time rundroop branch
% ticBytes(gcp);
cc=cell(size(conKeys)); %Cell to save SigmaP post contingences with violations, after force limits
mpcOPF_P = parallel.pool.Constant(mpcOPF);
% mpcOPF_P = mpcOPF;
conKeys_aux = round(conKeys / 3);
aux_loop = 1;
contEnd = 1;
tRACloop = tic;
%conKeys
%limitTime
while (toc(tRACtot)+toc(tRACloop) < limitTime) && (contEnd < conKeys)
    %fprintf('.')
    tRACloop = tic;
    contIni = aux_loop;
    contEnd = min(aux_loop+conKeys_aux, conKeys);
    aux_loop = contEnd+1;
    parfor i=contIni:contEnd
        
        pfsT = [];
        
        %SI SE RECORTAN VARIOS CAMPOS A mpcOPF que pasar�a?
        %contingency = struct; %Se crea una por cada iteraci�n? o se puede crear atras ?
        
        % Initialize temp slices
        noConvsTemp = false;
        voltViolTemp = false;
        
        % Initialize contingency struct
        %contingency.branch = conBranch_int(i);  %Al parecer guarda la contingencua pertinente
        %contingency.gen = [];
        
        % Apply contingency
             mpcPF = runpf_droop(mpcOPF_P.Value,contingency(i), mpoption('model',pf_model,'pf.enforce_q_lims',1,'verbose',0,'out.all',0,'pf.tol',1E-4)); %Entr� el MPCOPF como argumento a la funcion
%         mpcPF = runpf_droop(mpcOPF,contingency(i), mpoption('model',pf_model,'pf.enforce_q_lims',1,'verbose',0,'out.all',0,'pf.tol',1E-4)); %Entr� el MPCOPF como argumento a la funcion
        
        %mpcPF = runpf_droop(mpcOPF,contingency(i), mpoption('model',pf_model,...
        %    'pf.enforce_q_lims',1,'verbose',0,'out.all',0,'pf.tol',1E-4));
        % NOTE: For larger systems, DC PF could be required at this point in
        % order to obtain reasonable computation times.
        
        % Check convergence
        if ~mpcPF.success
            noConvsTemp = true;
        end
        
        % Check constraint violation
        [~, details, sigma_th, ~] = checkConstraints(mpcPF,1,mpcLim); %It works in intern index
        
        % Find branch with maximum overload and store for sigmaTh.val
        % Store data about all overloads for sigma.sigma_c
        [sMax,sMax_i] = max(max(sigma_th,[],2));
        sigma_aux=zeros(size(sigma_th,1),1);
        sigma_aux(:,1)=max(sigma_th,[],2);
        sigmaTh_d{i}=sigma_aux(sigma_aux>0);
        sigmaTh_d_i{i}=find(sigma_aux>0);
        if sMax>0
            sigmaTh_c{i}=i;
        else
            sigmaTh_c{i}=[];
        end
        
        sigmaTh(i,1) = sMax;
        sigmaTh_i(i,1) = sMax_i;
        
        if ~isempty(details) && mpcPF.success
            % If there were violations, check if they were voltage violations
            voltViolTemp = ~isempty(details.voltageLim);
        end
        resultPF = checkRigidConstraints(mpcPF,forceV); %Handles all notation in internal index
        
        % Function CHECKRIGIDCONSTRAINTS is always called on the resulting
        % case.
        [~,~,~,sigma_mis] = checkConstraints(resultPF,1,mpcLim); %Handles all notation in internal index
        
        % If required, handle pfs output
        if pfs_flag
            if forceV
                pfsT = forceVoltages(resultPF);
                %pfsT = force_pvpqViol(pfsT , mpcLim);
            else
                pfsT = resultPF;
            end
            %A POSSIBLE SOLUTION IS
            % Remove unnecessary fields
            %  pfsT = rmfield(pfsT,{'userfcn','softlims','om','x','mu','var','lin','nle','nli','qdc','raw'});
            if isfield(pfsT,'qflow')
                pfsT = rmfield(pfsT,'qflow');
            end
        end
        
        % Update output arrays
        noConvs(i) = noConvsTemp;
        voltViol(i) = voltViolTemp;
        sigmaMis(i,:) = [sum(evalPenalty(abs(real(sigma_mis)))),...
            sum(evalPenalty(abs(imag(sigma_mis))))];
        
        
        %Edit pfs struct
        if(voltViolTemp)  %Equivalent to   mpcFArr = arrayfun(@forceVoltages,pfs(logical(result.voltViol)));
            mpc_vv=forceVoltages(pfsT);
            [~,~,~,cc{i}] = checkConstraints(mpc_vv,1,mpc_vv);
        end
        
        %Save just index, vm, va, vmax, vmin ,BS,
        
        [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
            VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
        [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
            TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
            ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
        [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
            MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
            QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
        %%- Constants for reactive compensation matrix
        
        
        
        if forceV==1
            pfsT.order=order_save;
            pfsT=int2ext(pfsT);
            %pfsT.bus=[pfsT.bus(:,BUS_I) pfsT.bus(:,VM) pfsT.bus(:,VA) pfsT.bus(:,PD) pfsT.bus(:,QD)  pfsT.bus(:,BS) pfsT.bus(:,VMAX) pfsT.bus(:,VMIN)];
            pfsT.bus=[pfsT.bus(:,BUS_I) pfsT.bus(:,VM) pfsT.bus(:,VA) pfsT.bus(:,VMAX) pfsT.bus(:,VMIN) pfsT.bus(:,BS)];
            %pfsT.branch=[pfsT.branch(:,F_BUS) pfsT.branch(:,T_BUS) pfsT.branch(:,PF) pfsT.branch(:,QF) pfsT.branch(:,PT) pfsT.branch(:,QT)];
            pfsT.gen=[pfsT.gen(:,GEN_BUS) pfsT.gen(:,PG) pfsT.gen(:,QG) pfsT.gen(:,VG)];
            
            if ~isfield(pfsT,'iterations')
                pfsT.iterations=0;
            end
            pfsT=rmfield(pfsT,'order');
            pfsT=rmfield(pfsT,'branch');
            pfs(i)=pfsT;  %Results are saved
        else
            
            if ~isfield(pfsT,'iterations')
                pfsT.iterations=0;
            end
            if voltViol(i)==1
                pfsT.order=order_save;
                pfsT=int2ext(pfsT);
                %pfsT.bus=[pfsT.bus(:,BUS_I) pfsT.bus(:,VM) pfsT.bus(:,VA) pfsT.bus(:,PD) pfsT.bus(:,QD)  pfsT.bus(:,BS) pfsT.bus(:,VMAX) pfsT.bus(:,VMIN)];
                pfsT.bus=[pfsT.bus(:,BUS_I) pfsT.bus(:,VM) pfsT.bus(:,VA) pfsT.bus(:,VMAX) pfsT.bus(:,VMIN) pfsT.bus(:,BS)];
                %pfsT.branch=[pfsT.branch(:,F_BUS) pfsT.branch(:,T_BUS) pfsT.branch(:,PF) pfsT.branch(:,QF) pfsT.branch(:,PT) pfsT.branch(:,QT)];
                pfsT.gen=[];
                pfsT=rmfield(pfsT,'order');
                pfsT=rmfield(pfsT,'branch');
                pfs(i)=pfsT;  %Results are saved
            else
                pfsT.bus=[];
                pfsT.gen=[];
                if isfield(pfsT,'order')
                    pfsT=rmfield(pfsT,'order');
                end
                pfsT=rmfield(pfsT,'branch');
                pfs(i)=pfsT;  %Results are saved
            end            
        end
    end
end
%disp('.')
%toc(tRACtot)
%contEnd
%% Sigma
% Summary of overloads. For each line with at least one overload, find
% worst overload and store in sigmaTh.summary
uniqueLines = unique(sigmaTh_i(sigmaTh(:,1)>0 & ~noConvs(:,1),1));
sigma_uL = zeros(size(uniqueLines));
cont_i = zeros(size(uniqueLines));
for j=1:length(uniqueLines)
    uLine = uniqueLines(j);
    [sigma_uL(j), ii] = max(sigmaTh(sigmaTh_i(:,1)==uLine,1));
    sigma_uL(j) = sigma_uL(j)./... %Convert to p.u. in line base
        (mpcOPF.branch(uLine,RATE_C)/mpcOPF.baseMVA);
    iiv = find(sigmaTh_i(:,1)==uLine);
    cont_i(j)=iiv(ii);
end


%-- Convert indices back to external numbering
mpcOPF.order=order_save;
sigmaTh_i(:,1) = mpcOPF.order.branch.status.on(sigmaTh_i(:,1));
uniqueLines(:,1) = mpcOPF.order.branch.status.on(uniqueLines(:,1));

sigma.sigmaTh = struct('val',sigmaTh);
sigma.sigmaTh.sigmaTh_i = sigmaTh_i;
sigma.sigmaTh.summary = [uniqueLines, sigma_uL, cont_i];
sigma.sigmaMis = sigmaMis;


%-- Individual contingency information
% Erase empty data
for i=1:size(sigmaTh_d,2)
    if(noConvs(i,1))
        sigmaTh_d{i}=[];
        sigmaTh_d_i{i}=[];
        sigmaTh_c{i}=[];
    end
end
sigmaTh_d(cellfun('isempty',sigmaTh_d))=[];
sigmaTh_d_i(cellfun('isempty',sigmaTh_d_i))=[];
sigmaTh_c(cellfun('isempty',sigmaTh_c))=[];

%-- Output struct
sigma.sigma_c.sigmaTh = sigmaTh_d;
sigma.sigma_c.sigmaTh_i = sigmaTh_d_i;
sigma.sigma_c.sigmaTh_c = sigmaTh_c;
% sigma.sigma_c.H_c = H_c;
% sigma.sigma_c.H = H;

%% Output struct
result.voltViol = voltViol;
result.noConvs = noConvs;
result.sigma = sigma;
voltPenalty = sum(cellfun(@evalPenalty,cc));
%fprintf('Extra time: %4.3f\n',toc(ts))
end


