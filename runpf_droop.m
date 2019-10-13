function [ results ] = runpf_droop( mpc, contingency, mpopt, enableRecursive, DeltaK)
%RUNPF_DROOP Run power flow taking participation factors into account
%   This function runs a power flow of the Matpower case MPC in which the
%   generators' real power injections respect their participation factors,
%   according to the Nov 30 2018 Problem Formulation.
%
%   results = RUNPF_DROOP(mpc, contingency, mpopt)
%   results = RUNPF_DROOP(mpc, contingency, mpopt, enableRecursive)
%
%   MPC must be the result of an OPF, in which one or more branches will
%   be disabled as the result of a contingency. MPC must contain the
%   original values of the OPF result.
%
%   CONTINGENCY is a struct with up to two fields: GEN and/or BRANCH. Only
%   one of them may be non-empty and must contain the unique id of the 
%   element that shall be set out of service.
%
%   ENABLERECURSIVE is a boolean that enables (TRUE) or disables (FALSE)
%   the recursive loop that enforces PV/PQ conditions at the generation
%   busses. Default: TRUE.
%
%   DELTAK input can be used for specifying an estimated value of DeltaK,
%   (e.g. the one found on a calling recursive function), which can improve
%   convergence of the case on a recursive call.
%
%   Recursive Loop
%   A recursive loop has been added to the function to try to solve all
%   PV/PQ violations at all generation busses. The recursive call blocks
%   subsequent recursive calls through the ENABLERECURSIVE input argument
%   to avoid infinite loop situations.

%   Current version: 12
%   The  swing bus is choosen according the area generator with highest reserve.
%   The runpf function was changed to runpf_edit to prevent the output of order field 
%   Contingency now enter in internal notarion, this is necessary because 
%   the field order in mpc is removed before parFor loop.
%   Contingency.gen has three fields to respect the logic of last version.

%   Current version: 11
%   Changed options for first PF run. Now Q limits are disregarded in the
%   first run.
%   Fixed bug about missing swing bus in some areas after PV/PQ conversion.

%    About version: 10
%   Changed options for first PF run. Now Q limits are disregarded in the
%   first run.
%   Fixed bug about missing swing bus in some areas after PV/PQ conversion.
%
%   About version 9B:
%   Corrected function behavior for cases with multiple areas. Also minor
%   fixes to improve time of execution.
%   9B: Fixed bug in PV/PQ enforcement
%
%   About version 8:
%   Added DELTAK input for improving convergence of recursive calls. Also
%   fixed bug in indexing for fixed generators and bug in indexing for some
%   recursive calls of generator contingencies.
%
%   About version 7:
%   Fixed two bugs with calculation of Delta: wrong delta taken if 
%   recursive droop entered and wrong APF total taken in generator 
%   contingencies.
%
%   About version 6:
%   Fixes intedified bugs related with usage of internal indexing. Adjusted
%   for cases that contain converted fixed generators.
%
%   About version 5:
%   This function now uses internal indexing. Data validation was added,
%   and some code was adapted accordingly.
%
%   About version 4:
%   contingency.branch and contingency.gen no longer contain the index of
%   the outaged element but its PSS/E id, which is converted to an index
%   internally using indexMap.
%   Generation contingencies included in formulation
%   Added recursive loop
%
%   Copyright (c) 2019, Gers USA
%   by Tomas Valencia tvalencia@gersusa.com and
%   Andres Rios andres.rios@gers.com.co
%   edited by Dario Arango dario.arango@gers.com.co

%% Validate input data
if numel(fieldnames(contingency))>2
    error('Too many fields in struct')
elseif numel(fieldnames(contingency))>1
    if isfield(contingency,'gen') && isfield(contingency,'branch')
        if ~isempty(contingency.gen) && ~isempty(contingency.branch)
            error('Only one field may be non-empty')
        elseif isempty(contingency.gen) && isempty(contingency.branch)
            error('One field must be non-empty')
        else
            isGenCON = ~isempty(contingency.gen);
        end
    else
        error('Non recognized field in struct')
    end
else
    if ~isfield(contingency,'gen') && ~isfield(contingency,'branch')
        error('Non recognized field in struct')
    else
        if isfield(contingency,'gen')
            isGenCON = ~isempty(contingency.gen);
        else
            isGenCON = 0;
        end
    end
end

if isGenCON
    if isempty(contingency.gen)
        error('Index array must be non-empty')
    end
else
    if isempty(contingency.branch)
        error('Index array must be non-empty')
    end
end

if nargin<3
    mpopt = mpoption('verbose',0, 'out.all',0,'pf.enforce_q_lims',1);
elseif isempty(mpopt)
    mpopt = mpoption('verbose',0, 'out.all',0,'pf.enforce_q_lims',1);
else
    mpopt = mpoption(mpopt);
end

%---- By default allow one recursive loop
if nargin<4
    enableRecursive = 1;
end

%%---- Check internal indexing
% if ~isfield(mpc,'order')
%     error('Case must be in internal indexing')
% end
% if mpc.order.state ~= 'i'
%     error('Case must be in internal indexing')
% end

%%--- If no DeltaK specified, set it to 0
if nargin<5
    DeltaK = 0;
    inDeltaK = false;
else
    inDeltaK = true;
end

%% Constants
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
%% Get impact areas and outaged element index in internal indexing

% Area under contingency
if ~isGenCON
    conBranch= contingency.branch;
    From = mpc.branch(conBranch,F_BUS);
    To = mpc.branch(conBranch,T_BUS);
    Area=unique([mpc.bus(From,BUS_AREA) mpc.bus(To,BUS_AREA)]);
else
    %Generator contingency
    conGenE= contingency.gen.e; %External 
    conGen=contingency.gen.ns;   %No sorted (compute before ParFor)
    if conGen==0
    %Means generator is already OFF (probably fixed gen)
     busFG =contingency.gen.bf;
     Area = mpc.bus(busFG,BUS_AREA);
    else
    conGen = contingency.gen.o; %int (sorted)
     Area = mpc.bus(mpc.bus(:,BUS_I)==mpc.gen(conGen,GEN_BUS),BUS_AREA);
    end
end
%The output in until this part are  conBranch,conGen, BusFG
%% Apply contingency
% For storing power flow result

    tempPF = mpc;

if ~isGenCON
    % Branch contingency
    tempPF.branch(conBranch,BR_STATUS)=0;
else
    % Generator contingency
    % If it is a fixed generator, application of contingency is different
    % Remove converted load
    if mpc.fixedGen(conGenE)
        tempPF.bus(busFG,PD) = tempPF.bus(busFG,PD) + contingency.gen.SMAX(1);
        tempPF.bus(busFG,QD) = tempPF.bus(busFG,QD) + contingency.gen.SMAX(2);
    else
        %Normal generator
        tempPF.gen(conGen,GEN_STATUS)=0;
    end
end

%tempPF = ext2int(int2ext(tempPF)); % Internal indexing just changed
                        % because of the outage that was just applied
%% Get participating generators
GenInArea=ismember(tempPF.gen(:,GEN_BUS),...
    tempPF.bus(ismember(tempPF.bus(:,BUS_AREA),Area),BUS_I)) & ...
    tempPF.gen(:,GEN_STATUS);
%%-- Make sure swing bus is in area of impact
tempPF = checkSwingArea(tempPF,GenInArea);

%%-Get Pg vector and APF for gens in impact area
% Dispatch in base case
Pg = tempPF.gen(GenInArea,PG);

% Normalize participation factors of generators in the impact areas
aNorm = tempPF.gen(GenInArea,APF)./sum(tempPF.gen(GenInArea,APF));
%% Redistribution
%------------------------Generation Redistribution-------------------------
% Run first PF to estimate Delta value. For this, ignore Q limits
if ~inDeltaK
    tempPF = runpf_edit(tempPF,mpoption(mpopt,'pf.enforce_q_lims',0));
    %tempPF = ext2int(tempPF); % runpf_editoutput is in external indexing
    Pgc=tempPF.gen(GenInArea,PG);
    Pgk=1.2*Pgc; % Initialize so that the while loop is entered
    DeltaK=0;
else
   % If DeltaK specified, no need for initial power flow. 
end

it = 1;
% Outer loop calculates a proxy of \Delta_k, the incremental generation
% that must be covered by all participating generators
while (inDeltaK || abs(sum(Pgc)-sum(Pgk))>1e-3) && it<5 && tempPF.success
    if ~inDeltaK
        % If DeltaK specified through input, no need for initial guess
        DeltaK=sum(Pgc)-sum(Pg);
        Pgk=1.2*Pgc;
    end
    % The inner loop increases \Delta_k to make up for the incremental
    % generation that is not covered by participating generators that reach
    % their limits.
    it_i = 1;
    while (inDeltaK || abs(sum(Pgc)-sum(Pgk))>1e-6) && it_i<300
        
        Pgk_teo=Pg+aNorm*DeltaK;
        Pgk=max(min(Pgk_teo,tempPF.gen(GenInArea,PMAX))...
            ,tempPF.gen(GenInArea,PMIN));
        if ~inDeltaK
% % %             DeltaK = DeltaK +(sum(Pgc)-sum(Pgk))./sum(aNorm(Pgk==Pgk_teo));
            if sum(aNorm(Pgk==Pgk_teo))~=0
                DeltaK = DeltaK +(sum(Pgc)-sum(Pgk))./sum(aNorm(Pgk==Pgk_teo)); %Algunos no se han saturado
            elseif max(abs(Pgk_teo-Pgk))<1e-4
                DeltaK = DeltaK +(sum(Pgc)-sum(Pgk))./sum(aNorm(abs(Pgk_teo-Pgk)<1e-4)); %Se permite tolerancia
            else
                break  %Sino, quiere decir que todos estn saturados, no se podra abastecer la demanda sin violar
                %En este caso el slack asumira la demanda faltante
            end
        else
            % If DeltaK was specified through input, do not calculate it
            % for the first iteration: take given value.
            inDeltaK = false;
            break
        end
        it_i=it_i+1;
        
    end
    tempPF.gen(GenInArea,PG)=Pgk;
    
    % Run PF and change Q enforce strategy if an error occurs
    if mpopt.pf.enforce_q_lims==1
        try
            tempPF = runpf_edit(tempPF,mpopt);
        catch ME
            switch ME.identifier
                case 'MATLAB:badsubscript'
                    % Change Q enforce strategy
                    try
                        tempPF = runpf_edit(tempPF,...
                            mpoption(mpopt,'pf.enforce_q_lims',2));
                    catch ME
                        switch ME.identifier
                            case 'MATLAB:badsubscript'
                                % Do not enforce Q
                                tempPF = runpf_edit(tempPF,...
                                    mpoption(mpopt,'pf.enforce_q_lims',0));
                            otherwise
                                rethrow(ME)
                        end
                    end
                otherwise
                    rethrow(ME)
            end
        end
    else
        tempPF = runpf_edit(tempPF, mpopt);
    end
    % runpf_editoutput is in external indexing and swing bus might have been
    % modified by Qlims enforcement
    tempPF = checkSwingArea(tempPF,GenInArea);
    Pgc=tempPF.gen(GenInArea,PG);
    it = it+1;
end

% Recursive loop to solve PV/PQ violations
% The idea here is to identify generators that have been converted to PQ
% busses and are violating PV/PQ conditions and re-convert them to PV
% busses and re-run PF.
if tempPF.success && enableRecursive
    it = 0;
    b_ind = eval_pvpqViol(tempPF);
    while it<10 && ~isempty(b_ind) && tempPF.success
        if it==0 || it>3
            mpopt = mpoption(mpopt,'pf.enforce_q_lims',0);
        else
            mpopt = mpoption(mpopt,'pf.enforce_q_lims',2);
        end
        %mpopt = mpoption(mpopt,'pf.enforce_q_lims',0);
        tempPF.bus(b_ind,BUS_TYPE) = PV;
        mpcOPFaux = mpc;
        mpcOPFaux.bus = tempPF.bus; %Le pasa las caracteristicas de los buses 
        %tempPF tiene los mismos generadores de 
        %gOn = ismember(mpcOPFaux.order.gen.status.on(mpcOPFaux.order.gen.i2e),tempPF.order.gen.status.on(tempPF.order.gen.i2e));
        gOn=(tempPF.gen(:,GEN_STATUS)==1);
        % Necessary to assign values to all online generators
        mpcOPFaux.gen(:,QG) = tempPF.gen(:,QG);
        %disp('entr\E9 a la recursividad: ')
        tempPF = runpf_droop(mpcOPFaux,contingency,mpopt,0,DeltaK);
        b_ind = eval_pvpqViol(tempPF);
    
        it=it+1;
    end
end
%results = checkRigidConstraints(tempPF);
    results = tempPF;
if isfield(tempPF,'delta')
    % If entered in recursive runpf_droop, use that delta
    results.delta = tempPF.delta;
else
    results.delta = DeltaK/sum(tempPF.gen(GenInArea,APF));
end
end

function mpc_out = checkSwingArea(mpc,GenInArea)
%CHECKSWINGAREA Checks there is only one swing bus and that it is located
%in the area of impact. If it is not, the largest eligible generator is
%taken as swing generator.
%
%   MPC_OUT = CHECKSWINGAREA(MPC, GENINAREA)
%
%   MPC must be a MatPower struct in internal ordering after applying the
%   contingency, which corresponds to TEMPF in the code above.
%
%   GENINAREA must be a column vector with same number of rows as MPC.GEN,
%   with boolean values indicating whether each generator in MPC.GEN
%   belongs to the area of impact. It is defined in the code above.


%% Constants
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

%%
aGB = mpc.gen(GenInArea,GEN_BUS);%Gen busses in the area of impact
refBus = mpc.bus(aGB,BUS_TYPE)==REF;%Ref busses in the area of impact
if sum(refBus)>1 || sum(mpc.gen(ismember(mpc.gen(:,GEN_BUS),...
            mpc.bus(aGB(refBus),BUS_I)),GEN_STATUS))==0
        % If swing bus outside area of impact or in area of impact but
        % all generators there are offline, or there is more than one
        % swing bus,define new swing bus at largest
        % generator in the area

        % Set all swing busses to PV busses
        mpc.bus(mpc.bus(:,BUS_TYPE)==REF,BUS_TYPE)=PV;
        % Eligible busses: in area and non-PQ
        eB = mpc.bus(mpc.gen(:,GEN_BUS),BUS_TYPE)==PV & GenInArea;
        if  ~any(eB)
            % Forces largest gen in area to PV node if no more PV nodes are
            % available
%             aux=[find(GenInArea),mpc.gen(GenInArea,PMAX)];
%             aux2=aux(find(mpc.gen(GenInArea,PMAX)==max(mpc.gen(GenInArea,PMAX)),1),1);
            aux=[find(GenInArea),mpc.gen(GenInArea,PMAX)-mpc.gen(GenInArea,PG)];
            aux2 = aux(aux(:,2)==max(aux(:,2)),1);
            mpc.bus(mpc.gen(aux2,GEN_BUS),BUS_TYPE)=PV;
            eB = mpc.bus(mpc.gen(:,GEN_BUS),BUS_TYPE)==PV & GenInArea;
        end
%         [~,largestGen] = max(mpc.gen(GenInArea & eB,PMAX));
        [~,largestGen] = max(mpc.gen(GenInArea & eB,PMAX)-mpc.gen(GenInArea & eB,PG));
        gA = find(GenInArea & eB);
        mpc.bus(mpc.gen(gA(largestGen(1)),GEN_BUS),BUS_TYPE)=REF;
end
% Set all swing busses outside impact area as regular PV busses
% noAreaBusses = ismember(mpc.bus(:,BUS_I), mpc.gen(~GenInArea,GEN_BUS));
% mpc.bus(noAreaBusses & mpc.bus(:,BUS_TYPE)==REF,BUS_TYPE)=PV;
mpc_out = mpc;
end
