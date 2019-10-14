function [conIndex] = selectContingencies_nf_par(rtCL,rtCG,mpcOPF,contingencies,option)
%SELECTCONTINGENCIES select a subset of contingencies
%   [CONINDEX] = SELECTCONTINGENCIES(RTCL, RTCG,
%   MPCOPF,CONTINGENCIES,OPTION)
%
%   CONINDEX contains the indices of selected contingencies
%
%   RTCL is the selection rate in p.u. for contingencies.branch
%
%   RTCG is the selection rate in p.u. for contingencies.gen
%
%   MPCOPF is the OPF result.
%
%   CONTINGENCIES must be the CONTINGENCIES struct returned by the
%   CONVERT2MPC function.
%
%   OPTION allows to select different strategies for contingencies
%   selection
%
%   Current version: 6
%
%   About this version:
%   Parallelization is used on sorted lists computation.
%
%   About version 5:
%   'For' loops were changed by nested functions. contIndex function now
%   converts from element index to contingency index.
%
%   About version 4:
%   Sorted lists of elements are created using additional cirteria:
%   Pbr(t,f), Qbr(t,f), Pgn, Qgn
%
%   About version 3:
%   The empty option (' ') is included in main switch. This option allows
%   to select 100% of contigencies.
%
%   About version 2:
%   Check for empty fields in contingencies strcuture and retrieve
%   contingencies length
%
%   Copyright (c) 2019, Gers USA
%   by Daniel Agudelo-Martinez daniel.agudelo@gers.com.co
%% Constants
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
BASE_KV=10;
BUS_I=1;
%%

% Check input variables
if ~isempty(contingencies.branch)
    conKeysB = keys(contingencies.branch);
    lconKB = length(conKeysB);
else
    lconKB = 0 ;
end
if ~isempty(contingencies.gen)
    conKeysG = keys(contingencies.gen);
    lconKG = length(conKeysG);
else
    lconKG = 0;
end

%Different criteria for contingencies selection
switch option
    %100% of contingencies
    case{''}
        %Branch contingencies
        if lconKB>0
            conBV=values(contingencies.branch);
            [~, sLines_ix]=sortrows(abs(mpcOPF.branch(mpcOPF.branch(:,11)>0,:)),-14);
            conIndex_l= contIndex( mpcOPF.indexMap.branch, sLines_ix, conBV );
        else
            conIndex_l=0;
        end
        %Gen contingencies
        if lconKG>0
            ShuntConDPosition = mpcOPF.gen(:,10)==0 & mpcOPF.gen(:,9)==0;
            OnlyGen=mpcOPF.gen(~ShuntConDPosition,[2,3,9,4]);
            conGV=values(contingencies.gen);
            [~, sGens_ix]=sortrows(abs(OnlyGen),1);
            conIndex_g= contIndex( mpcOPF.indexMap.gen, sGens_ix, conGV )+lconKB;
        else
            conIndex_g=0;
        end
        
        %Power criteria
    case {'p'}
        %Selecting ordered contigencies
        %Branch contingencies
        if lconKB>0
            conBV=values(contingencies.branch);
            getFlow = @(x) max([sqrt(x(:,1).^2+x(:,2).^2),...
                sqrt(x(:,3).^2+x(:,4).^2)],[],2);
            powfact = @(x) cos(atan2(x(:,2),x(:,1)));
                        spmd(6)
%             sList_l=[];
            sList_BRs=[];
            sList_BRq=[];
            sList_BRqD=[];
            sList_BRv=[];
            sList_BRpf=[];
            sList_BRpfD=[];
            sList_BRconn=[];
                            switch labindex
                                case 1
            % Sorted lines according to real power flow (from)
%             getFlow = @(x) max([sqrt(x(:,1).^2+x(:,2).^2),...
%                 sqrt(x(:,3).^2+x(:,4).^2)],[],2);

%             BrS=[find(mpcOPF.branch(:,1)) getFlow(mpcOPF.branch(:,[PF,QF,PT,QT]))];
%             XfS=sortrows([BrS(mpcOPF.isxfmr,1) BrS(mpcOPF.isxfmr,2)],-2);
%             sList_xs= contIndex( mpcOPF.indexMap.branch, XfS(:,1), conBV );
%             LnS=sortrows([BrS(~mpcOPF.isxfmr,1) BrS(~mpcOPF.isxfmr,2)],-2);
%             sList_ls= contIndex( mpcOPF.indexMap.branch, LnS(:,1), conBV );
            [~, sBR_ix_S]=...
                sort(getFlow(mpcOPF.branch(:,[PF,QF,PT,QT])),'descend');
            sList_BRs= contIndex( mpcOPF.indexMap.branch, sBR_ix_S, conBV );
%             ------------------------------------------------------------------
            
            %             % Sorted lines according to real power flow (to)
            %             [~, sLines_ix_pt]=sortrows(abs(mpcOPF.branch(mpcOPF.branch(:,11)>0,:)),-16);
            %             conIndex_l_pt= contIndex( mpcOPF.indexMap.branch, sLines_ix_pt, conBV );
            % 			conIndex_l_pt=conIndex_l_pt(1:round(length(conIndex_l_pt)*rtCL));
                                case 2
            
                                    % Sorted lines according to reactive power (from)
                                    %[~, sLines_ix_qf]=sort(abs(mpcOPF.branch(:,QF)),'descend');
                                    %sList_lqf= contIndex( mpcOPF.indexMap.branch, sLines_ix_qf, conBV );
                                    %sList_l=sList_l(1:round(length(sList_lq)*rtCL));
%                                     -------------------------------------------------
%                                     BrQ=[find(mpcOPF.branch(:,1)) max(mpcOPF.branch(:,[QF,QT]),[],2)];
%                                     XfQ=sortrows([BrQ(mpcOPF.isxfmr,1) BrQ(mpcOPF.isxfmr,2)],-2);
%                                     sList_xq= contIndex( mpcOPF.indexMap.branch, XfQ(:,1), conBV );
%                                     LnQ=sortrows([BrQ(~mpcOPF.isxfmr,1) BrQ(~mpcOPF.isxfmr,2)],-2);
%                                     sList_lq= contIndex( mpcOPF.indexMap.branch, LnQ(:,1), conBV );
                                    [~, sBR_ix_Q]=...
                                                    sort(max(abs(mpcOPF.branch(:,[QF,QT])),[],2),'descend');
                                    sList_BRq= contIndex( mpcOPF.indexMap.branch, sBR_ix_Q, conBV );
                                case 3
                                    % Sorted lines according to Q diff
                                    %[~, sLines_ix_qt]=sort(abs(mpcOPF.branch(:,QT)),'descend');
                                    %sList_lqt= contIndex( mpcOPF.indexMap.branch, sLines_ix_qt, conBV );
                                    %sList_l=sList_l(1:round(length(sList_l)*rtCL));
                                    [~, sBR_ix_QD]=...
                                                    sort(abs(diff(abs(mpcOPF.branch(:,[QF,QT])),1,2)),'descend');
                                    sList_BRqD= contIndex( mpcOPF.indexMap.branch, sBR_ix_QD, conBV );
                                case 4
            %                         % Sorted lines according to power rate
            %                         [~, sLines_ix_r]=sort(abs(mpcOPF.branch(:,RATE_A)),'descend');
            %                         sList_l= contIndex( mpcOPF.indexMap.branch, sLines_ix_r, conBV );
            %                         sList_l=sList_l(1:round(length(sList_l)*rtCL));
                                    % Sorted lines according to voltage
%                                     sLines_ix_v=[];
%                                     sBus_ix_v=unique([mpcOPF.bus(mpcOPF.bus(:,10)>100,1); mpcOPF.bus(mpcOPF.bus(:,10)>100,2)]);%from&to
%                                     for ii=1:length(sBus_ix_v)
%                                     sLines_ix_v=unique([sLines_ix_v;find(mpcOPF.branch(mpcOPF.branch(:,11)==1,1)==sBus_ix_v(ii))]);
%                                     end
%                                     conIndex_l_v= contIndex( mpcOPF.indexMap.branch, sLines_ix_v, conBV );

% % %                                 for i=1:size(mpcOPF.branch,1)
% % %                                     voltFBus(i,2)=mpcOPF.bus(find(mpcOPF.branch(i,2)==mpcOPF.bus(:,1)),10);
% % %                                 end
                                %---Tests-----------------------
%                                 voltFBus=arrayfun(@(x) mpcOPF.bus(find(x==mpcOPF.bus(:,BUS_I),1),BASE_KV), mpcOPF.branch(:,F_BUS));
%                                 voltTBus=arrayfun(@(x) mpcOPF.bus(find(x==mpcOPF.bus(:,BUS_I),1),BASE_KV), mpcOPF.branch(:,T_BUS));
                                [voltFBus, voltTBus] = voltBus_function(mpcOPF.branch, mpcOPF.bus, F_BUS, T_BUS, BUS_I, BASE_KV);
%                                 -------Working only for netowrks having
%                                 less or equal branches than buses-------
%                                 [~,sBR_ix_V]=sort(max([mpcOPF.bus(mpcOPF.branch(:, F_BUS),BASE_KV),...
%                                    mpcOPF.bus(mpcOPF.branch(:,T_BUS),BASE_KV)],[],2),'descend');
%                                 ----------------------------------------
                                [~,sBR_ix_V]=sort(max([voltFBus,voltTBus],[],2),'descend');
%                                 [~,sBR_ix_V]=sort(100*(max([mpcOPF.bus(mpcOPF.branch(:,F_BUS),BASE_KV), ...
%                                     mpcOPF.bus(mpcOPF.branch(:,T_BUS),BASE_KV)],[],2)<300),'ascend');
                                 sList_BRv= contIndex( mpcOPF.indexMap.branch, sBR_ix_V, conBV );
                                case 5
%                                  Sorted lines according to power factor
                                pfT=abs(powfact(mpcOPF.branch(:,[PT,QT])));
                                pfF=abs(powfact(mpcOPF.branch(:,[PF,QF])));
                                [~,sBR_ix_pf]=sort(min([pfF,pfT],[],2),'ascend');%pf value
                                sList_BRpf= contIndex( mpcOPF.indexMap.branch, sBR_ix_pf, conBV );
                                [~,sBR_ix_pfD]=sort(abs(pfF-pfT),'descend');%pf differences
                                sList_BRpfD= contIndex( mpcOPF.indexMap.branch, sBR_ix_pfD, conBV );
                                case 6
%                                 Sorted lines according to connectivity
% % %                                 for i=1:size(mpcOPF.branch,1)
% % %                                     connBRF(i,1)=sum(mpcOPF.branch(:,1)==mpcOPF.branch(i,1));
% % %                                 end
                                %------Tests--------
%                                 connBRF=arrayfun(@(x) sum(mpcOPF.branch(:,F_BUS)==x), mpcOPF.branch(:,F_BUS));
%                                 connBRT=arrayfun(@(x) sum(mpcOPF.branch(:,T_BUS)==x), mpcOPF.branch(:,T_BUS));
                                [connBRF, connBRT] = connBR_function(mpcOPF.branch, F_BUS, T_BUS);
                                [~,sBR_ix_conn]=sort(max([connBRF connBRT],[],2),'ascend');%connectivity                                
                                sList_BRconn= contIndex( mpcOPF.indexMap.branch, sBR_ix_conn, conBV );
%                                 otherwise
            %                         %sList_l=[];
%             sList_BRs=[];
%             sList_BRq=[];
%             sList_BRqD=[];
%             sList_BRv=[];
%             sList_BRpf=[];
%             sList_BRpfD=[];
%             sList_BRconn=[];
                            end
end
            %Branch contingency indices
%             conIndex_l=unique([sList_l{1:4}]);
% % % % % % % % %             for d=1:size(sList_ls,1)
% % % % % % % % %                 aux_list(d,2)=find(d==sList_ls);
% % % % % % % % %                 aux_list(d,3)=find(d==sList_lqf);
% % % % % % % % %             end
% % %             rList_ls=1e6.*ones(size(conBV,2),1);
% % %             rList_xs=rList_ls;
% % %             rList_lqf=rList_ls;
% % %             rList_lqt=rList_ls;
% % % % % %             %rList_ls(:,1)=arrayfun(@(x) find(x==sList_ls),aux_list_ls);
% % % % % %             [~,rList_ls]=sort(sList_ls);
% % % % % %             %rList_xs(:,1)=arrayfun(@(x) find(x==sList_xs),aux_list_xs);
% % % % % %             [~,rList_xs]=sort(sList_xs);
% % % % % %             Ns = min(numel(rList_ls),numel(rList_xs));
% % % % % % 			rList_BRs=[reshape([rList_ls(1:Ns)';rList_xs(1:Ns)'],[],1);rList_ls(Ns+1:end);rList_xs(Ns+1:end)];
% % % % % %             %rList_lqf(:,1)=arrayfun(@(x) find(x==sList_lqf),aux_list);
% % % % % %             [~,rList_lq]=sort(sList_lq);
% % % % % %             %rList_lqt(:,1)=arrayfun(@(x) find(x==sList_lqt),aux_list);
% % % % % %             [~,rList_xq]=sort(sList_xq);
% % % % % %             Nq = min(numel(rList_lq),numel(rList_xq));
% % % % % % 			rList_BRq=[reshape([rList_lq(1:Nq)';rList_xq(1:Nq)'],[],1);rList_lq(Nq+1:end);rList_xq(Nq+1:end)];
% % %             %aux_list(:,4)=mean([aux_list(:,2) aux_list(:,3)],2
            [~,rList_BRs]=sort(sList_BRs{1});
            [~,rList_BRq]=sort(sList_BRq{2});
            [~,rList_BRv]=sort(sList_BRv{4});
            [~,rList_BRpf]=sort(sList_BRpf{5});
            [~,rList_BRpfD]=sort(sList_BRpfD{5});
            [~,rList_BRconn]=sort(sList_BRconn{6});
            [~,rList_BRqD]=sort(sList_BRqD{3});
% % %             [~,aux_list2]=sort(mean([...
% % %                 %min([...
% % %                 %0.1.*rList_ls...
% % %                 %,rList_xs],[],2)...
% % % 				0.1*rList_BRs
% % %                 , rList_BRq...
% % %                 ],2),'ascend');
% % % % % % 790
% % % % % % 792
% % % % % % 1744
% -----Summaries
% Connectivity:
% connF=arrayfun(@(x) sum(mpcOPF.branch(:,1)==x), unique(mpcOPF.branch(:,1)));
% connT=arrayfun(@(x) sum(mpcOPF.branch(:,2)==x), unique(mpcOPF.branch(:,2)));
% conn_summaryF=[unique(mpcOPF.branch(:,1)) connF];
% conn_summaryT=[unique(mpcOPF.branch(:,2)) connT];
% Voltage:
% conn=arrayfun(@(x) sum(mpcOPF.bus(:,10)==x), unique(mpcOPF.bus(:,10)));
% conn_summary=[unique(mpcOPF.bus(:,10)) conn]
% BusTye:
% busType=arrayfun(@(x) sum(mpcOPF.bus(:,2)==x), unique(mpcOPF.bus(:,2)));
% busType_summary=[unique(mpcOPF.bus(:,2)) busType]
% -----Criteria
% voltage:
% contIX=1012;
% max([mpcOPF.bus(mpcOPF.branch(mpcOPF.indexMap.branch.psse2mpc(conBV{contIX}),F_BUS),BASE_KV),mpcOPF.bus(mpcOPF.branch(mpcOPF.indexMap.branch.psse2mpc(conBV{contIX}),T_BUS),BASE_KV)],[],2)
% Bustype:
% contIX=1012;
% mpcOPF.bus(mpcOPF.branch(mpcOPF.indexMap.branch.psse2mpc(conBV{contIX}),1),2)
% Connectivity:
% contIX=1012;
% connF=length(find(mpcOPF.branch(mpcOPF.indexMap.branch.psse2mpc(conBV{contIX}),1)==mpcOPF.branch(:,1)))
% connT=length(find(mpcOPF.branch(mpcOPF.indexMap.branch.psse2mpc(conBV{contIX}),2)==mpcOPF.branch(:,2)))
% -----Coincidence
% rate=0.5;
% length(find(~ismember(rlvntCont.scenario_1.noConvs{1},sList_l(1:round(length(sList_l)*0.5)))))
% -----Position per list
% contIX=1574;
% a=[find(sList_BRs==contIX)/length(sList_BRs)...
% ,find(sList_BRq==contIX)/length(sList_BRq)...
% ,find(sList_BRpf==contIX)/length(sList_BRpf)...
% ,find(sList_BRpfD==contIX)/length(sList_BRpfD)...
% ,find(sList_BRconn==contIX)/length(sList_BRconn)...
% ,find(sList_BRv==contIX)/length(sList_BRv)...
% ,find(sList_BRqD==contIX)/length(sList_BRqD)...
% ]
%             aux=[rList_BRs, rList_BRq, rList_BRpf, rList_BRpfD, rList_BRconn, rList_BRv];
            w=[0.01,0.01,1,2,1, 1, 2]';
            wMean=w/sum(w);
% % % %             [~,aux_list2]=sort(mean([...
% % % %                 rList_BRs...
% % % %                 ,rList_BRq...
% % % %                 ,rList_BRpf...
% % % %                 ,rList_BRpfD...
% % % %                 ,rList_BRconn...
% % % %                 rList_BRv...
% % % %                 ],2),'ascend');
            [~,aux_list2]=sort([...
                rList_BRs...
                ,rList_BRq...
                ,rList_BRpf...
                ,rList_BRpfD...
                ,rList_BRconn...
                ,rList_BRv...
                ,rList_BRqD...
                ]*wMean,'ascend');
            %List_xfmr=find(mpcOPF.isxfmr);
            %sList_l=List_xfmr;
            %sList_l=unique([sList_l;aux_list2]);
            sList_l=aux_list2;
            
            conIndex_l=sList_l(1:round(length(sList_l)*rtCL));
        else
            conIndex_l=0;
        end
        
        %Gen contingencies
        if lconKG>0        
            conGV=values(contingencies.gen);
            
            getS = @(x) sqrt(x(:,1).^2+x(:,2).^2);
            sGen = getS(mpcOPF.gen(:,[QG,PG]));
            [~,sGEN_ix_s] = sort(sGen,'descend');
            i_r = ismember(sGEN_ix_s,mpcOPF.reactgen);
            sList_GENs= contIndex(mpcOPF.indexMap.gen, sGEN_ix_s(~i_r), conGV )+lconKB;
            
%             Sorted lines according to power factor
            powfact = @(x) cos(atan2(x(:,2),x(:,1)));
            pfGen=abs(powfact(mpcOPF.gen(:,[PG,QG])));
            [~,sGEN_ix_pf]=sort(pfGen,'ascend');
            i_r = ismember(sGEN_ix_pf,mpcOPF.reactgen);
            sList_GENpf= contIndex(mpcOPF.indexMap.gen, sGEN_ix_pf(~i_r), conGV )+lconKB;
            
%             ShuntConDPosition = mpcOPF.gen(:,10)==0 & mpcOPF.gen(:,9)==0 &...
%                 mpcOPF.gen(:,8)>0;
%             OnlyGen=mpcOPF.gen(~ShuntConDPosition,[2,3,9,4]);

%             spmd
%                 sList_g=[];
%                 switch labindex
%                     case 1
%                         % Sorted gens according to real power
%                         [~, sGens_ix_p]=sortrows(abs(OnlyGen),-1);%2
%                         sList_g= contIndex( mpcOPF.indexMap.gen, sGens_ix_p, conGV )+lconKB;
%                         sList_g=sList_g(1:round(length(sList_g)*rtCG));
%                     case 2
%                         % Sorted gens according to reactive power
%                         [~, sGens_ix_q]=sortrows(abs(OnlyGen),-2);%3
%                         sList_g= contIndex( mpcOPF.indexMap.gen, sGens_ix_q, conGV )+lconKB;
%                         sList_g=sList_g(1:round(length(sList_g)*rtCG));
%                     case 3
%                         % Sorted gens according to maximum real power
%                         [~, sGens_ix_maxp]=sortrows(abs(OnlyGen),-3);%9
%                         sList_g = contIndex( mpcOPF.indexMap.gen, sGens_ix_maxp, conGV )+lconKB;
%                         sList_g=sList_g(1:round(length(sList_g)*rtCG));
%                     case 4
%                         % Sorted gens according to maximum reactive power
%                         [~, sGens_ix_maxq]=sortrows(abs(OnlyGen),-4);%4
%                         sList_g = contIndex( mpcOPF.indexMap.gen, sGens_ix_maxq, conGV )+lconKB;
%                         sList_g=sList_g(1:round(length(sList_g)*rtCG));
%                     otherwise
%                         %sList_g=[];
%                 end
%             end
            
            %Gen contingency indices
%             conIndex_g=unique([sList_g{1:4}]);
% % %             aux_list=1:size(sList_ls,1);
% % %             rList_ls(:,1)=arrayfun(@(x) find(x==sList_ls),aux_list);
% % %             rList_lqf(:,1)=arrayfun(@(x) find(x==sList_lqf),aux_list);
            
            [~,rList_GENs]=sort(sList_GENs);
            [~,rList_GENpf]=sort(sList_GENpf);
            
             [~,aux_list3]=sort(min([...
                rList_GENs...
                ,rList_GENpf...
                ],[],2),'ascend');

            %sList_g=sList_gs;
            sList_g=aux_list3+lconKB;
            conIndex_g=sList_g(1:round(length(sList_g)*rtCG));
        else
            conIndex_g=0;
        end
        
end
%Branch and Gen contingency indices
conIndex=[conIndex_l;...
    conIndex_g];

end

function [voltFBus, voltTBus] = voltBus_function(branch, bus, F_BUS, T_BUS, BUS_I, BASE_KV)
voltFBus=arrayfun(@(x) bus(find(x==bus(:,BUS_I),1),BASE_KV), branch(:,F_BUS));
voltTBus=arrayfun(@(x) bus(find(x==bus(:,BUS_I),1),BASE_KV), branch(:,T_BUS));    

%     for i=1:size(branch,1)
%         voltFBus(i,1)=bus(find(branch(i,F_BUS)==bus(:,BUS_I),1),BASE_KV);
%         voltTBus(i,1)=bus(find(branch(i,T_BUS)==bus(:,BUS_I),1),BASE_KV);
%     end
end

function [connBRF, connBRT] = connBR_function(branch, F_BUS, T_BUS)
% for i=1:size(branch,1)
%     connBRF(i,1)=sum(branch(:,F_BUS)==branch(i,F_BUS));
%     connBRT(i,1)=sum(branch(:,T_BUS)==branch(i,T_BUS));
% end
connBRF=arrayfun(@(x) sum(branch(:,F_BUS)==x), branch(:,F_BUS));
connBRT=arrayfun(@(x) sum(branch(:,T_BUS)==x), branch(:,T_BUS));
end