function MyMatlab1(InFile1, InFile2, InFile3, InFile4,...
    TimeLimitInSeconds, ScoringMethod, NetworkModel)
%MYMATLAB1 Version1 Script. This script generates solution1.txt
%
% Input Arguments:
% TimeLimitInSeconds and ScoringMethod are not used yet, the algorithm acts
% as it would for ScoringMethod = 0.
% InFile1: CON file
% InFile2: INL file
% InFile3: RAW file
% InFile4: ROP file


fprintf('Started\n')
% if isempty(gcp('nocreate'))
% parpool('local',4);
% end
startTime = tic;
[mpc,contingencies] = convert2mpc(InFile3,InFile4,InFile2,InFile1);
% [mpc,contingencies] = convert2mpc_par(InFile3,InFile4,InFile2,InFile1);
%TimeLimitInSeconds=600;
%ScoringMethod=2;

[mpcOPF, ~, mpcOPF_or] = solveSCOPF(mpc,contingencies,...
    true,TimeLimitInSeconds, ScoringMethod, startTime);
save('mpc.mat','mpcOPF','mpcOPF_or','contingencies');
disp('Creating solution1.txt')
a=tic;
create_solution1(fixGen2Normal(gen2shunts(mpcOPF)));
toc(a)
end
