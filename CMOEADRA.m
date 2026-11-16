classdef CMOEADRA < ALGORITHM
% <2026> <multi/many> <real/binary/permutation><constrained/none>
% Constrained MOEA with dynamic resource allocation for problems with large infeasible regions
% e1 --- 1 --- Type of environmental selection for forward exploration(1. SPEA2 2. NSGA-II 3. modified NSGA-III)
% e2 --- 1 --- Type of environmental selection for feasible exploitation(1. SPEA2 2. NSGA-II 3. modified NSGA-III)
% adapt --- 1 --- Enable the soft-constraint annealing in the diversity archive (0/1)
% faHi --- 50 --- Parent quota of the forward exploration archive at the start of the run
% faFloor --- 30 --- Parent quota of the forward exploration archive at the end of the run
% earlyFrac --- 0.5 --- Fraction of the budget during which the one-way merge is active
% feasGate --- 25 --- Merge runs only while FEA holds fewer than this many feasible members

%------------------------------- Reference --------------------------------
% F. Luo, Y. Hu, Y. Tang, Z. Qian, S. Liu, H. Yuan, and W. Deng. A
% constrained multiobjective evolutionary algorithm with dynamic resource
% allocation for problems with large infeasible regions. Submitted, 2026.
% The three-archive framework is adapted from CMOEA-CD (Z. Liu, F. Han,
% Q. Ling, H. Han, J. Jiang, IEEE Trans. Evol. Comput. 29(6) (2025) 2771-2784).
%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
% Modifications for CMOEA-DRA (c) 2026 F. Luo, Y. Hu, Y. Tang, Z. Qian,
% S. Liu, H. Yuan, and W. Deng.

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            % earlyFrac : fraction of the budget during which the one-way merge is active
            % feasGate  : the merge runs only while FEA holds fewer than this many
            %             feasible members
            [e1,e2,adapt,faHi,faFloor,earlyFrac,feasGate] = ...
                Algorithm.ParameterSet(1,1,1,50,30,0.5,25);
            N  = Problem.N;
            Ns = floor(N/3);

            %% Generate random population
            Population = Problem.Initialization();
            FA  = [];
            DA  = [];
            FEA = Population;
            Offspring = Population;
            zmin      = min(Population.objs,[],1) - 1e-6;

            %% Optimization
            while Algorithm.NotTerminated(FEA)
                zmin = min(zmin, min(Offspring.objs, [], 1) - 1e-6);
                % Parent counts (archives keep their capacities; only the number of
                % parents drawn from each archive changes):
                %   nFA = min(faHi, max(faFloor, round(faFloor + (faHi-faFloor)*w)))
                %   nDA = Ns,   nFE = 3*Ns - nFA - nDA
                if adapt
                    rate = min(1, Algorithm.pro.FE/Algorithm.pro.maxFE);
                    w    = 0.5*(1 + cos(pi*rate));
                    % w DECREASES 1 -> 0, so the count must be floor + (hi-floor)*w.
                    % Written as (hi - (hi-floor)*w) it would run BACKWARDS: the
                    % declared maximum would be reached at the END of the run.
                    nFA  = min(faHi, max(faFloor, round(faFloor + (faHi-faFloor)*w)));
                    nDA  = Ns;
                    nFE  = 3*Ns - nFA - nDA;
                    daRate = rate;
                else
                    nFA = Ns;  nDA = Ns;  nFE = Ns;
                    daRate = 1;
                end
                % Archive updates. DA carries the quadratic soft-constraint annealing:
                % an infeasible candidate is scored alpha*h + CV with alpha = rate^2, so
                % early in the run a low-violation infeasible solution can outrank a
                % feasible one; daRate = 1 (adapt = 0) reproduces the baseline.
                FA   = ForwardExplorationArchiveDRA(FA, Offspring, zmin, Ns, e1);
                % DA is updated first and is never modified afterwards, which is what
                % makes the merge one-way: DA never sees FEA.
                DA   = DiversityEnhancementArchiveDRA(DA, Offspring, zmin, Ns, daRate);
                % FEA draws on DA's members through its candidate pool only: its own
                % selection, capacity and output are unchanged. The merge is applied in
                % the early part of the run (t < earlyFrac) and only while FEA holds
                % fewer than feasGate feasible members.
                if isempty(FEA.cons)
                    nFeasNow = N;      % unconstrained: feasibility is never scarce
                else
                    nFeasNow = sum(sum(max(0,FEA.cons),2) <= 0);
                end
                if adapt && min(1,Algorithm.pro.FE/Algorithm.pro.maxFE) < earlyFrac && nFeasNow < feasGate
                    PoolFEA = [FEA, DA, Offspring];
                else
                    PoolFEA = [FEA, Offspring];
                end
                FEA  = FeasibilityExploitationArchiveDRA(PoolFEA, [], N, e2);
                Pop1 = FA;
                Pop2 = DA;
                % The draw is with replacement when nFA exceeds |FA| = Ns and without
                % replacement otherwise; the M >= 3 clamp on DA is applied before nFE
                % is derived, so the slots it frees go to FEA.
                if adapt
                    if numel(Pop1) >= nFA
                        if numel(Pop1) > nFA, Pop1 = Pop1(randperm(numel(Pop1),nFA)); end
                    else
                        Pop1 = Pop1(unidrnd(numel(Pop1), [1, nFA]));
                    end
                    if numel(Pop2) > nDA, Pop2 = Pop2(randperm(numel(Pop2),nDA)); end
                    nFA = numel(Pop1);
                    nDA = numel(Pop2);
                    nFE = 3*Ns - nFA - nDA;
                    assert(nFE >= 2,'CMOEADRA:budget','nFE fell below 2');
                end
                Pop3 = FEA(unidrnd(length(FEA), [1, nFE]));
                MatingPool_Pop1_1 = randperm(length(Pop1));
                MatingPool_Pop1_2 = randperm(length(Pop1));
                MatingPool_Pop2_1 = randperm(length(Pop2));
                MatingPool_Pop2_2 = randperm(length(Pop2));
                MatingPool_Pop3_1 = randperm(length(Pop3));
                MatingPool_Pop3_2 = randperm(length(Pop3));
                if rand() < 0.5
                    Offspring1 = OperatorDE(Problem, Pop1, Pop1(MatingPool_Pop1_1), Pop1(MatingPool_Pop1_2),{1,0.5,1,1});
                    Offspring2 = OperatorDE(Problem, Pop2, Pop2(MatingPool_Pop2_1), Pop2(MatingPool_Pop2_2),{1,0.5,1,1});
                    Offspring3 = OperatorDE(Problem, Pop3, Pop3(MatingPool_Pop3_1), Pop3(MatingPool_Pop3_2),{1,0.5,1,1});
                else
                    Offspring1 = OperatorGA(Problem, Pop1(MatingPool_Pop1_1), {1,20,1,1});
                    Offspring2 = OperatorGA(Problem, Pop2(MatingPool_Pop2_1), {1,20,1,1});
                    Offspring3 = OperatorGA(Problem, Pop3(MatingPool_Pop3_1), {1,20,1,1});
                end
                Offspring = [Offspring1, Offspring2, Offspring3];
            end
        end
    end
end
