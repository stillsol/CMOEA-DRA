function DA = DiversityEnhancementArchiveDRA(DA, Offspring, zmin, Ns, rate)
% The Diversity Enhancement Archive of CMOEA-CD with a soft (annealed) feasibility
% preference. The baseline scores an infeasible candidate as h + CV and a feasible
% one by its angle, so an infeasible candidate can never enter while a feasible one
% exists. Here the infeasible score is alpha*h + CV with alpha = (FE/maxFE)^2:
% alpha = 1 reproduces the baseline, while alpha < 1 lets a low-CV infeasible
% solution survive as a stepping stone. Feasible candidates keep their angular term
% and the threshold h is not annealed.

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

    %% Annealing coefficient  (quadratic soft-constraint schedule)
    alpha = max(0, min(1, rate^2));

    %% Dominance relation calculation
    S  = [DA, Offspring];
    DA = [];
    NonDominated = DominationCalDRA(S, 1);
    S   = S(NonDominated);
    Obj = S.objs;
    Con = S.cons;
    CV  = sum(max(0,Con),2);
    %% Minimum angular distance calculation
    [~, M]  = size(Obj);
    [W, Ns] = UniformPoint(Ns, M);
    Angle_W_to_W = acos(1 - pdist2(W,W,'cosine'));
    Angle_W_to_W(eye(Ns) == 1) = inf;
    h = mean(min(Angle_W_to_W));

    %% Enviornmental selection
    if sum(CV <= 0) > 0
        zmax = max(Obj(CV <= 0, :), [], 1);
    else
        [~,index] = min(CV);
        zmax = Obj(index,:);
    end
    Obj          = (Obj - zmin) ./ (zmax - zmin);
    Angle_S_to_W = sin(acos(1 - pdist2(W,Obj,'cosine')));
    for i = 1 : Ns
        Angle = Angle_S_to_W(i,:);
        list  = Angle <= h;
        if sum(list) == 0
            [~,index]   = min(Angle_S_to_W(i,:));
            list(index) = true;
        end
        T        = inf(length(S), 1);
        T(list)  = CV(list);
        feasible = T <= 0;
        if sum(feasible) <= 0
            % every candidate of this subspace is infeasible: the annealed
            % penalty only rescales CV, hence still the least violating one
            [~,index] = min(T);
        else
            % a feasible candidate exists: the annealed penalty decides whether
            % an infeasible candidate can still enter (alpha < 1 -> yes)
            T           = inf(length(S), 1);
            T(list)     = alpha * h + CV(list);
            T(feasible) = Angle_S_to_W(i,feasible);
            [~,index]   = min(T);
        end
        DA = [DA,S(index)];
    end
end








