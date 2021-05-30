simConst = SimulationConst;
estConst = EstimatorConst;

nP = 500:500:10000;
Ks = 0.01:0.005:0.1;
nPmin = 0;
kmin = 0;
rmsemin = 100000;
for i = 1:numel(nP)
    for j = 1:numel(Ks)
        estConst.N_particles = nP(i);
        estConst.K = Ks(j);
        rmse = 0;
        for sim = 1:5
            rmse = rmse + run(simConst, estConst, false, 0);
        end
        if rmse < rmsemin
            rmsemin = rmse;
            kmin = Ks(j);
            nPmin = nP(i);
        end
    end
end