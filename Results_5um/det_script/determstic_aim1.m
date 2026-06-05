%% Run sweepCVEnergy.m Deterministic sweep --> outputs: [ATP cost, CV] for each internode length

intrnode_ls = [30:5:160];
results = sweepCVEnergy(intrnode_ls);

writetable(results.table);
writetable(results.table, "CVEnergy_table_every_5um.csv")