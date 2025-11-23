mkdir -p rivet-plots
rivet-mkhtml --mc-errs -c test/plots.conf cpu-lep.yoda:"CPU" gpu-lep.yoda:"GPU" test/Herwig/herwig-lep.yoda:"Herwig" -o rivet-plots/lep
rivet-mkhtml --mc-errs -c test/plots.conf cpu-lepnlo.yoda:"CPU" gpu-lepnlo.yoda:"GPU" test/Herwig/herwig-lepnlo.yoda:"Herwig" -o rivet-plots/lep-nlo
rivet-mkhtml --mc-errs -c test/plots.conf cpu-lhc.yoda:"CPU" gpu-lhc.yoda:"GPU" test/Herwig/herwig-lhc.yoda:"Herwig" -o rivet-plots/lhc
rivet-mkhtml --mc-errs -c test/plots.conf cpu-lhcnlo.yoda:"CPU" gpu-lhcnlo.yoda:"GPU" test/Herwig/herwig-lhcnlo.yoda:"Herwig" -o rivet-plots/lhc-nlo