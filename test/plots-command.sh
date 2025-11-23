mkdir -p rivet-plots
rivet-mkhtml --mc-errs -c test/plots.conf cpu-lep.yoda gpu-lep.yoda test/Herwig/herwig-lep.yoda -o rivet-plots/lep
rivet-mkhtml --mc-errs -c test/plots.conf cpu-lepnlo.yoda gpu-lepnlo.yoda test/Herwig/herwig-lepnlo.yoda -o rivet-plots/lep-nlo
rivet-mkhtml --mc-errs -c test/plots.conf cpu-lhc.yoda gpu-lhc.yoda test/Herwig/herwig-lhc.yoda -o rivet-plots/lhc
rivet-mkhtml --mc-errs -c test/plots.conf cpu-lhcnlo.yoda gpu-lhcnlo.yoda test/Herwig/herwig-lhcnlo.yoda -o rivet-plots/lhc-nlo