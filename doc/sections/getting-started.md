# Getting Started

The code simulates just one experiment, so this should take a little time.

The file `rungaps` can be used to build and operate both the CPU and GPU generators. It has been coded with all the routines (including the one for paper results).

NB: If you get a permission denied error, please run `chmod +x rungaps`.

This should build the program and generate 10000 events on the GPU. The output should look something like this:

```bash
--------------------------------------------------------
  ##########      ######     ###########    ##########  
 ##        ##    ##    ##    ##        ##  ##        ## 
 ##             ##      ##   ##        ##  ##           
 ##    ######  ##        ##  ###########    ##########  
 ##        ##  ############  ##                      ## 
 ##        ##  ##        ##  ##            ##        ## 
  ##########   ##        ##  ##             ##########  
                                                        
          a GPU-Amplified Parton Shower, v2.0.0         
--------------------------------------------------------
Process: LEP, E_cms: 91.2 GeV

             e+ \                       / q
                 \                     /
                  \        Z/γ        /
                   /\/\/\/\/\/\/\/\/\/
                  /                   \
                 /                     \
             e- /                       \ qbar
    
Number of Events: 10000
Running GPU only
--------------------------------------------------------
Compiling gpu-shower
Running gpu-shower...
Initialising...
 - Using 40 blocks and 256 threads per block.
 - Using LHAPDF with CT14lo set

Generating matrix elements...
Showering partons...
Completed Events: 4978/10000
partition at 5207/10000
Completed Events: 7390/10000
partition at 7550/10000
Completed Events: 8682/10000
partition at 8766/10000
Completed Events: 9358/10000
partition at 9409/10000
Completed Events: 9681/10000
partition at 9704/10000
Completed Events: 9836/10000
partition at 9850/10000
Completed Events: 9920/10000
partition at 9927/10000
Completed Events: 9958/10000
partition at 9962/10000
Completed Events: 9976/10000
partition at 9980/10000
Completed Events: 10000/10000
Veto Alg Trials: 338160
Veto Alg Errors: 0
Analysing events...

Total cross-section: 4.25e+01 nb

EVENT GENERATION COMPLETE

ME Time: 0.000324121 s
Sh Time: 0.248078 s
An Time: 0.0249141 s

Total Time: 0.273316 s

Histograms written to gpu.yoda
Timing data written to gpu-time.dat
------------------------------------------------
```

Then you have free rein over what you need; see [README.md](../../README.md) for the options. Moreover, that is all there is to it!

---

## Navigation

- [📚 Documentation Home](../README.md)
- [🏠 Repository Home](https://gitlab.com/siddharthsule/gaps)
