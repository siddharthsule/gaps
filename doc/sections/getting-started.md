# Getting Started

The code simulates just one experiment, so this should take a little time.

The file `rungaps` can be used to build and operate both the CPU and GPU generators. It has been coded with all the routines (including the one for paper results).

NB: If you get a permission denied error, please run ```chmod +x rungaps```.

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
                                                        
          a GPU-Amplified Parton Shower, v1.3.0         
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

Generating matrix elements...
Showering partons...
Completed Events: 10000/10000
Analysing events...

EVENT GENERATION COMPLETE

ME Time: 0.000327246 s
Sh Time: 0.0425515 s
An Time: 0.0215414 s

Total Time: 0.0644201 s

Histograms written to gpu.yoda
Timing data written to gpu-time.dat
------------------------------------------------
```

Then you have free reign over what you need, just see [README.md](../../README.md) for the options. And that's all there is to it!
