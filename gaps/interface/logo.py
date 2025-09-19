#!/usr/bin/env python3

import os
import subprocess

# ------------------------------------------------------------------------------
# GAPS LOGO


def logo():

    # ANSI escape codes for colors
    blue_colour = "\033[34m"
    green_colour = "\033[32m"
    reset_col = "\033[0m"

    def print_blue(text):
        print(f"{blue_colour}{text}{reset_col}")

    def print_green(text):
        print(f"{green_colour}{text}{reset_col}")

    # Print the ASCII art
    print("--------------------------------------------------------")
    print_blue("  ##########      ######     ###########    ##########  ")
    print_blue(" ##        ##    ##    ##    ##        ##  ##        ## ")
    print_blue(" ##             ##      ##   ##        ##  ##           ")
    print_blue(" ##    ######  ##        ##  ###########    ##########  ")
    print_blue(" ##        ##  ############  ##                      ## ")
    print_blue(" ##        ##  ##        ##  ##            ##        ## ")
    print_blue("  ##########   ##        ##  ##             ##########  ")
    print_blue("                                                        ")
    print_green("          a GPU-Amplified Parton Shower, v2.0.0         ")
    print("--------------------------------------------------------")

# ------------------------------------------------------------------------------
# LEP LO


def lep_lo():
    print(r"""
             e+ \                       / q
                 \                     /
                  \        Z/γ        /
                   /\/\/\/\/\/\/\/\/\/
                  /                   \
                 /                     \
             e- /                       \ qbar
    """)


# ------------------------------------------------------------------------------
# LHC LO


def lhc_lo():
    print(r"""
              q \                       / e+
                 \                     /
                  \        Z/γ        /
                   /\/\/\/\/\/\/\/\/\/
                  /                   \
                 /                     \
           qbar /                       \ e-
    """)

# ------------------------------------------------------------------------------
# LEP NLO


def lep_nlo():
    print(r"""
 e+ \               / q        e+ \               / q
     \             /               \             /
      \    Z/γ    /                 \    Z/γ    /-oo
       /\/\/\/\/\/                   /\/\/\/\/\/   o
      /           \                 /           \-oo
     /             \               /             \
 e- /               \ qbar     e- /               \ qbar
          

 e+ \           ooo-/ q        e+ \               / q
     \          o  /               \             /
      \    Z/γ  o-/                 \    Z/γ    / 
       /\/\/\/\/\/                   /\/\/\/\/\/   
      /           \                 /         o-\
     /             \               /          o  \
 e- /               \ qbar     e- /           ooo-\ qbar
          
           
 e+ \               / q        e+ \               / q
     \             /               \             /
      \    Z/γ    /-ooo- g          \    Z/γ    /
       /\/\/\/\/\/                   /\/\/\/\/\/    
      /           \                 /           \-ooo- g
     /             \               /             \
 e- /               \ qbar     e- /               \ qbar  
    """)

# ------------------------------------------------------------------------------
# LHC NLO


def lhc_nlo():
    print(r"""
    q \                         q \                   
       \                           \              
        \     Z                  oo-\     Z       
         /\/\/\/\/\/             o   /\/\/\/\/\/
        /                        oo-/            
       /                           /              
 qbar /                      qbar /                   
          

    q \-ooo                      q \                 
       \  o                         \             
        \-o   Z                      \    Z        
         /\/\/\/\/\/                 /\/\/\/\/\/   
        /                           /-o          
       /                           /  o           
 qbar /                      qbar /-ooo               
          
           
    q \                          q \              
       \-ooo- g                     \             
        \     Z                      \    Z        
         /\/\/\/\/\/                 /\/\/\/\/\/   
        /                           /           
       /                           /-ooo- g        
 qbar /                      qbar /                   
                                     
                                     
         / qbar                  q \                  
 g -ooo-/                           \             
        \     Z                      \    Z        
         /\/\/\/\/\/                 /\/\/\/\/\/   
        /                           /            
       /                     g -ooo-\             
 qbar /                              \ q              
    """)

# ------------------------------------------------------------------------------
# Print the logo and the process


def print_logo(args):

    # Print the logo
    logo()

    # Print the process and parameters
    if args.process == 'LEP':
        print(f"Process: LEP, E_cms: {args.root_s} GeV")
        if args.nlo:
            lep_nlo()
        else:
            lep_lo()

    elif args.process == 'LHC':
        print(f"Process: LHC: Root S = {args.root_s / 1000} TeV")
        if args.nlo:
            lhc_nlo()
        else:
            lhc_lo()

    # Print the number of events
    print(f"Number of Events: {args.nevents}")

    # Print the run type
    if args.runtype == 'gpu':
        print("Running GPU only")
    elif args.runtype == 'cpu':
        print("Running CPU Shower only")
    elif args.runtype == 'compare':
        print("Running GPU and CPU Shower and comparing the results")
    elif args.runtype == 'full':
        print("Running GPU and CPU Shower for a range of event numbers")

    print("--------------------------------------------------------")
