# Command to Make Plots of Jet and Event Shapes
# Write something similar in the terminal
rivet-mkhtml -s --mc-errs -c plots.conf \
    ../results-times/cpp.yoda:"CPU" \  
    ../results-times/gaps.yoda:"GPU" \
    SH-Tutorial.yoda:"S. H."