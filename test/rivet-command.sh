# Command to Make Plots of Jet and Event Shapes
# Write something similar in the terminal
rivet-mkhtml -s --mc-errs -c plots.conf \
    ../results-times/cpp.yoda:"C++" \  
    ../results-times/gaps.yoda:"CUDA" \
    SH-Tutorial.yoda:"S. H."