# Command to Make Plots of Jet and Event Shapes
# Can run it from this directory
rivet-mkhtml -s --mc-errs -c plots.conf \
    ../results/cpp-static.yoda:"C++" \
    ../results/cuda-static.yoda:"CUDA" \
    SH-Tutorial.yoda:"Original" \