# Command to Make Plots of Jet and Event Shapes
# Write something similar in the terminal
rivet-mkhtml -s --mc-errs -c plots.conf \
    ../cpp.yoda:"C++" \
    ../gaps.yoda:"GAPS" \
    SH-Tutorial.yoda:"S. H." \