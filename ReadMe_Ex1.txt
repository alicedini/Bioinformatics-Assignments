To run the files it is sufficient to call them from the command line from the directory that stores them, by writing:
python [file name with .py extension] [sequence 1] [sequence 2] [match score] [mismatch score] [gap penalty]

The output of both the programs will contain:
-the best global or local alignment, from which the traceback starts
-the dynamic programming matrix
-one of the possible alignments of the two imput sequences gathered from the traceback
