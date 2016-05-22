# Strassen's Algorithm
Efficiently multiplies two input matrices by using an implementation of Strassen's algorithm in C. Switches to use the ordinary matrix multiplication when matrices hit a certain threshold of size (around 2^6). On small matrices the ordinary matrix multiplication algorithm is found to be more efficient because of lower overhead.
