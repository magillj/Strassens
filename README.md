Extra credit: Due Monday, May 9. In describing and analyzing Strassen's algorithm we assumed that we used divide and conquer all the way down to tiny matrices.
However, on small matrices the ordinary matrix multiplication algorithm will be faster because of lower overhead. This is a common issue with 
divide and conquer algorithms. The best way to run these algorithms typically is to test the input size n at the start to see if it is big 
enough to make using divide and conquer worthwhile; if n is larger than some threshold t then the algorithm would do a level of recursion, 
if n is below that threshold then it would do the non-recursive algorithm.

Your job in this question is to figure out the best choice for that threshold value for a version of Strassen's algorithm over the integers 
based on your implementation. (See the class slides for the description of the recursion used in Strassen's algorithm and for the code for 
the basic non-recursive algorithm for matrix multiplication.)

You should code up the pure algorithms first and then create the final hybrid algorithm. For simplicity you can assume that the size n 
of the matrix is a power of 2 and figure out the matrix size t=2i below which it is better to switch to the ordinary algorithm.
Your goal is to beat the ordinary algorithm by as much as possible and so find the smallest cross-over point you can. 
The language you choose to implement this in is somewhat up to you. However, the object-oriented implementation of two-dimensional arrays 
in Java with most of its standard class libraries is not great for working with two-dimensional sub-arrays. Using a language such as C 
that has more efficient array implementations and using integer arithmetic on the array indices to let you identify submatrices without 
copying them will give you better results. Check your answer for correctness against the naive algorithm. 
For your solution print out your code, the timings that you found, and the choice of t that you found works best.