/* Test whether omp is setup correctly
 *
 *  Compile with:
 *  g++ -fopenmp testomp.cpp -o testomp
 *  clang++ -Xclang -fopenmp -lomp testomp.cpp -o testomp
 *
 *  Run with:
 *  ./testomp
 *
 *  You should see output from multiple threads in interleaved order.
 */

#include <omp.h>
#include <stdio.h>

int main(int argc, char *argv[]) {

#pragma omp parallel
    {
        printf("Hello World... from thread = %d\n",
               omp_get_thread_num());
    }
}