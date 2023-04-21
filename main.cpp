#define _GNU_SOURCE
#include <sched.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>

using namespace std;

cpu_set_t  mask;
inline void assignToThisCore(int core_id)
{
    CPU_ZERO(&mask);
    CPU_SET(core_id, &mask);
    sched_setaffinity(0, sizeof(mask), &mask);
}

void hnsw_impl(string stage, string using_dataset, size_t data_size, size_t attr_size, size_t label_dim,  size_t M1, size_t M2);

int main(int argc, char **argv) {

    if (argc != 8){
        printf("Usage: ./main [stage: build or search or both] [dataset] [datasize] [labelsize] [label_dim]\n");
        exit(1);
    } else {
        if (string(argv[1]) != "build" && string(argv[1]) != "search" && string(argv[1]) != "both") {
            printf("[stage: build or search or both]\n");
            exit(1);
        }
    }
    // if (string(argv[1]) == "search")
    //     assignToThisCore(0);

    hnsw_impl(string(argv[1]), string(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), atoi(argv[7]));

    return 0;
};