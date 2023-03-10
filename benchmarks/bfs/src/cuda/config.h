#define MAX_THREADS_PER_BLOCK 512

// The number of Streaming Multiprocessors; 15 for Fermi architecture 30 for G280 at the moment of this document
#define NUM_SM 14

// The number of duplicated frontiers used in BFS_kernel_multi_blk_inGPU
#define NUM_BIN 8

// EXP = log(NUM_BIN), assuming NUM_BIN is still power of 2 in the future architecture
// Using EXP and shifting can speed up division operation
#define EXP 3

// This variable is also related with NUM_BIN; may change in the future architecture;
// Using MOD_OP and "bitwise and" can speed up mod operation
#define MOD_OP 7

#define INFINITY 2147483647 // 2^31-1
#define UP_LIMIT 16677216   // 2^24
#define WHITE 16677217
#define GRAY 16677218
#define GRAY0 16677219
#define GRAY1 16677220
#define BLACK 16677221
#define W_QUEUE_SIZE 400
