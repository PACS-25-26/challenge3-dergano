#!/bin/bash

make -C ..

# Input parameters
MATRIX_SIZE=256 
FUNC="8*pi^2*sin(2*pi*x)*sin(2*pi*y)"
BOUND="0"
EXEC="../laplace"

# output file
OUT_FILE="scalability_report.txt"

echo "Start scalability test" > $OUT_FILE
echo "" >> $OUT_FILE
# Test 1: Serial (1 MPI process, 1 OpenMP thread)
echo "Serial execution: 1 MPI rank, 1 Thread">>$OUT_FILE
mpirun -np 1 $EXEC $MATRIX_SIZE 1 "$FUNC" "$BOUND" >> $OUT_FILE
echo "" >> $OUT_FILE
# Test 2: Solo MPI (Es: 2 e 4 MPI processes, 1 Thread)
echo "Testing MPI scalability with 1 Thread">>$OUT_FILE
for procs in 2 4; do
    echo "Execution with $procs MPI ranks, 1 Thread">>$OUT_FILE
    mpirun -np $procs $EXEC $MATRIX_SIZE 1 "$FUNC" "$BOUND" >> $OUT_FILE
done
echo "" >> $OUT_FILE
# Test 3: Solo OpenMP ( 1 MPI process, 2 and 4 Threads)
echo "Testing OpenMP scalability with 1 MPI rank">>$OUT_FILE    
for threads in 2 4; do
    echo "Execution with 1 MPI rank, $threads Threads">>$OUT_FILE
    mpirun -np 1 $EXEC $MATRIX_SIZE $threads "$FUNC" "$BOUND" >> $OUT_FILE
done
echo "" >> $OUT_FILE
# Test 4: Ibrido (Es: 2 MPI processes, 2 Threads)
echo "Hybrid execution: 2 MPI ranks, 2 Threads">>$OUT_FILE
mpirun -np 2 $EXEC $MATRIX_SIZE 2 "$FUNC" "$BOUND" >> $OUT_FILE
echo "" >> $OUT_FILE
echo "End of tests. Scalability report saved to $OUT_FILE"