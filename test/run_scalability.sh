#!/bin/bash

# Compila il codice 
make -C ..

# Parametri di input per il tuo eseguibile
MATRIX_SIZE=256 # Nota: il tuo main attuale legge argv[1] in 'k' ma poi cicla su un vettore fisso.
FUNC="8*pi^2*sin(2*pi*x)*sin(2*pi*y)"
BOUND="0"
EXEC="../main"

# output file
OUT_FILE="scalability_report.txt"

echo "Inizio test di scalabilità..." > $OUT_FILE

# Test 1: Serial (1 MPI process, 1 OpenMP thread)
echo "Serial execution: 1 MPI rank, 1 Thread"
mpirun -np 1 $EXEC $MATRIX_SIZE 1 "$FUNC" "$BOUND" >> $OUT_FILE

# Test 2: Solo MPI (Es: 2 e 4 MPI processes, 1 Thread)
for procs in 2 4; do
    echo "Execution with $procs MPI ranks, 1 Thread"
    mpirun -np $procs $EXEC $MATRIX_SIZE 1 "$FUNC" "$BOUND" >> $OUT_FILE
done

# Test 3: Solo OpenMP ( 1 MPI process, 2 and 4 Threads)
for threads in 2 4; do
    echo "Execution with 1 MPI rank, $threads Threads"
    mpirun -np 1 $EXEC $MATRIX_SIZE $threads "$FUNC" "$BOUND" >> $OUT_FILE
done

# Test 4: Ibrido (Es: 2 MPI processes, 2 Threads)
echo "Hybrid execution: 2 MPI ranks, 2 Threads"
mpirun -np 2 $EXEC $MATRIX_SIZE 2 "$FUNC" "$BOUND" >> $OUT_FILE

echo "End of tests. Scalability report saved to $OUT_FILE"