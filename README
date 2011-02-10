1. Introduction
===============

This package is a proof-of-concept release of the Compressed Sparse eXtended
format for sparse matrices. This format seeks to minimize the memory footprint
of the column index array of the typical Compressed Sparse Row (CSR) format by
exploiting dense substructures inside the sparse matrix. Instead of storing a
single index for every nonzero element of the sparse matrix, CSX stores a short
description for each substructure found in the matrix (and selected for
encoding). This technique can save significant amount of main memory storage and
minimize the bandwidth requirements of the Sparse Matrix-Vector Multiplication
(SpMV) kernel. Finally, the CSX format employes runtime code generation (using
the LLVM compiler infrastructure) for emitting optimized SpMV routines for each
encoded pattern.

More information about the CSX format can be found in

K. Kourtis, V. Karakasis, G. Goumas, and N. Koziris,
"CSX: An extended compression format for SpMV on shared memory systems,"
16th ACM SIGPLAN Annual Symposium on Principles and Practice of Parallel
Programming (PPoPP'11) San Antonio, TX, USA, February 12-16, 2011.


2. Prerequisites and external dependencies
==========================================

This package depends on LLVM 2.7 and the LLVM GCC front-end compiler (llvm-gcc).

The package will not compile with LLVM 2.5, because the LLVM API has changed
since then. However, the package has not been tested against LLVM 2.6, but it
might compile without errors for this release, too.


2. Compiling the package
========================

This package is written in C and C++, so `g++' is required for compilation.

If you have properly installed and configured LLVM 2.7 and `g++', simply run
`make' in the top-level source directory of the package.


3. Running CSX
==============

After you have successfully compiled CSX, the final executable, named `spmv',
will be placed inside the `csx' directory. The spmv executable may be invoked as
follows:

[ENV=<value>] ... ./spmv <matrix_file> [...]

The execution of CSX is controlled by a set of environment variables and takes
as input one or more files describing sparse matrices. For each matrix, `spmv'
will perform an analysis of the patterns found and print the results, the
preprocessing time and the wall-clock time (and performance in Mflop/s) for the
execution of 128 consecutive SpMV operations.

3.1 Input matrix format
-----------------------

The `spmv' executable requires the input matrices to be in a variation of the
Matrix Market Exchange format. Specifically, the first line of the file must
contain the number of rows, columns and nonzeros, respectively, separated by
whitespace. The rest of the lines contain the nonzero elements of the matrix
(one per line) in the form of `<row> <column> <value>', sorted lexicographically
(i.e., row-major). To facilitate the conversion to the desired format, the
utility script `scripts/sort-mtx.sh' is supplied.

See also http://math.nist.gov/MatrixMarket/formats.html

3.2 Environment
---------------

The execution of `spmv' is solely controlled by environment variables, which are
described in detail below.

** Variables controlling multithreaded execution

MT_CONF         Set the cpu affinity for running `spmv'. This is the only way to
                setup a multithreaded execution for `spmv'. You should supply
                the cpu numbers (see /sys/devices/) of the desired cpus as a
                comma-separated list. If MT_CONF is not specified, it is set to
                `0', thus assuming single-threaded execution.

** Variables controlling the construction of CSX

XFORM_CONF      Set the pattern types to search for in the matrix. This is a
                comma-separated list of the pattern type IDs as specified in the
                enumeration `SpmIterOrder' in `spm.h'. The mapping between type
                IDs and types is the following:

                           0 -> NONE
                           1 -> HORIZONTAL
                           2 -> VERTICAL
                           3 -> DIAGONAL
                           4 -> REV_DIAGONAL
                           5 -> BLOCK_TYPE_START (not used)
                           6 -> BLOCK_R1 (not used)
                           7 -> BLOCK_R2
                           8 -> BLOCK_R3
                           9 -> BLOCK_R4
                          10 -> BLOCK_R5
                          11 -> BLOCK_R6
                          12 -> BLOCK_R7
                          13 -> BLOCK_R8
                          14 -> BLOCK_COL_START (not used)
                          15 -> BLOCK_C1 (not used)
                          16 -> BLOCK_C2
                          17 -> BLOCK_C3
                          18 -> BLOCK_C4
                          19 -> BLOCK_C5
                          20 -> BLOCK_C6
                          21 -> BLOCK_C7
                          22 -> BLOCK_C8
                          23 -> BLOCK_TYPE_END (not used)
                          24 -> XFORM_MAX
                
                The default value of XFORM_CONF is `0', i.e., NONE.

ENCODE_DELTAS   Set specific deltas to encode for one-dimensional patterns or
                block dimensions for two-dimensional patterns. ENCODE_DELTAS is
                a comma-separated list of "delta descriptors". A delta
                descriptor is a tuple in the form of `{d1,d2,...,dn}', where
                `di' is a specific delta to be encoded. For two-dimensional
                patterns types, the `di' signifies the second dimension of the
                block pattern. The ENCODE_DELTAS variable must be specified
                along with XFORM_CONF, in which case there is a one-to-one
                mapping between the pattern types and delta descriptors. For
                example, using XFORM_CONF=16,1,3 and
                ENCODE_DELTAS={4,7},{1},{11}, first `BLOCK_C2' patterns with a
                second dimension of 4 and 7, respectively (i.e., block patterns
                4x2 and 7x2) will be encoded into the matrix. Second,
                `HORIZONTAL' patterns with delta 1 will be encoded and, finally,
                `DIAGONAL' patterns with delta 11. The number of pattern types
                specified in XFORM_CONF must equal the number of delta
                descriptors specified in ENCODE_DELTAS.

SPLIT_BLOCKS    Enable splitting of larger blocks into smaller ones, when
                searching for block patterns. This option usually leads to
                better compression, since it improves the statistics (number of
                nonzero elements covered by a pattern) of smaller blocks making
                them eligible for encoding. However, this option should be
                avoided, if the improved preprocessing using sampling windows is
                to be used, since it can lead to inflated statistics, which may
                not correspond to the actual distribution of the patterns in the
                matrix. The value of this variable is ignored and it might be
                null as well.

** Variables controlling preprocessing

WINDOW_SIZE     If set and nonzero, enable preprocessing using sampling
                windows. This variable specifies the size of the sampling
                windows either in terms of number of rows or number of nonzero
                elements. The default policy is to extract sampling windows
                based on the number of nonzero elements, in which case
                WINDOW_SIZE denotes how many nonzero elements each sampling
                window will have. The size of the window will be "rounded" up to
                the next full row. The matrix could also be sampled row-wise, in
                which case WINDOW_SIZE is the number or rows of each
                window. However, to enable this policy, you should alter
                `spmv.cc' and recompile. Sampling based on the nonzero elements
                (default policy) is recommended, since it can provide better
                coverage of the matrix, especially if used in conjunction with
                the SAMPLES variable (see below), and lower preprocessing cost.

SAMPLING_PROB   If window preprocessing is enabled, specify the probability of
                selecting a window for pattern searching. The original matrix is
                statically partitioned into preprocessing windows and
                SAMPLING_PROB specifies the probability of considering each
                window for pattern search. The statistics gathered from all the
                selected windows are projected to the whole matrix and the
                selection of patterns types to encode proceeds as in the normal
                case (without sampling windows).

SAMPLES         Specify the maximum number of sampling windows to consider for
                searching. This variable can be used either in conjunction with
                the SAMPLING_PROB or alone. In the latter case, the sampling
                probability is computed automatically, so that the selected
                windows span the whole matrix. This variable in conjunction with
                a reasonable window size can provide a very good coverage of the
                nonzero elements of the original matrix. For our experiments, we
                usually use 5 sampling windows with 32K nonzero elements each.

4. Licence
==========

This package is distributed under the BSD Licence. See `LICENCE.txt' for more
information.


5. Authors
==========

This package is maintained by 

Kornilios Kourtis          <kkourt@cslab.ece.ntua.gr>
Vasileios Karakasis        <bkk@cslab.ece.ntua.gr>
Theodoros Goudouvas        <thgoud@cslab.ece.ntua.gr>
