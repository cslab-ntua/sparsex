#!/bin/sh

RM="/bin/rm -f"
CC=gcc
numa_test_src=$(mktemp)
numa_test=$(mktemp)

cleanup()
{
    $RM $numa_test_src $numa_test
}

exit_success()
{
    cleanup
    exit 0
}    

exit_failure()
{
    cleanup
    exit 1
}

cat > $numa_test_src <<EOF
#include <numa.h>
int main(void)
{
        if (numa_available() < 0)
                return 1;
        else
                return 0;
}
EOF
chmod u+x $numa_test

$CC -xc -lnuma $numa_test_src -o $numa_test 2> /dev/null && \
    eval "$numa_test" && echo -n $1 && exit_success

exit_failure
