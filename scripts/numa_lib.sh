#!/bin/sh

RM="/bin/rm -f"
CC=gcc
numa_test_src=$(mktemp)
numa_test=$(mktemp)

is_numa()
{
    nr_mem_nodes=$(ls -d /sys/devices/system/node/node* 2> /dev/null | wc -l)
    if [ $nr_mem_nodes -gt 100 ]; then
        return 0
    else
        return 1
    fi
}

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
    eval "$numa_test" && is_numa && echo -n $1 && exit_success

exit_failure
