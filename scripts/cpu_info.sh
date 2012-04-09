#!/bin/bash

vendor=$(cat /proc/cpuinfo | grep 'vendor_id' | head -n1 | awk -F': ' '{print $2}')
if [ "$vendor" = "AuthenticAMD" ]; then
    # AMD CPU
    echo "AMD"
elif [ "$vendor" = "GenuineIntel" ]; then
    # Intel CPU
    # Query the cpu family /proc/cpuinfo
    #   6: Intel Core
    #  15: Intel Netburst
    family=$(cat /proc/cpuinfo | grep 'cpu family' | head -n1 | awk -F': ' '{print $2}')
    if [ "$family" = "6" ]; then
        echo "CORE"
    elif [ "$family" = "15" ]; then
        echo "XEON"
    fi
else
    # Other CPU; capitalize the vendor id
    echo $(echo $vendor | tr a-z A-Z | tr - _)
fi
