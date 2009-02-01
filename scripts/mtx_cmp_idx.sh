#!/bin/zsh
if [ -z "$2" ]; then
	echo "Usage: $0 <mmf_file> <mmf_file>"
	exit 1
fi
diff =(awk {'print $1 " " $2'} < $1) =(awk {'print $1 " " $2'} < $2)
