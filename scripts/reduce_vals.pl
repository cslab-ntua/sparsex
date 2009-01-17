#!/usr/bin/perl -w
# replace the values of a matrix file with a limited number of vals 

@vals = ("1.00", "-1.00");
while (<>){
	print;
	if (/^#/){
		next;
	} else {
		last;
	}
}

while (<>){
	if ( !/(\S+\s+\S+\s+)(\S+)\s*$/ ){
		die "error in line $_\n"
	}
	print "$1" . $vals[int(rand(@vals))] . "\n";
}
