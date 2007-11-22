for type in "double" #"float"
do
	for file in $*
	do
		#echo python2.5 crs.py $file $type 32
		python2.5 crs.py $file $type 32
		#echo python2.5 crs.py $file $type 64
		python2.5 crs.py $file $type 64
		for i in $(seq 2 13)
		do
			#echo python2.5 delta_parse.py $file $type $i
			python2.5 delta_parse.py $file $type $i
		done
		#echo python2.5 delta_parse.py $file $type 100000
		python2.5 delta_parse.py $file $type 100000
	done
done | tee delta-$(hostname)-$(date -I)
