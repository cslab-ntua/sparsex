for file in $*
do
	SPM_DELTA_SEQ_LIMIT=1000000 python2.5 python/delta_mt.py $file
done | tee ~/spm-delta-MT-results.noseq.mt__${MT_CONF}__.$(date +%Y%m%d.%H:%M:%S).$(hostname)
