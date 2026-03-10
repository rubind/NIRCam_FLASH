for NRCFILT in F090W F200W F335M F444W
do
    for NTERMS in 0 6
    do
	python make_flat.py photo_unflat.txt $NRCFILT 10000 $NTERMS
    done
done

