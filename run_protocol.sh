#! /bin/bash
for i in `seq $2 -1 $1`;
do
	./bin/frontend.exe -n $3 -t $4 -m $5 -p $i -F $6 &
	echo "Running $i..." &
done

