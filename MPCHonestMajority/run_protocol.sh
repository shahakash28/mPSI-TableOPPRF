#! /bin/bash
for i in `seq $2 -1 $1`;
do	
	./build/MPCHonestMajority -partyID $i -partiesNumber $3 -numBins $4 -inputsFile ../../in_party_$i.txt -outputsFile output.txt -circuitFile ic.txt -fieldType $5 -genRandomSharesType HIM -multType DN -verifyType Single -partiesFile Parties.txt -internalIterationsNumber 1 > test/$3/$4/op_party_$i.txt&
	echo "Running $i..." &
done
