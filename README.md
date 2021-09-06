# Programmable Oblivious PRF & multi-party PSI
This is a the implementation of our CCS 2021 paper: **Efficient Linear Multiparty PSI and Extensions to Circuit/Quorum PSI** [[ePrint](https://eprint.iacr.org/2021/172)], using table OPPRF as implemented in: [[MultipartyPSI](https://github.com/osu-crypto/MultiPartyPSI)] and arithmetic circuit as implemented in: [[MPCHonestMajority](https://github.com/cryptobiu/MPC-Benchmark/tree/master/MPCHonestMajority)]. The implementation is built almost entirely on the implementation in (https://github.com/osu-crypto/MultiPartyPSI).  

## Installations

### Required libraries
 C++ compiler with C++14 support. There are several library dependencies including [`Boost`](https://sourceforge.net/projects/boost/), [`Miracl`](https://github.com/miracl/MIRACL), [`NTL`](http://www.shoup.net/ntl/) , and [`libOTe`](https://github.com/osu-crypto/libOTe). For `libOTe`, it requires CPU supporting `PCLMUL`, `AES-NI`, and `SSE4.1`. Optional: `nasm` for improved SHA1 performance.   Our code has been tested on both Windows (Microsoft Visual Studio) and Linux. To install the required libraries:
  * windows: open PowerShell,  `cd ./thirdparty`, and `.\all_win.ps1` (the script works with Visual Studio 2015. For other version, you should modify [`MSBuild`](https://github.com/osu-crypto/MultipartyPSI/blob/implement/thirdparty/win/getNTL.ps1#L3) at several places in the script.)
  * linux: `cd ./thirdparty`, and `bash .\all_linux.get`.   
Also install [`libscapi`](https://github.com/cryptobiu/libscapi).

NOTE: If you meet problem with `all_win.ps1` or `all_linux.get` which builds boost, miracl and libOTe, please follow the more manual instructions at [`libOTe`](https://github.com/osu-crypto/libOTe)

### Building the Project
After cloning project from git,
##### Windows:
1. build cryptoTools,libOTe, and libOPRF projects in order.
2. add argument for bOPRFmain project (for example: -u)
3. run bOPRFmain project

##### Linux:
make (requirements: `CMake`, `Make`, `g++` or similar)


## Running the code
The database is generated randomly. The outputs include the average online/offline/total runtime that displayed on the screen and output.txt.
#### Flags:
    -u		unit test which computes PSI of 5 paries, 2 dishonestly colluding, each with set size 2^12 in semihonest setting
	-n		number of parties
	-p		party ID
	-m		set size as a power of 2 (i.e for set size 2^12 = 4096 enter 12)
	-t		number of corrupted parties (in semihonest setting)
	-a		run in augmented semihonest model. Table-based OPPRF is by default.
				0: Table-based; 1: POLY-seperated; 2-POLY-combined; 3-BloomFilter
	-r		optimized 3PSI when r = 1
	-F	output file			
#### Examples:
##### 1. Unit test:
	./bin/frontend.exe -u

##### 2. nPSI:
Compute PSI of 5 parties, 2 dishonestly colluding, each with set size 2^12 in semihonest setting

	./bin/frontend.exe -n 5 -t 2 -m 12 -p 0 -F output.txt
	& ./bin/frontend.exe -n 5 -t 2 -m 12 -p 1 -F output.txt
	& ./bin/frontend.exe -n 5 -t 2 -m 12 -p 2 -F output.txt
	& ./bin/frontend.exe -n 5 -t 2 -m 12 -p 3 -F output.txt
	& ./bin/frontend.exe -n 5 -t 2 -m 12 -p 4 -F output.txt

Alternatively, you can run

	./run\_protocol.sh 0 4 5 2 12 output.txt

## Summary

      1. git clone https://github.com/osu-crypto/MultipartyPSI.git  
      2. cd thirdparty/
      3. bash all_linux.get
      4. cd ..
      5. cmake .
      6.  make -j
      7. ./run\_protocol.sh 0 4 5 2 12 output.txt


## Help
For any questions on building or running the library, please contact the following as applicable:
	[`Ni Trieu`](http://people.oregonstate.edu/~trieun/) at trieun at oregonstate dot edu
	Akash Shah or Nishka Dasgupta at akashshah08 at outlook dot com or nishka dot dasgupta at yahoo dot com for queries regarding the circuit code
