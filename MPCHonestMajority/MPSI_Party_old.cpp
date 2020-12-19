#include "MPSI_Party.h"
#include "ZpKaratsubaElement.h"

template <class FieldType> MPSI_Party<FieldType>::MPSI_Party(int argc, char* argv[]) : ProtocolParty<FieldType>(argc, argv) {
	//The call to ProtocolParty constructor initializes inherited variables 
	// N, T, m_partyID, as well as the VDM matrix and related vectors.

	//Initialize global variables that have not been inherited.
        CmdParser parser = this->getParser();

        this->num_bins = stoi(parser.getValueByKey(this->arguments, "numBins"));
        this->myInputFile = parser.getValueByKey(this->arguments, "inputsFile");
        this->myOutputFile = parser.getValueByKey(this->arguments, "outputsFile");

	masks.resize(num_bins);
	a_vals.resize(num_bins);
	add_a.resize(num_bins);
	mult_outs.resize(num_bins);
	outputs.resize(num_bins);

	cout << "mpsi_print" << mpsi_print;

        readMPSIInputs();

	//Generation of shared values, such as triples, must be done later.
}

//read num_bins MPSI inputs
template <class FieldType> void MPSI_Party<FieldType>::readMPSIInputs() {
	ifstream myfile;
	int input;
	int i = 0;
	myfile.open(myInputFile);
	do {
		myfile >> input;
		add_a[i] = this->field->GetElement(input);
		i++;
	} while(!(myfile.eof()));
	myfile.close();

	if (mpsi_print == true) {
		cout << "Party " << this->m_partyID << "has read inputs: ";
		for(i=0; i<this->num_bins; i++) {
			cout << add_a[i] << " ";
		}
		cout << "\n";
	}
}

//perform MPSI
template <class FieldType> void MPSI_Party<FieldType>::runMPSI() {
	this->honestMult->invokeOffline();
	
	//Generate additive and T-sharings
	modDoubleRandom(this->num_bins, this->randomTAndAddShares);

	//Generate T-sharings
	//this->generateRandomShares(this->num_bins, this->masks);

	//Generate T and 2T sharings for multiplication
	//this->generateRandom2TAndTShares(this->num_bins, this->randomTAnd2TShares);

	//Evaluate the circuit
	//evaluateCircuit();
}

//prepare additive and T-threshold sharings of secret random value r_j using DN07's protocol
template <class FieldType> void MPSI_Party<FieldType>::modDoubleRandom(int no_random, vector<FieldType>& randomElementsToFill) {
	int index = 0;
	int N = this->N;
	int T = this->T; 
	TemplateField<FieldType> &field = this->field;
	
	vector<FieldType> x1(N), y1(N), y2(N), t1(N), r1(N), t2(N), r2(N);

	vector<vector<FieldType>> sendBufsElements(N);

	vector<vector<byte>> sendBufsBytes(N);
	vector<vector<byte>> recBufsBytes(N);
	// the number of buckets (each bucket requires one double-sharing
	// from each party and gives N-2T random double-sharings)
	int no_buckets = (no_random / (N-T))+1;

	//maybe add some elements if a partial bucket is needed
	randomElementsToFill.resize(no_buckets*(N-T)*2);

	for(int i=0; i < N; i++)
	{
		sendBufsElements[i].resize(no_buckets*2);
		sendBufsBytes[i].resize(no_buckets*field->getElementSizeInBytes()*2);
		recBufsBytes[i].resize(no_buckets*field->getElementSizeInBytes()*2);
	}

	/**
	 *  generate random sharings.
	 *  first degree T, then additive
	 *
	 */
	for(int k=0; k < no_buckets; k++)
	{
		// generate random degree-T polynomial
		for(int i = 0; i < T+1; i++)
 		{
			// A random field element, uniform distribution, 
			// note that x1[0] is the secret which is also random
			x1[i] = field->Random();
		}

		this->matrix_vand.MatrixMult(x1, y1,T+1); // eval poly at alpha-positions

		y2[0] = x1[0];
		// generate N-1 random elements
		for(int i = 1; i < N; i++)
		{
			// A random field element, uniform distribution
			y2[i] = field->Random();
			//all y2[i] generated so far are additive shares of the secret x1[0]
			y2[0] = y2[0] - y2[i];
		}

		// prepare shares to be sent
		for(int i=0; i < N; i++)
		{
			//cout << "y1[ " <<i<< "]" <<y1[i] << endl;
			sendBufsElements[i][2*k] = y1[i];
			sendBufsElements[i][2*k + 1] = y2[i];
		}
	}

	int fieldByteSize = field->getElementSizeInBytes();
	for(int i=0; i < N; i++)
	{
		for(int j=0; j<sendBufsElements[i].size();j++) {
			field->elementToBytes(sendBufsBytes[i].data() + (j * fieldByteSize), sendBufsElements[i][j]);
		}
	}

	this->roundFunctionSync(sendBufsBytes, recBufsBytes, 1);

	for(int k=0; k < no_buckets; k++) {
		for (int i = 0; i < N; i++) {
			t1[i] = field->bytesToElement(recBufsBytes[i].data() + (2*k * fieldByteSize));
			t2[i] = field->bytesToElement(recBufsBytes[i].data() + ((2*k +1) * fieldByteSize));
		}
		this->matrix_vand_transpose.MatrixMult(t1, r1,N-T);
		this->matrix_vand_transpose.MatrixMult(t2, r2,N-T);

	        //copy the resulting vector to the array of randoms
		for (int i = 0; i < (N - T); i++) {
			randomElementsToFill[index*2] = r1[i];
			randomElementsToFill[index*2 +1] = r2[i];
			index++;
		}
	}

	if (mpsi_print == true) {
		cout << "Party " << this->m_partyID << "has shares: ";
		cout << randomElementsToFill[0] << " " << randomElementsToFill[1] << "\n";
	}
}

//reshare values in parameter as t-sharings
//save own share in a_vals
template <class FieldType> void MPSI_Party<FieldType>::reshare(vector<FieldType>& vals) {
	int index = 0;
	int N = this->N;
	int T = this->T;
	int no_vals = this->num_bins;

	vector<vector<byte>> recBufsBytes(N);
	
	TemplateField<FieldType> &field = this->field;

	vector<FieldType> x1(N), y1(N);

	vector<vector<FieldType>> sendBufsElements(N);
	vector<vector<byte>> sendBufsBytes(N);

	int fieldByteSize = field->getElementSizeInBytes();

	for(int i=0; i < N; i++)
	{
		sendBufsElements[i].resize(no_vals);
		sendBufsBytes[i].resize(no_vals*fieldByteSize);
		recBufsBytes[i].resize(no_vals*fieldByteSize);
	}

	if(this->m_partyID == 0) {
		//generate T-sharings of the values in vals
		for(int k = 0; k < no_vals; k++)
		{
			//set x1[0] as the secret to be shared
			x1[0] = vals[k];
			// generate random degree-T polynomial
			for(int i = 1; i < T+1; i++)
			{
				// A random field element, uniform distribution
				x1[i] = field->Random();
			}

			this->matrix_vand.MatrixMult(x1, y1,T+1); // eval poly at alpha-positions

			// prepare shares to be sent
			for(int i=0; i < N; i++)
			{
				//cout << "y1[ " <<i<< "]" <<y1[i] << endl;
				sendBufsElements[i][k] = y1[i];
			}
		
			//this->a_vals[k] = y1[0];
		}
	}

	for(int i=0; i < N; i++)
	{
		for(int j=0; j<sendBufsElements[i].size();j++) {
		field->elementToBytes(sendBufsBytes[i].data() + (j * fieldByteSize), sendBufsElements[i][j]);
		}
	}

	this->roundFunctionSync(sendBufsBytes, recBufsBytes, 2);
	
	for(int k=0; k < no_vals; k++) {
		this->a_vals[k] = field->bytesToElement(recBufsBytes[0].data() + (k * fieldByteSize));
        }

	if (mpsi_print == true) {
		cout << "Party " << this->m_partyID << " has received t-sharings: ";
		for(int k=0; k < no_vals; k++) {
			cout << this->a_vals[k];
		}
		cout << "\n";
	}

}

template <class FieldType> void MPSI_Party<FieldType>::evaluateCircuit() {
	add_rj();
	//subtract_rj();
	//mult_sj();
	//if(this->m_partyID == 0) {
	//	this->openShare(this->num_bins, this->mult_outs, this->outputs);
	//	outputPrint();
	//}
}

template <class FieldType> void MPSI_Party<FieldType>::add_rj() {
	int j;
	vector<FieldType> reconar; // reconstructed aj+rj
	reconar.resize(num_bins);

	//add additive share of rj to corresponding share of aj
	for(j=0; j<num_bins; j++) {
		add_a[j] = add_a[j] + this->randomTAnd2TShares[j*2+1];
	}

	//reconstruct additive shares, store in reconar
	addShareOpen(num_bins, add_a, reconar);	
	
	//reshare and save in a_vals;
	reshare(reconar);
}

template <class FieldType> void MPSI_Party<FieldType>::subtract_rj() {
	int j;

	for(j=0; j<num_bins; j++) {
		a_vals[j] = a_vals[j] - this->randomTAnd2TShares[j*2];
	}
}

template <class FieldType> void MPSI_Party<FieldType>::mult_sj() {
	this->DNHonestMultiplication(this->masks, this->a_vals, this->mult_outs, this->num_bins); 
}

//open an additive sharing
//code identical to openShare() except for final reconstruction step
template <class FieldType>
void MPSI_Party<FieldType>::addShareOpen(int numShares, vector<FieldType> &Shares, vector<FieldType> &secrets) {
	int N = this->N;
	int T = this->T;
	TemplateField<FieldType>& field = this->field;

	vector<vector<byte>> sendBufsBytes(N);
	vector<vector<byte>> recBufsBytes(N);
	vector<FieldType> x1(N);
	int fieldByteSize = field->getElementSizeInBytes();

	//resize vectors
	for(int i=0; i < N; i++)
	{
		sendBufsBytes[i].resize(numShares*fieldByteSize);
		recBufsBytes[i].resize(numShares*fieldByteSize);
	}

	//set the first sending data buffer
	for(int j=0; j<numShares;j++) {
		field->elementToBytes(sendBufsBytes[0].data() + (j * fieldByteSize), Shares[j]);
	}

	//copy the same data for all parties
	for(int i=1; i<N; i++){
		sendBufsBytes[i] = sendBufsBytes[0];
	}

	//call the round function to send the shares to all the users and get the other parties share
	this->roundFunctionSync(sendBufsBytes, recBufsBytes,12);

	//reconstruct each set of shares to get the secret
	//only leader does this step
	if(m_partyId == 0) {
		for(int k=0; k<numShares; k++){
			//get the set of shares for each element
			for(int i=0; i < N; i++) {
				x1[i] = field->bytesToElement(recBufsBytes[i].data() + (k*fieldByteSize));
			}
		
			//reconstruct the shares
			secrets[k] = *(field->GetZero());
			for(int i=0; i < N; i++) {
				secrets[k] = secrets[k] + x1[i];
			}
      		}
	}
}

//print output results to file
template <class FieldType> void MPSI_Party<FieldType>::outputPrint() {
	vector<int> matches;
	TemplateField<FieldType> &field = this->field;
	ofstream myfile; 
	int counter=0;
	int i;

	for(i=0; i < this->num_bins; i++) {
		if(outputs[i] == *field->GetZero()) {
			matches.push_back(i);
			counter++;
		}
	}
	
	myfile.open(myOutputFile);
	for(i=0; i < counter; i++) {
		myfile << matches[i] << "\n";
	}
	myfile.close();
}

template <class FieldType> MPSI_Party<FieldType>::~MPSI_Party() {
	
}
