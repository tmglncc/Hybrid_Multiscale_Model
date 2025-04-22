#include "MarkovChain.hpp"

using namespace std;

MarkovChain::MarkovChain(int N) {
	this->N = N;

	this->initial = new float[this->N];
	this->transition = new float*[this->N];
	for (int i = 0; i < this->N; i++) {
		this->transition[i] = new float[this->N];

		this->initial[i] = 0.5;
		for (int j = 0; j < this->N; j++)
			this->transition[i][j] = 0.5;
	}
}

MarkovChain::MarkovChain() {
	this->N = 0;
	this->initial = NULL;
	this->transition = NULL;
}

MarkovChain::~MarkovChain() {
	if (this->transition != NULL) {
		for (int i = 0; i < this->N; i++)
			delete this->transition[i];
		delete this->transition;
	}

	if (this->initial != NULL)
		delete this->initial;

	this->N = 0;
}

bool MarkovChain::isStochasticMatrix() {
	float sum = 0.0;
	for (int i = 0; i < this->N; i++)
		sum += this->initial[i];
	if (sum != 1.0)
		return false;

	for (int i = 0; i < this->N; i++) {
		sum = 0.0;
		for (int j = 0; j < this->N; j++)
			sum += this->transition[i][j];
		if (sum != 1.0)
			return false;
	}

	return true;
}

/**
 * @brief setModel : function to set the parameters of the model
 * @param transition : NXN transition matrix
 * @param initial : 1XN initial probability matrix
 */
void MarkovChain::setModel(int N, float* initial, float* transition) {
	this->N = N;

	this->initial = new float[this->N];
	this->transition = new float*[this->N];
	for (int i = 0; i < this->N; i++) {
		this->transition[i] = new float[this->N];

		this->initial[i] = initial[i];
		for (int j = 0; j < this->N; j++)
			this->transition[i][j] = transition[i*this->N + j];
	}
}

/**
 * @brief computeProbability : compute the probability that the sequence
 *        is generate from the markov chain
 * @param sequence : is a vector of integral sequence starting from 0
 * @return         : is probability
 */
float MarkovChain::computeProbability(vector<int> sequence) {
	float prob = this->initial[sequence[0]];
	for (int i = 0; i < sequence.size()-1; i++)
		prob *= this->transition[sequence[i]][sequence[i+1]];
	
	return prob;
}

/**
 * @brief initialRand : function to generate a radom state
 * @param matrix      : input matrix
 * @param index       ; row of matrix to consider
 * @return
 */
int MarkovChain::generateInitialState(Ran *ran) {
	float randomValue = (float) ran->doub();
	float prob = this->initial[0];
	int state = 0;

	//select the index corresponding to the highest probability
	//or if all the cols of matrix have transitioned
	while (randomValue > prob && state < this->N) {
		state++;
		prob += this->initial[state];
	}

	return state;
}

/**
 * @brief initialRand : function to generate a radom state
 * @param matrix      : input matrix
 * @param index       ; row of matrix to consider
 * @return
 */
int MarkovChain::generateCurrentState(Ran *ran, int previousState) {
	float randomValue = (float) ran->doub();
	float prob = this->transition[previousState][0];
	int state = 0;

	//select the index corresponding to the highest probability
	//or if all the cols of matrix have transitioned
	while (randomValue > prob && state < this->N) {
		state++;
		prob += this->transition[previousState][state];
	}

	return state;
}

/**
 * @brief generateSequence is a function that generates a sequence of specified length
 * @param n : is the length of the sequence
 * @return : is a vector of integers representing generated sequence
 */
vector<int> MarkovChain::generateSequence(Ran *ran, int time) {
	vector<int> sequence(time);

	//select a random initial value of sequence
	sequence[0] = generateInitialState(ran);
	for (int i = 1; i < time; i++) {
		//select a random transition to next sequence state
		sequence[i] = generateCurrentState(ran, sequence[i-1]);
	}

	return sequence;
}

void MarkovChain::printSequence(vector<int> sequence) {
	cout << "[";
	for (int i = 0; i < sequence.size(); i++)
		cout << sequence[i] << "\t";
	cout << "]" << endl;
}
