/***************************************************
channel coding course work: conolutional codes
this program template has given the message generator, bpsk modulation, awgn channel model and bpsk demodulation,
you should first determine the encoder structure, then define the message and codeword length, generate the state table, write the convolutional encoder and decoder.

if you have any question, please contact me via e-mail: xingjyue@mail2.sysu.edu.cn
***************************************************/

#define  _crt_secure_no_warnings
//#define _CRT_SECURE_NO_WARNINGS
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <limits>

using namespace std;

#define infinity INT_MAX
#define in_bit 1
#define out_bit 2
#define state_num 4 //the number of state
#define constraint_length 3 // the constraint length of conv. code
#define memory_num 2 //the number of memory units (constraint length substract 1)

#define message_length (int)1e4 //the length of message (frame)
#define codeword_length (int)2e4 //the length of codeword

// Code Vector
// (7,5)octal conv code
int generator_table[out_bit][constraint_length] = {1,1,1,1,0,1 };
//{ 1,1,1,1,0,0,1,
//1,0,1,1,0,1,1 };
float code_rate = (float)message_length / (float)codeword_length;
int branch_num = (int)pow(2, in_bit)*state_num;
// input num is not equal to the input bit num
int input_num = (int)pow(2, in_bit);
//determine the decoding algorithm
//0 for hard decision viterbi, 1 for soft decision viterbi, 2 for BCJR algorithm, 3 for uncoded system
char type = '2';

// channel coefficient
#define pi 3.1415926
double n0, sgm;

//state table, the size should be defined yourself, the column of table :input,current state,next state,output
//int state_table[2 * state_num][in_bit+out_bit+2];
int **state_table;
double **branch_table = new double*[branch_num];

//int branch_table[2 * state_num][message_length];//branch metrics table
//int path_table[state_num][message_length + 1];//path metrics table
//double branch_table[2 * state_num][message_length];//branch metrics table

double path_table[state_num][message_length/in_bit + 1];//path metrics table
int parentID_table[state_num][message_length/in_bit];//trellis transition id table

int message[message_length], codeword[codeword_length];//message and codeword
int re_codeword[codeword_length];//the received codeword
int de_message[message_length];//the decoding message

double tx_symbol[codeword_length][2];//the transmitted symbols
double rx_symbol[codeword_length][2];//the received symbols

int bin2dec(int[], int);
int* dec2bin(int, int);

void statetable();
void encoder();
void modulation();
void demodulation();
void channel();
double gaussrand(double,double);
void branchtable();
void pathtable();
void decoder(double);
void bcjr_decoder(double);

void main()
{
	// gateway function
	int i;
	float snr, start, finish;
	long int bit_error, seq, seq_num;
	double ber;
	double progress;

	cout << "generator table" << endl;
	for (int i = 0; i < out_bit; i++)
	{
		for (int j = 0; j < constraint_length; j++)
			cout << generator_table[i][j] << " ";
		cout << endl;
	}

	switch (type)
	{
	case '0':
		for (int branch = 0; branch < branch_num; branch++)
			branch_table[branch] = new double[message_length / in_bit];
		break;
	case '1':
		for (int branch = 0; branch < branch_num; branch++)
			branch_table[branch] = new double[message_length / in_bit];
		break;
	default:break;
	}

	//generate state table
	statetable();
	//cout << "statetable function is good!" << endl;
	//random seed
	srand((int)time(0));

	//input the snr and frame number
	printf("\nenter start snr: ");
	scanf_s("%f", &start);
	printf("\nenter finish snr: ");
	scanf_s("%f", &finish);
	printf("\nplease input the frame number: ");
	scanf_s("%d", &seq_num);

	ofstream file;
	file.open("(7,5)_bcjr.csv");
	//file.open("(7,5)_SoftViterbi.csv");
	//file.open("(7,5)_HardViterbi.csv");
	//file.open("(7,5)_uncode.csv");
	//file.open("(15,13)_bcjr.csv");
	//file.open("(23,35)_bcjr.csv");
	//file.open("(171,133)_bcjr.csv");
	//file.open("(7,7,5)_bcjr.csv");
	//file.open("(7,5)_uncode.csv");

	for (snr = start; snr <= finish; snr += 1)
	{
		//channel noise
		n0 = (1.0 / code_rate) / pow(10.0, (float)(snr) / 10.0);
		sgm = sqrt(n0/2);

		bit_error = 0;

		for (seq = 1; seq <= seq_num; seq++)
		{
			if (type == '3')
			{// type = '3' (uncode)
				if (message_length != codeword_length)
					cout << "Error!The message length is not equal to the codeword length!" << endl;
				for (i = 0; i < message_length; i++)
					message[i] = rand() % 2;
				for (int i = 0; i < message_length; i++)
					codeword[i] = message[i];
				modulation(); channel(); demodulation();
				for (int i = 0; i < codeword_length; i++)
					de_message[i] = re_codeword[i];
			}
			else {
				// pay attention that message is appended by 0 whose number is equal to the memory number (constraint length - 1) of encoder structure.
				for (i = 0; i < message_length - memory_num*in_bit; i++)
					message[i] = rand() % 2;
				for (i = message_length - memory_num*in_bit; i < message_length; i++)
					message[i] = 0;
				//convolutional encoder
				encoder();
				//cout << "encoder functikon is good!" << endl;
				//bpsk modulation
				modulation();
				//awgn channel
				channel();
				//bpsk demodulation, it's needed in hard-decision viterbi decoder
				if (type == '0') demodulation();

				//convolutional decoder
				decoder(n0);
			}
			//calculate the number of bit error
			for (i = 0; i < message_length; i++)
			{
				if (message[i] != de_message[i])
					bit_error++;
			}
			progress = (double)(seq * 100) / (double)seq_num;
			//calculate the ber
			//ber = (double)bit_error / (double)(message_length*seq);
			//print the intermediate result
			//printf("progress=%2.1f, snr=%2.1f, bit errors=%2.1d, ber=%e\r", progress, snr, bit_error, ber);
		}
		//calculate the ber
		ber = (double)bit_error / (double)(message_length*seq_num);
		//print the final result
		printf("progress=%2.1f, snr=%2.1f, bit errors=%2.1d, ber=%e\n", progress, snr, bit_error, ber);
		/*cout << "transmitted message" << endl;
		for (int i = 0; i < message_length; i++)
			cout << message[i] << " ";
		cout << endl;
		cout << "transmitted codeword" << endl;
		for (int i = 0; i < codeword_length; i++)
			cout << codeword[i] << " ";
		cout << endl;
		cout << "received codeword" << endl;
		for (int i = 0; i < codeword_length; i++)
			cout << re_codeword[i] << " ";
		cout << endl;
		cout << "demodulated message" << endl;
		for (int i = 0; i < message_length; i++)
			cout << de_message[i] << " ";
		cout << endl;*/
		file << snr << "," << ber << endl;
	}
	file.close();
	// delete the point
	if (type == '0'||type=='1')
	{
		for (int branch = 0; branch < branch_num; branch++)
		{
			delete[] branch_table[branch];
			delete[] state_table[branch];
		}
		delete[] branch_table;
		delete[] state_table;
	}
	system("pause");
}
int* dec2bin(int dec_val, int size)
{
	int count = 0, current = dec_val;
	int *p = new int[size];
	int *q = new int[size];
	for (int i = 0; i < size; i++)
		*(p + i) = 0;
	while (current != 0)
	{
		p[count] = current % 2;
		current = current / 2;
		count++;
	}
	// reverse
	for (int i = 0; i < size; i++)
		*(q + i) = *(p + size - 1 - i);
	return q;
}
void statetable()
{
	state_table = new int*[branch_num];
	for (int branch = 0; branch < branch_num; branch++)
		state_table[branch] = new int[in_bit + out_bit + 2];
	//generator polynomial 
	/*g[0][0] = 1; g[0][1] = 1; g[0][2] = 1;
	g[1][0] = 1; g[1][1] = 0; g[1][2] = 1;*/
	for (int i = 0; i < state_num; i++)
	{
		//int current_state[memory_num*in_bit];//binary representation of current state
		int next_state[memory_num*in_bit];//binary representation of next state
		// convert decimal to binary
		int *current_state = dec2bin(i, memory_num*in_bit);
		// determine the last (K-1)*k, shift k input message into register
		for (int j = memory_num * in_bit - 1; j > in_bit - 1; j--)
			next_state[j] = current_state[j - in_bit];
		// input message
		for (int index = 0; index < input_num; index++)
		{
			for (int in_offset = 0; in_offset < in_bit; in_offset++)
			{
				int *temp = dec2bin(index, in_bit);
				state_table[input_num * i + index][in_offset] = *(temp + in_offset);
			}
			// current state
			state_table[input_num * i + index][in_bit] = i;
			//next state
			for (int j = 0; j < in_bit; j++)
				next_state[j] = state_table[input_num * i + index][j];
			state_table[input_num*i + index][in_bit + 1] = bin2dec(next_state, memory_num*in_bit);
			// output codeword
			for (int out_offset = 0; out_offset < out_bit; out_offset++)
			{
				int temp = 0;
				for (int in_offset = 0; in_offset < in_bit; in_offset++)
					temp += state_table[input_num * i + index][in_offset] * generator_table[out_offset][in_offset];
				for (int state_offset = 0; state_offset < memory_num*in_bit; state_offset++)
					temp += current_state[state_offset] * generator_table[out_offset][state_offset + in_bit];
				state_table[input_num * i + index][out_offset + in_bit + 2] = temp % 2;
			}
		}
	}
	cout << "state table" << endl;
	cout << right << setw(in_bit * 7) << "input" << " ";
	cout << right << setw(memory_num * 7) << "current state" << " ";
	cout << right << setw(memory_num * 7) << "next state" << " ";
	cout << right << setw(out_bit * 7) << "output" << " ";
	cout << right << setw(3) << "ID" << endl;
	for (int i = 0; i < (in_bit * 7 + memory_num * 14 + out_bit * 7 + 7); i++)
		cout << "-";
	cout << endl;
	for (int i = 0; i < branch_num; i++)
	{
		// input
		for (int j = 0; j < in_bit; j++)
			cout << right << setw(7) << state_table[i][j];
		cout << "|";
		// current state (binary representation)
		for (int j = 0; j < memory_num*in_bit; j++)
			cout << right << setw(7) << *(dec2bin(state_table[i][in_bit], memory_num) + j);
		cout << "|";
		// next state (binary representation)
		for (int j = 0; j < memory_num*in_bit; j++)
			cout << right << setw(7) << *(dec2bin(state_table[i][in_bit + 1], memory_num) + j);
		cout << "|";
		// output
		for (int j = in_bit + 2; j < in_bit + out_bit + 2; j++)
			cout << right << setw(7) << state_table[i][j];
		cout << "|";
		cout << right << setw(3) << i << endl;
	}
	for (int i = 0; i < (in_bit * 7 + memory_num * 14 + out_bit * 7 + 7); i++) cout << "-";
	cout << endl;
}
int bin2dec(int binary[], int length)
{
	int demical_value = 0;
	for (int i = 0; i < length; i++)
	{
		demical_value += (int)(binary[length - 1 - i] * pow(2, i));
	}
	return demical_value;
}

void encoder()
{
	//convolution encoder, the input is message[] and the output is codeword[]
	// generator matrix initialization
	/*int G[message_length][codeword_length];
	for (int i = 0; i < message_length; i++)
	{
		for (int j = 0; j < codeword_length; j++)
		{
			G[i][j] = 0;
			G[i][2 * i] = 1; G[i][2 * i + 1] = 1;
			if (2 * i + 5 < codeword_length)
			{
				G[i][2 * i + 2] = 1; G[i][2 * i + 3] = 0;
				G[i][2 * i + 4] = 1; G[i][2 * i + 5] = 1;
			}
			else if (2 * i + 3 < codeword_length)
			{
				G[i][2 * i + 2] = 1; G[i][2 * i + 3] = 0;
			}
		}
	}
	for (int k = 0; k < codeword_length; k++)
	{
		int tmp = 0;
		for (int i = 0; i < message_length; i++)
		{
			tmp += message[i] * G[i][k];
		}
		codeword[k] = tmp % 2;
	}*/
	// convolution operation
	for (int msg = in_bit - 1; msg < message_length; msg += in_bit)
		for (int offset = 0; offset < out_bit; offset++)
		{
			int tmp = 0;
			for (int reg = 0; reg < constraint_length*in_bit; reg++)
			{
				if (msg >= reg)
					tmp += generator_table[offset][reg] * message[msg - reg];
			}
			codeword[msg*out_bit + offset] = tmp % 2;
		}
}

void modulation()
{
	//bpsk modulation
	int i;

	//0 is mapped to (1,0) and 1 is mapped tp (-1,0)
	for (i = 0; i < codeword_length; i++)
	{
		tx_symbol[i][0] = -1 * (2 * codeword[i] - 1);
		tx_symbol[i][1] = 0;
	}
}
void channel()
{
	//awgn channel
	int i, j;
	double u, r, g;

	for (i = 0; i < codeword_length; i++)
	{
		for (j = 0; j < 2; j++)
		{
			u = (float)rand() / (float)RAND_MAX;
			if (u == 1.0)
				u = 0.999999;
			r = sgm*sqrt(2.0*log(1.0 / (1.0 - u)));

			u = (float)rand() / (float)RAND_MAX;
			if (u == 1.0)
				u = 0.999999;
			g = (float)r*cos(2 * pi*u);

			rx_symbol[i][j] = tx_symbol[i][j] + g;
		}
		/*rx_symbol[i][0] = tx_symbol[i][0] + gaussrand(0, sgm);
		rx_symbol[i][1] = tx_symbol[i][1];*/
	}
}
double gaussrand(double mu, double sigma)
{
	const double epsilon = std::numeric_limits<double>::min();
	const double two_pi = 2.0*3.14159265358979323846;

	static double z0, z1;
	static bool generate;
	generate = !generate;

	if (!generate)
		return z1 * sigma + mu;

	double u1, u2;
	do
	{
		u1 = rand() * (1.0 / RAND_MAX);
		u2 = rand() * (1.0 / RAND_MAX);
	} while (u1 <= epsilon);

	z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
	z1 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);
	return z0 * sigma + mu;
}
void demodulation()
{
	int i;
	double d1, d2;
	for (i = 0; i < codeword_length; i++)
	{
		d1 = (rx_symbol[i][0] - 1)*(rx_symbol[i][0] - 1) + rx_symbol[i][1] * rx_symbol[i][1];
		d2 = (rx_symbol[i][0] + 1)*(rx_symbol[i][0] + 1) + rx_symbol[i][1] * rx_symbol[i][1];
		if (d1 < d2)
			re_codeword[i] = 0;
		else
			re_codeword[i] = 1;
	}
}
void branchtable()
{
	bool **flag = new bool*[branch_num];
	for (int branch = 0; branch < branch_num; branch++)
		flag[branch] = new bool[message_length / in_bit];
	// initialization 
	for (int r = 0; r < branch_num; r++)
	{
		for (int c = 0; c < message_length / in_bit; c++)
		{
			if (c < memory_num || c >= message_length / in_bit - memory_num)
				flag[r][c] = false;
			else flag[r][c] = true;
		}
	}
	for (int j = 0; j < message_length / in_bit; j++)
	{
		for (int i = 0; i < branch_num; i++)
		{
			branch_table[i][j] = infinity;
		}
	}
	// initial condition
	//flag[0][0] = true; flag[1][0] = true;
	for (int i = 0; i < input_num; i++)
		flag[i][0] = true;
	for (int i = 0; i < branch_num; i += input_num)
		flag[i][message_length / in_bit - memory_num] = true;
	for (int col = 1; col < message_length / in_bit; col++)
	{
		if (col < memory_num)
		{
			// search the state table
			for (int row = 0; row < branch_num; row++)
			{
				if (flag[row][col - 1])
				{
					// next state
					int nextstate = state_table[row][in_bit + 1];
					for (int id = 0; id <branch_num; id++)
					{
						// find out whose starting state is nextstate
						if (state_table[id][in_bit] == nextstate)
							flag[id][col] = true;
					}
				}
			}
		}
		else if (col > message_length / in_bit - memory_num)
		{
			// search the state table
			for (int row = 0; row <branch_num; row++)
			{
				if (flag[row][col - 1])
				{
					// next state
					int nextstate = state_table[row][in_bit + 1];
					for (int id = 0; id < branch_num; id++)
					{
						// find out whose starting state is nextstate
						// pay attention to judging if id is even or odd
						if (state_table[id][in_bit] == nextstate && id % input_num == 0)
							flag[id][col] = true;
					}
				}
			}
		}
	}
	for (int column = 0; column < message_length / in_bit; column++)
	{
		for (int row = 0; row < branch_num; row++)
		{
			if (column >= memory_num && column < message_length / in_bit - memory_num)
			{
				double dist = 0;
				for (int offset = 0; offset < out_bit; offset++)
				{
					if (type == '0')
						// ? 
						dist += abs(re_codeword[out_bit* column + offset] - state_table[row][in_bit + 2 + offset]);
					else if (type == '1')
					{
						double tmp = (rx_symbol[out_bit * column + offset][0] - (-1 * (2 * state_table[row][in_bit + 2 + offset] - 1)))*(rx_symbol[out_bit* column + offset][0] - (-1 * (2 * state_table[row][in_bit + 2 + offset] - 1)))
							+ (rx_symbol[out_bit * column + offset][1] - 0)*(rx_symbol[out_bit* column + offset][1] - 0);
						//dist += sqrt(tmp);
						dist += tmp;
					}
				}
				branch_table[row][column] = dist;
				//branch_table[row][column] = sqrt(dist);
			}
			else
			{
				if (flag[row][column])
				{
					double dist = 0;
					for (int offset = 0; offset < out_bit; offset++)
					{
						if (type == '0')
							dist += abs(re_codeword[out_bit * column + offset] - state_table[row][in_bit + 2 + offset]);
						else if (type == '1')
						{
							double tmp = (rx_symbol[out_bit * column + offset][0] - (-1 * (2 * state_table[row][in_bit + 2 + offset] - 1)))*(rx_symbol[out_bit* column + offset][0] - (-1 * (2 * state_table[row][in_bit + 2 + offset] - 1)))
								+ (rx_symbol[out_bit * column + offset][1] - 0)*(rx_symbol[out_bit* column + offset][1] - 0);
							//dist += sqrt(tmp);
							dist += tmp;
						}
					}
					branch_table[row][column] = dist;
					//branch_table[row][column] = sqrt(dist);
				}
			}
		}
	}
	/*cout << "branch table" << endl;
	for (int i = 0; i < branch_num; i++)
	{
		for (int j = 0; j < message_length/in_bit; j++)
		{
			if (branch_table[i][j] != infinity)
				cout << right << setw(8) << branch_table[i][j] << " ";
			else
				cout << right << setw(8) << "*" << " ";
		}
		cout << endl;
	}
	cout << "branch exist flag" << endl;
	for (int r = 0; r < branch_num; r++)
	{
		for (int c = 0; c < message_length/in_bit; c++)
		{
			if (flag[r][c]) cout << "-" << " ";
			else cout << "*" << " ";
		}
		cout << endl;
	}*/
	for (int branch = 0; branch < branch_num; branch++)
		delete[] flag[branch];
	delete[] flag;
}
void pathtable()
{
	// initialize path table
	for (int j = 0; j < message_length / in_bit + 1; j++)
		for (int i = 0; i < state_num; i++)
			path_table[i][j] = infinity;
	// initiaize trellis id table
	for (int i = 0; i < state_num; i++)
		for (int j = 0; j < message_length / in_bit; j++)
			parentID_table[i][j] = -1;// illegal ID
	// initial condition
	path_table[0][0] = 0;
	for (int j = 1; j < message_length / in_bit + 1; j++)
	{
		for (int i = 0; i < state_num; i++)
		{
			int pre_state;
			double min = path_table[i][j];
			//search the branch connected to i-th state in the state table
			for (int id = 0; id < branch_num; id++)
			{
				if (state_table[id][in_bit + 1] == i)// if the next state is i
				{
					pre_state = state_table[id][in_bit];
					if ((path_table[pre_state][j - 1] + branch_table[id][j - 1]) <= min)
					{
						//update the path metrics of i-th state at time j
						min = path_table[pre_state][j - 1] + branch_table[id][j - 1];
						parentID_table[i][j - 1] = id;
					}
				}
			}
			path_table[i][j] = min;
		}
	}
	//cout << "path table" << endl;
	//for (int i = 0; i < state_num; i++)
	//{
	//	for (int j = 0; j < message_length/in_bit + 1; j++)
	//	{
	//		if(path_table[i][j]!=infinity)
	//		cout << right << setw(8) << path_table[i][j] << " ";
	//		else
	//			cout << right << setw(8) << "*" << " ";
	//	}
	//	cout << endl;
	//}
	//cout << "trellis ID table" << endl;
	//for (int i = 0; i < state_num; i++)
	//{
	//	for (int j = 0; j < message_length/in_bit; j++)
	//	{
	//		cout << right << setw(2) << parentID_table[i][j] << " ";
	//	}
	//	cout << endl;
	//}
}
void decoder(double N0)
{
	// Viterbi decoding (0 for hard, 1 for soft)
	if (type == '0' || type == '1')
	{
		branchtable();
		pathtable();
		// determine the minimum path metrics at time message_length/in_bit+1
		double min = path_table[0][message_length / in_bit];
		int end_state = 0;
		// backtracking
		int branch_id;
		for (int j = message_length / in_bit - 1; j > -1; j--)
		{
			branch_id = parentID_table[end_state][j];
			for (int in_offset = 0; in_offset < in_bit; in_offset++)
				de_message[j*in_bit + in_offset] = state_table[branch_id][in_offset];
			//update the end state
			end_state = state_table[branch_id][in_bit];
		}
	}
	// BCJR decoding algorithm
	else if (type == '2')
		bcjr_decoder(N0);
}
void bcjr_decoder(double N0)
{
	double alpha[state_num][message_length / in_bit + 1];
	double beta[state_num][message_length / in_bit + 1];
	for (int i = 0; i < state_num; i++)
		for (int j = 0; j < message_length / in_bit + 1; j++)
		{
			alpha[i][j] = 0;
			beta[i][j] = 0;
		}
	double **gamma = new double*[branch_num];
	for (int branch = 0; branch < branch_num; branch++)
		gamma[branch] = new double[message_length / in_bit];
	// initial condition
	alpha[0][0] = 1; beta[0][message_length / in_bit] = 1;
	//cout << "n0: " << N0 << endl;
	// determine the state transition probabilities
	for (int column = 0; column < message_length / in_bit; column++)
	{
		for (int branch = 0; branch < branch_num; branch++)
		{
			// priori prob. of info.
			gamma[branch][column] = 0.5;
			for (int offset = 0; offset < out_bit; offset++)
			{
				double temp,temp1,temp2, prob,prob1,prob2;
				// channel observation
				temp=(rx_symbol[out_bit*column + offset][0] + (2 * state_table[branch][in_bit + 2 + offset] - 1))*(rx_symbol[out_bit*column + offset][0] + (2 * state_table[branch][in_bit + 2 + offset] - 1))
					+ (rx_symbol[out_bit*column + offset][1] - 0)*(rx_symbol[out_bit*column + offset][1] - 0);
				temp1 = (rx_symbol[out_bit*column + offset][0] - 1)*(rx_symbol[out_bit*column + offset][0] - 1)
					+ (rx_symbol[out_bit*column + offset][1] - 0)*(rx_symbol[out_bit*column + offset][1] - 0);
				temp2= (rx_symbol[out_bit*column + offset][0] +1)*(rx_symbol[out_bit*column + offset][0] + 1)
					+ (rx_symbol[out_bit*column + offset][1] - 0)*(rx_symbol[out_bit*column + offset][1] - 0);
				prob1 = (1/(pi*N0))*exp(-temp1 / N0);
				prob2= (1 / (pi*N0))*exp(-temp2 / N0);
				if (state_table[branch][in_bit + 2 + offset] == 0) prob = prob1 / (prob1 + prob2);
				else prob = prob2 / (prob1 + prob2);
				/*cout << "prob:" << prob << endl;*/
				//prob = (1 / (pi*N0))*exp(-temp / N0);
				gamma[branch][column] *= prob;
			}
		}
	}
	// determine the prob. of each begining state
	// forward recursion
	for (int column = 1; column < message_length / in_bit + 1; column++)
	{		
		// normalization factor
		double NA = 0;
		for (int cur_state = 0; cur_state < state_num; cur_state++)
		{
			for (int pre_state = 0; pre_state < state_num; pre_state++)
			{
				if (alpha[pre_state][column - 1]!=0)
				{
					int branch_id; bool flag = false;// if there exists a branch connecting pre state to cur state, then flag is set true
					 // search the state table to check if previous state node is connected to current state node
					for (int i = 0; i < input_num; i++)
						if (state_table[input_num*pre_state + i][in_bit + 1] == cur_state)
						{
							branch_id = input_num*pre_state + i; flag = true; break;
						}
					if(flag) alpha[cur_state][column] += alpha[pre_state][column - 1] * gamma[branch_id][column - 1];
					else continue;
				}
			}
			NA += alpha[cur_state][column];
		}
		// normalization process
		for (int cur_state = 0; cur_state < state_num; cur_state++)
			alpha[cur_state][column] = (1 / NA)*alpha[cur_state][column];
	}
	// determine the prob. of each ending state
	// backward recursion
	for (int column = message_length / in_bit - 1; column > -1; column--)
	{
		double NB = 0;
		for (int cur_state = 0; cur_state < state_num; cur_state++)
		{
			for (int next_state = 0; next_state < state_num; next_state++)
			{
				if (beta[next_state][column + 1]!=0)
				{
					int branch_id; bool flag = false;
					// search the state table to check if current state node is connected to next state node
					for (int i = 0; i < input_num;i++)
						if (state_table[input_num*cur_state + i][in_bit + 1] == next_state)
						{
							branch_id = input_num*cur_state + i; flag = true; break;
						}
					if (flag) beta[cur_state][column] += beta[next_state][column + 1] * gamma[branch_id][column];
					else continue;
				}
			}
			NB += beta[cur_state][column];
		}
		// normalization process
		for (int cur_state = 0; cur_state < state_num; cur_state++)
			beta[cur_state][column] = (1 / NB)*beta[cur_state][column];
	}
	/*cout << "gamma" << endl;
	for (int i = 0; i < branch_num; i++)
	{
		for (int j = 0; j < message_length / in_bit; j++)
			cout << right << setw(10) << gamma[i][j] << " ";
		cout << endl;
	}
	cout << "alpha" << endl;
	for (int i = 0; i < state_num; i++)
	{
		for (int j = 0; j < message_length / in_bit + 1; j++)
			cout << right << setw(10) << alpha[i][j] << " ";
		cout << endl;
	}
	cout << "beta" << endl;
	for (int i = 0; i < state_num; i++)
	{
		for (int j = 0; j < message_length / in_bit+1; j++)
			cout <<right<<setw(10)<<beta[i][j] << " ";
		cout << endl;
	}*/
	// determine the a posteriori prob. of each info. bit
	double *prob = new double[input_num];
	
	for (int time = 0; time < message_length / in_bit; time++)
	{
		for (int i = 0; i < input_num; i++)
		prob[i] = 0;
		double NP = 0;
		for (int msg = 0; msg < input_num; msg++)
		{
			for (int id = msg; id < branch_num; id+=input_num)
			{
				int current_state = state_table[id][in_bit];
				int next_state = state_table[id][in_bit + 1];
				prob[msg] += alpha[current_state][time] * gamma[id][time] * beta[next_state][time + 1];
			}
			NP += prob[msg];
		}
		for (int msg = 0; msg < input_num; msg++)
		{
			prob[msg] = (1 / NP)*prob[msg];
			/*cout << "msg=" << msg << " " << "prob: " << prob[msg] << endl;*/
		}
		// search the maximum prob.
		double max = prob[0]; int index = 0;
		for (int msg = 1; msg < input_num; msg++)
		{
			if (prob[msg] > max) {
				index = msg; max = prob[msg];
			}
		}
		for (int i = 0; i < in_bit; i++)
			de_message[time*in_bit + i] = state_table[index][i];
	}
	delete[] prob;
	for (int branch = 0; branch < branch_num; branch++)
		delete[] gamma[branch];
	delete[] gamma;
}