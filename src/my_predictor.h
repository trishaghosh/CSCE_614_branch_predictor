/* 
This file contains code for simulating the SNP predictor as described
in Renee St. Amant, Daniel A. Jimenez and Doug Burger, 
Low-Power, High-Performance Analog Neural Branch Prediction, 
Proceedings of the 41st Annual International Symposium on Microarchitecture 
(MICRO-41), Lake Como, Italy, November 2008.  

To use this code, download the CBP2 infrastructure 
(see http://cava.cs.utsa.edu/camino/cbp2/cbp2-infrastructure-v2/)
and use this code as "my_predictor.h" in the src directory.

*/

#include <math.h>

typedef char byte;

class snp_update : public branch_update {
public:
	double	yout;
	unsigned int *addresses;
	unsigned int set;
	byte *global_history;
	unsigned int address;
	snp_update (void) { }
};

class snp : public branch_predictor {
private:
	int	update_index, 		// next available update state entry
		num_entries,		// number of blocks in a column
		min_entries,		// minimum number of blocks in a column
		history_length,		// global history length
		path_length,		// number of addresses in path history
		theta,			// training threshold
		tc,			// counter for threshold adaptivity
		max_weight[1000], 	// i'th weight maximum value (run-time constants)
		min_weight[1000],	// i'th weight minimum value (run-time constants)
		weight_bits,		// number of bits per weight
		abits,			// number of address bits in hash
		nsets,			// number of sets
		num_bias_entries,	// number of bias weights
		history_modulus,	// controls redundant history
		modulo_type,		// type of redundant (or not) history
		cut,			// beyond this weight positions reduce weights by one bit
		block_size,		// size of a block
		virtual_block_size;	// size of a block for hashing purposes

	double 
		*coefficients, 		// coefficients from constructor
		bias_coefficient;	// bias coefficient

	bool ahead_pipelined;		// if true then simulate ahead pipelining

	int
		*lgnents,		// log2 (number of entries) in k'th column
		*nents;			// number of entries in k'th column

	snp_update
		update_buf[100];	// queue of update state entries

	byte
		***weights;		// array of weights vectors

	byte
		*bias_weights;		// bias weight array

	byte
		*global_history,	// global pattern histories
		*spec_global_history;	// speculative version

	unsigned int
		*addresses, 		// the sequence of recent branch PCs
		*spec_addresses;	// the sequence of recent branch PCs

	// return an array of weights vectors of the right sizes

	byte **make_weights (void) {
		byte **v = new byte*[num_entries];
		for (int i=0; i<num_entries; i++) {
			int n = 0;
			for (int j=0; j<history_length; j+=block_size) {
				int k = j / block_size;
				if (nents[k] > i) n += block_size;
			}
			v[i] = (byte *) malloc (n);
			memset (v[i], 0, n);
		}
		return v;
	}

	unsigned int bias_hash (unsigned int address) {
		return (address >> 1) % num_bias_entries;
	}

	// log base 2

	int lg (int x) {
		int c = 0;
		while (x) {
			x /= 2;
			if (x) c++;
		}
		return c;
	}

public:

	// constructor

	snp (	
		bool ahead = false,	// simulate ahead-pipelining
		int ct = 64,		// weights after this length are reduced by 1 bit
		int bs = 8,		// block size in weights
		int vbs = 2,		// "virtual" block size
		int be = 4096,		// number of bias entries
		int ne = 512, 		// number of weights vectors
		int me = 256, 		// min entries
		int hl = 128,	 	// global history length
		int ml = 8,		// history modulus
		int ns = 1,		// number of sets
		int ab = 9,		// # of bits of path address to use
		int b = 7, 		// number of bits per weight
		int mt = 3,		// "modulo type" - choose from 1 of 3 ways of choosing bits
		double *v = NULL)	// coefficients vector
		: branch_predictor () { 		
		int	i, k;

		history_length = hl;
		block_size = bs;
		ahead_pipelined = ahead;
		cut = ct;
		virtual_block_size = vbs;
		modulo_type = mt;
		history_modulus = ml;
		num_bias_entries = be;
		abits = ab;
		weight_bits = b;
		int maxw = (1 << (b - 1)) - 1;
		int minw = - maxw - 1;
		for (i=0; i<=history_length; i++) {
			if (i < cut) {
				max_weight[i] = maxw; 
				min_weight[i] = minw; 
			} else {
				max_weight[i] = maxw/2;
				min_weight[i] = minw/2;
			}
		}
		// initialize from params
		num_entries = ne;
		min_entries = me;
		coefficients = (double *) malloc (history_length * sizeof (double));
		if (v) {
			bias_coefficient = v[0];
			for (i=0; i<history_length; i++) {
				if (v[i] >= 1) coefficients[i] = v[i+1];
				else coefficients[i] = 1.0;
			}
		} else {
			bias_coefficient = 1.0;
			for (i=0; i<history_length; i++) coefficients[i] = 1.0;
		}
		lgnents = new int[history_length / block_size + 1];
		nents = new int[history_length / block_size + 1];
		int z = num_entries;
		for (k=0,i=0; i<history_length; i+=block_size, k++) {
			nents[k] = z;
			lgnents[k] = lg (z);
			z /= 2;
			if (z < min_entries) z = min_entries;
		}
		theta = (int) (2.14 * hl + 20.58);
		tc = 0;
		nsets = ns;

		// the history length must be a multiple of block_size

		while (history_length % block_size) history_length++;

		// a path element accounts for 8 history elements

		path_length = history_length / virtual_block_size + 100;
		addresses = new unsigned int[path_length];
		spec_addresses = new unsigned int[path_length];

		// initialize perceptron weights

		bias_weights = new byte[num_bias_entries];

		weights = new byte **[nsets];
		for (k=0; k<nsets; k++)
			weights[k] = make_weights ();

		// initialize histories

		global_history = (byte *) malloc (history_length);
		memset (global_history, 0, history_length);
		spec_global_history = (byte *) malloc (history_length);
		memset (spec_global_history, 0, history_length);
		for (i=0; i<100; i++) {
			update_buf[i].addresses = new unsigned int [path_length];
			update_buf[i].global_history = (byte *) malloc (history_length);
			memset (update_buf[i].global_history, 0, history_length);
		}
		update_index = 0;
	}

	void shift_history (byte *v, bool t) {
		memmove (&v[1], &v[0], history_length-1);
		v[0] = t;
	}

	// shift an address into the sequence of recent branch PCs

	void shift_address (unsigned int *v, unsigned int addr) {
		memmove (&v[1], &v[0], path_length * sizeof (unsigned int));
		v[0] = addr;
	}

	unsigned int shuffle (unsigned int *v, int n, int k) {
		int	i, j, x, mask, count;
	
		mask = 2;
		count = 0;
		i = k % n;
		x = 0;
		for (j=0; j<lgnents[k]; j++) {
			x <<= 1;
			x |= !!(mask & v[i]);
			count++;
			if (count == n) {
				mask <<= 1;
				count = 0;
			}
			i++;
			if (i == n) i = 0;
		}
		return x;
	}

	unsigned int index_entries (unsigned int x, unsigned int *v, int k) {
		int n = block_size / virtual_block_size;
		assert (n);
		unsigned int z = shuffle (v, n, k) % (1<<abits);
		if (!ahead_pipelined)
			z ^= x;
		return z % nents[k];
	}

	byte *hisb (snp_update *u, int i) {
		switch (modulo_type) {
			case 3:
			if ((i/history_modulus) & 1)
				return &u->global_history[(i % block_size) + (i / history_modulus) * virtual_block_size];
			else
				return &u->global_history[i % history_modulus];
			default:
				return &u->global_history[i];
		}
	}

	double compute_output (snp_update *u) {
		int	i, k;
		double sum_pos, sum;
		unsigned int *v = u->addresses;
		sum_pos = 0;
		sum = bias_coefficient * bias_weights[bias_hash (u->address)];
		sum_pos = 0;
		byte **s = weights[u->set];
		k = 0;
		int vk = 0;
		for (i=0; i<history_length; i+=block_size) {
			byte *bl = &s[index_entries(u->address, &v[vk], k)][i];
			byte *hp = hisb (u, i);
			double *co = &coefficients[i];
			for (int j=0; j<block_size; j++) {
				int h = hp[j] ? 1 : -1;
				sum_pos += h * bl[j] * co[j];
			} 
			k++;
			vk += block_size / virtual_block_size;
		}
		int y = (int) ((sum + sum_pos));
		return y;
	}

	// train the perceptron formed by looking at one weight from each
	// of the addresses in v[]

	void train (snp_update *u, bool taken, bool correct) {
		int	i, k, a;
		byte *c;
		unsigned int *v = u->addresses;

		a = fabs(u->yout);

		// dynamic threshold updating per Seznec

		if (!correct) {
			tc++;
			if (tc >= 1) {
				theta++;
				tc = 0;
			}
		}
		if (correct && a < theta) {
			tc--;
			if (tc <= -1) {
				theta--;
				tc = 0;
			}
		}
		if (correct && a >= theta) return;
		c = &bias_weights[bias_hash (u->address)];
		if (taken) {
			if (*c < max_weight[0]) (*c)++; 
		} else {
			if (*c > min_weight[0]) (*c)--; 
		}
		k = 0;
		int vk = 0;
		for (i=0; i<history_length; i+=block_size) {
			for (int j=0; j<block_size; j++) {
				bool h = (bool) *(hisb (u, i+j));
				byte *c = &weights[u->set][index_entries(u->address, &v[vk], k)][i+j];
				if (taken == h) {
					if (*c < max_weight[i+1]) (*c)++;
				} else {
					if (*c > min_weight[i+1]) (*c)--;
				}
			}
			k++;
			vk += block_size / virtual_block_size;
		}
	}

	// predict the branch at this address 

	branch_update *predict (unsigned int address) { 
		snp_update *u; 
		u = &update_buf[update_index++];
		if (update_index >= 100) update_index = 0;
		u->address = address;
		u->set = address % nsets;
		memcpy (u->global_history, spec_global_history, history_length);
		memcpy (u->addresses, spec_addresses,
			sizeof (unsigned int) * path_length);
		u->yout = compute_output (u);
		u->direction_prediction (u->yout >= 0);
		shift_address (spec_addresses, u->address);
		shift_history (spec_global_history, u->direction_prediction() ^ !!(u->address & 4));
		return u;
	}

	branch_update *predict (branch_info & b) {
		return predict (b.address);
	}

	// update an executed branch

	void update (branch_update *u, bool taken) {
		snp_update *hu = (snp_update *) u;
		
		bool correct = taken == u->direction_prediction();

		train (hu, taken, correct);
		shift_history (global_history, taken ^ !!(hu->address & 4));
                shift_address (addresses, hu->address);

		// recover from misprediction

		if (taken != hu->direction_prediction ()) {
			memcpy (spec_global_history, global_history, history_length);
			memcpy (spec_addresses, addresses,
				path_length * sizeof (unsigned int));
		}
	}
};


class my_update : public branch_update {
public:
	branch_update *p;
};

class my_predictor : public branch_predictor {
public:
	snp *pred;
	branch_info bi;
	my_update u;

	my_predictor (void) { 
		double a = 0.04;
		double b = 0.05;
		double *sp = new double[1024];
		for (int i=0; i<1024; i++) 
			sp[i] = 1.0 / (a + b * i);
		pred = new snp (
			false,  // ahead pipelining
			32,     // cut
			8,      // block size
			2,      // virtual block size
			4096,   // bias entries
			512,    // weights vectors
			256,    // min weights vectors
			128,    // global history length
			8,      // history modulus
			1,      // number of sets
			9,      // number of address bits in path component
			7,      // number of bits per weight
			3,      // modulo type
			sp);    // coefficients vector
	}

	branch_update *predict (branch_info & b) {
		bi = b;
		if (b.br_flags & BR_CONDITIONAL) {
			u.p = pred->predict (b);
			u.direction_prediction (u.p->direction_prediction());
		} else {
			u.direction_prediction (true);
		}
		u.target_prediction (0);
		return &u;
	}

	void update (branch_update *, bool taken, unsigned int target) {
		if (bi.br_flags & BR_CONDITIONAL) {
			pred->update (u.p, taken);
		}
	}
};
