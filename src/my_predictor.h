// my_predictor.h
// This file contains a sample my_predictor class.
// It is a simple 32,768-entry gshare with a history length of 15.
// Note that this predictor doesn't use the whole 32 kilobytes available
// for the CBP-2 contest; it is just an example.

#include "tage_predictor.h"
#include "snp_predictor.h"
#include "my_predictor_gshare.h"

class my_update : public branch_update {
public:
	unsigned int index;
};

class my_predictor : public branch_predictor {
public:
#define HISTORY_LENGTH	15
#define TABLE_BITS	15
#define HYSTERESIS	2
#define NUM_PCS_TRACKED	100000000
#define NUM_MISPRED_EPOCHS	20
#define MISPRED_THRES	100000
	branch_info bi;
	my_update u;
	branch_update *tage_update;
	branch_update *snp_update;
	branch_update *gshare_update;
	tage_predictor tage_pred;
	snp_predictor snp_pred;
	gshare_predictor gshare_pred;
	long long int tage_mispred, snp_mispred;
	long long int tot_tage_mispred, tot_snp_mispred;
	long long int tot_diff;
	long long int inst_count;
	long long int check;
	int hyst_counter;


	int mis_per_addr[NUM_PCS_TRACKED] = {0};
	bool is_difficult_pc[NUM_PCS_TRACKED] = {0};

	branch_update *predict (branch_info & b) {
		bi = b;
		tage_update = tage_pred.predict(b);
		if(is_difficult_pc[bi.address%NUM_PCS_TRACKED])
			snp_update = snp_pred.predict(b);
		// gshare_update = gshare_pred.predict(b);

		// printf("dir: %d", tage_update->direction_prediction());
		// printf("tar: %d", tage_update->target_prediction());

		// u.direction_prediction(snp_update->direction_prediction());
		// u.target_prediction(snp_update->target_prediction());
		// u.direction_prediction(tage_update->direction_prediction());
		// u.target_prediction(tage_update->target_prediction());
		// u.direction_prediction(gshare_update->direction_prediction());
		// u.target_prediction(gshare_update->target_prediction());

	if(is_difficult_pc[bi.address%NUM_PCS_TRACKED] == 0)
		return tage_update;
	else
		return snp_update;
	}

	void update (branch_update *u, bool taken, unsigned int target) {

		// tage_mispred = (tage_update->direction_prediction () != taken) + (tage_update->target_prediction () != target);
		// snp_mispred = (snp_update->direction_prediction () != taken) + (snp_update->target_prediction () != target);
		tage_mispred = (tage_update->direction_prediction () != taken);
		if(is_difficult_pc[bi.address%NUM_PCS_TRACKED] == 1)
			snp_mispred = (snp_update->direction_prediction () != taken);

		tot_snp_mispred += snp_mispred;
		tot_tage_mispred += tage_mispred;

		tot_diff += tage_mispred ^ snp_mispred;

		if(tage_mispred){
			mis_per_addr[bi.address%NUM_PCS_TRACKED]++;
			if(mis_per_addr[bi.address%NUM_PCS_TRACKED] > MISPRED_THRES)
				is_difficult_pc[bi.address%NUM_PCS_TRACKED] = 1;
			hyst_counter--;
			if(hyst_counter < -1*(HYSTERESIS-1))
				hyst_counter = -1*(HYSTERESIS-1);
		} else if (snp_mispred)
		{
			hyst_counter++;
			if(hyst_counter > HYSTERESIS)
				hyst_counter = HYSTERESIS;
		}
		if(hyst_counter != 0)
			check++;

		tage_pred.update(u,taken,target);
		if(is_difficult_pc[bi.address%NUM_PCS_TRACKED] == 1)
			snp_pred.update(u,taken,target);

		inst_count++;
		// if(inst_count == 1e7){
		// 	double avg_mispred = 0.0;
		// 	long long int mispred_arr[NUM_MISPRED_EPOCHS] = {0};
		// 	for(int i=0; i<NUM_PCS_TRACKED; i++){
		// 		avg_mispred += mis_per_addr[i];
		// 		if((mis_per_addr[i]/10) < NUM_MISPRED_EPOCHS)
		// 			mispred_arr[(mis_per_addr[i]/10)]++;
		// 		else
		// 			mispred_arr[NUM_MISPRED_EPOCHS-1]++;
		// 		}
		// 	for(int i=0; i<NUM_MISPRED_EPOCHS; i++){
		// 		printf ("%0.3f epoch - %d\n", 1.0 * (mispred_arr[i]), i);
		// 	}
		// 	avg_mispred /= NUM_PCS_TRACKED;
			// printf ("%0.3f avg_mispred\n", 1.0 * (avg_mispred));
			// printf ("%0.3f check\n", 1000.0 * (check / 1e8));
			// printf ("%0.3f tot_diff\n", 1000.0 * (tot_diff / 1e8));
			// }
			
	}
};
