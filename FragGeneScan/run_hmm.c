#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include "hmm.h"


int main (int argc, char **argv){

  int i, j, c, num_seq;
  HMM hmm;
  char *obs_seq, *obs_head;
  TRAIN train;
  int wholegenome;
  int format=0;
  FILE *fp_out, *fp_aa, *fp_dna, *fp;
  char hmm_file[4096] = "";
  char aa_file[4096] = "";
  char seq_file[4096] = "";
  char out_file[4096] = "";
  char dna_file[4096] = ""; 
  char train_file[4096] = "";
  char mstate_file[4096] = "";
  char rstate_file[4096] = "";
  char nstate_file[4096] = "";
  char sstate_file[4096] = "";
  char pstate_file[4096] = "";
  char s1state_file[4096] = "";     /* stop codon of gene in - stand */
  char p1state_file[4096] = "";
  char dstate_file[4096] = "";
  char train_dir[4096] = "";
  int count=0;
  char mystring[1000] = "";
  int *obs_seq_len;
  int bp_count;  /* count the length of each line in input file */

  strncpy(train_dir, argv[0], strlen(argv[0])-12);
  strcat(train_dir, "train/");
  strcpy(mstate_file, train_dir);
  strcat(mstate_file, "gene");
  strcpy(rstate_file, train_dir);
  strcat(rstate_file, "rgene");
  strcpy(nstate_file, train_dir);
  strcat(nstate_file, "noncoding");
  strcpy(sstate_file, train_dir);
  strcat(sstate_file, "start");
  strcpy(pstate_file, train_dir);
  strcat(pstate_file, "stop");
  strcpy(s1state_file, train_dir);
  strcat(s1state_file, "stop1");
  strcpy(p1state_file, train_dir);
  strcat(p1state_file, "start1");
  strcpy(dstate_file, train_dir);
  strcat(dstate_file, "pwm");


  /* read command line argument */
  if (argc <= 8){    
    fprintf(stderr, "ERROR: You missed some parameters for input\n");
    print_usage();
    exit(EXIT_FAILURE);
  }

  while ((c=getopt(argc, argv, "fs:o:w:t:")) != -1){

    switch (c){
    case 's':
      strcpy(seq_file, optarg);
      if (access(seq_file, F_OK)==-1){
	fprintf(stderr, "ERROR: Sequence file [%s] does not exist\n", seq_file);
	print_usage();
	exit(EXIT_FAILURE);
      }
      break;  
    case 'w':
      wholegenome = atoi(optarg);
      if (wholegenome != 0 && wholegenome != 1){
	fprintf(stderr, "ERROR: An incorrect value for the option -w was entered\n");
	print_usage();
	exit(EXIT_FAILURE);
      }
      break;
    case 'o':
      strcpy(out_file, optarg);
      break;
    case 't':
      strcpy(train_file, optarg);
      strcpy(hmm_file, train_dir);
      strcat(hmm_file, train_file);

      if (access(hmm_file, F_OK)==-1){
	fprintf(stderr, "ERROR: The file for model parameters [%s] does not exist\n", hmm_file);
	print_usage();
	exit(EXIT_FAILURE);
      }
      break;
    case 'f':
      format = 1;
      break;
    }
  }

  
  /* check whether the specified files exist */
  if (access(mstate_file, F_OK)==-1){
    fprintf(stderr, "Forward prob. file [%s] does not exist\n", mstate_file);
    exit(1);
  }
  if (access(rstate_file, F_OK)==-1){
    fprintf(stderr, "Backward prob. file [%s] does not exist\n", rstate_file);
    exit(1);
  }
  if (access(nstate_file, F_OK)==-1){
    fprintf(stderr, "noncoding prob. file [%s] does not exist\n", nstate_file);
    exit(1);
  }
  if (access(sstate_file, F_OK)==-1){
    fprintf(stderr, "start prob. file [%s] does not exist\n", sstate_file);
    exit(1);
  }
  if (access(pstate_file, F_OK)==-1){
    fprintf(stderr, "stop prob. file [%s] does not exist\n", pstate_file);
    exit(1);
  }
  if (access(s1state_file, F_OK)==-1){
    fprintf(stderr, "start1 prob. file [%s] does not exist\n", s1state_file);
    exit(1);
  }
  if (access(p1state_file, F_OK)==-1){
    fprintf(stderr, "stop1 prob. file [%s] does not exist\n", p1state_file);
    exit(1);
  }
  if (access(dstate_file, F_OK)==-1){
    fprintf(stderr, "pwm dist. file [%s] does not exist\n", dstate_file);
    exit(1);
  }
  if (access(hmm_file, F_OK)==-1){
    fprintf(stderr, "hmm file [%s] does not exist\n", hmm_file);
    exit(1);
  }

  /* read all initial model */
  hmm.N=NUM_STATE;
  get_train_from_file(hmm_file, &hmm, mstate_file, rstate_file, nstate_file, sstate_file, pstate_file,s1state_file, p1state_file, dstate_file, &train);
 
  /* create output file name */
  strcpy(aa_file, out_file);
  strcat(aa_file, ".faa");
  strcpy(dna_file, out_file);
  strcat(dna_file, ".ffn");

  remove (out_file);
  remove (aa_file);
  remove (dna_file);

  fp_aa = fopen (aa_file , "w");
  fp_out = fopen (out_file , "w");
  fp_dna = fopen (dna_file , "w");
  fp = fopen (seq_file, "r");

  while ( fgets (mystring , sizeof mystring , fp) ){
    if (mystring[0] == '>'){
      count++;
    }
  }
  num_seq = count;
  obs_seq_len = (int *)malloc(num_seq * sizeof(int));
  printf("no. of seqs: %d\n", num_seq);


  i = 0;
  count = 0;
  rewind(fp);
  while ( fgets (mystring , sizeof mystring , fp) ){
    if (mystring[0] == '>'){
      if (i>0){
        obs_seq_len[count] = i;
        count++;
      }
      i = 0;
    }else{
      bp_count = strlen(mystring)-1;
      while(mystring[bp_count-1] == 10 || mystring[bp_count-1]==13){
	bp_count --;
      }

      i += bp_count;
    }
  }
  obs_seq_len[count] = i;

  count=-1;
  rewind(fp);
  j ==0;

  while ( fgets (mystring , sizeof mystring  , fp) ){

    if (mystring[0] == '>'){

      if (count>=0 && j>0){

 	get_prob_from_cg(&hmm, &train, obs_seq);

 	if (strlen(obs_seq)>70){
	  viterbi(&hmm, obs_seq, fp_out, fp_aa, fp_dna, obs_head, wholegenome, format);
	}
      }

      bp_count = strlen(mystring)-1;
      while(mystring[bp_count-1] == 10 || mystring[bp_count-1]==13){
	bp_count --;
      }

      obs_head = (char *)malloc((bp_count+1) * sizeof(char));
      memset(obs_head, 0, (bp_count+1) * sizeof(char));
      memcpy(obs_head, mystring, bp_count);

      if (count== -1 || (count>=0 && j>0)){
	count++;
	obs_seq = (char *)malloc(obs_seq_len[count] * sizeof(char) + 1);
	memset(obs_seq, 0, obs_seq_len[count] * sizeof(char) + 1);
      }
      j = 0;

    }else{
      bp_count = strlen(mystring)-1;
      while(mystring[bp_count-1] == 10 || mystring[bp_count-1]==13){
	bp_count --;
      }
      memcpy(obs_seq+j, mystring, bp_count);
      j += bp_count;
    }
  }

  if (count>=0){

    get_prob_from_cg(&hmm, &train, obs_seq);

    if (strlen(obs_seq)>70){
      viterbi(&hmm, obs_seq, fp_out, fp_aa, fp_dna, obs_head, wholegenome, format);
    }
  }
  
  free(obs_seq_len);
  free(obs_head);
  fclose(fp_out);
  fclose(fp_aa);
  fclose(fp_dna);
  fclose(fp);
}


