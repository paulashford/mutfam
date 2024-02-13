# mutfam_permut.R
# Ash Nov 2015

# Functions for assessing mutation enrichment in a FunFam
# by permutation of mutations

# ------------------------------------------------------------------------------
# mf.mut_permut: return p-value (of no_trials trials) for permutations of:
#   no_muts over sequence seq_len given:
#   ff_observed_count muts in FunFam with amino span [ff_aa_start, ff_aa_end]

# Allows multiple muts per amino.
# mut_permutate (mutate>permutations>permut - haha...groan)
# ------------------------------------------------------------------------------
# For n trials (typically 1000), distribute the no_muts randomly over sequence
# seq_len allowing multiple muts per amino.
# Count number of trials where observed mut count for FunFam boundary exceeds permuted version
mfp.mut_permut <- function(no_trials, no_muts, seq_len, ff_observed_count, ff_aa_start, ff_aa_end){

  mut_p_val <- 0;

  # Unmutated sequence
  #seq <- rep(0, seq_len);
  for (i in 1:no_trials){
    # Unmutated sequence
    seq <- rep(0, seq_len);

    # Get permutation of mutated aminos in this trial
    mutated_aminos <- mfp.mut_distrib(no_muts,seq_len);

    # Get a frequency table (to account for multiple muts at any given position)
    mutated_aminos_counts <- count(mutated_aminos);

    # Update sequence's mutations
    # **** Not correct - doesn't sum muts at same amino so under counts random muts!
    ##seq[mutated_aminos] <- seq[mutated_aminos] + 1;

    # Simple function to return muts at given pos, or 0 if non (no NAs!)
    freqMut <- function(aNum){
      freq <- mutated_aminos_counts[mutated_aminos_counts$x==aNum,"freq"];
      if (length(freq)==0){freq <- 0};
      return(freq);
      #print(paste("Amino num:", aNum, " Freq:",freq, sep=""))
    }
    # Now update all aminos at once... for given mutated_aminos indices
    seq[mutated_aminos] <- sapply(mutated_aminos, function(aNum) seq[aNum] <- seq[aNum] + freqMut(aNum))


    # Freq table of muts (maybe watch for non-integer treatment if large vals?)
    dt_mutated_aminos_list <- data.table(table(mutated_aminos));


    # Number of muts within FunFam boundary in this trial
    r_ffmuts <- sum(seq[ff_aa_start:ff_aa_end]);
    seqstring <- paste(seq, collapse='');
    #print (paste("Trial: ",i ," seq: ",seqstring, " sum seq: ",sum(seq), "muts in ff: ",r_ffmuts,collapse=''));

    # If random greater or eq to observed, increase p-val (i.e. worsen it...)
      if (r_ffmuts >= ff_observed_count){
        mut_p_val <- mut_p_val + 1;
      }
    }

  # Normalise p-val by no of random trials
  mut_p_val <- mut_p_val / no_trials;
  return(mut_p_val);
}

# Distribute no_muts over sequence of length seq_len
# allow.dupes=T: Multiple muts at given position are possible
mfp.mut_distrib <- function(no_muts, seq_len, allow.dupes=T){
max_trials = 100000;
    if ((no_muts > seq_len) && (allow.dupes == F)){
      print("Number of mutations should not exceed sequence length when allow.dupes == F");
      return(0);
    }

    if (allow.dupes == T){
        return (sort(floor(runif(no_muts, min=1, max=seq_len+1))));
    } else {
      # Bit crude - if *don't* want any dupes, just keep trying until distinct mut amino nos returned
      trial_num <- 0;
      # With certain no_muts / seq_len may not find solution so don't get stuck
      while (trial_num < max_trials){
        trial_num <- trial_num + 1;
        trial_distrib <- floor(runif(no_muts, min=1, max=seq_len + 1));
        # If the union of set of aminos from trial_distrib with itself is equal
        # to number of elements, then they are all unique...
        if (length(union(trial_distrib, trial_distrib)) == length(trial_distrib)){
          #print(trial_num)
          return(sort(trial_distrib));
        }
      }
      print(paste("No unique amino acid pemutations found before max. trial nums reached (" , max_trials , ")"), sep="");
    }
}

# For sequence positions from 1:seq_len...
# If we know which positions are 'correct' (correct_amino_vec) and
# we have some predicted seq positions (pred_amino_vec) then
# p-value is fraction of permutation trials where intersection of correct positions
# and trial positions is *at least* as good as intersection of predictions and correct positons...
mfp.seq_sel_pval <- function(seq_len, correct_amino_vec, pred_amino_vec, no_trials = 10000){

  # How many correct preds?
  correct_preds <- length(intersect(correct_amino_vec, pred_amino_vec));

  if (correct_preds == 0){
    print("No predicted sequence positions match the correct ones!");
    return(0);
  }

  # Number of predicted positions
  num_preds <- length(pred_amino_vec);

  mut_p_val <- 0;

  for (i in 1:no_trials){
    # Unmutated sequence
    seq <- rep(0, seq_len);

    # Get permutation of mutated aminos in this trial
    seq_permut <- mfp.mut_distrib(num_preds, seq_len, allow.dupes=F);

    # If we get at least the same number of correct positions in this random
    # trial then inc (worsen) p-val
    if (length(intersect(seq_permut, correct_amino_vec)) >= correct_preds){
      mut_p_val <- mut_p_val +1;
    }

  }

  # Normalise p-val by no of random trials
  mut_p_val <- mut_p_val / no_trials;
  return(mut_p_val);
}
