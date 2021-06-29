#include "seed.h"

std::ostream& operator<<(std::ostream& os, seed& s){
  os << s.str();
  return os;
}


seed::seed(const int nbcyclepos, const int maxcyclepos){
  _seedNbCyclePos  = nbcyclepos;
  _seedMaxCyclePos = maxcyclepos;

  if (_seedNbCyclePos > _seedMaxCyclePos) {
    _ERROR("Seed Cycle Pos","cannot set " << _seedNbCyclePos << " positions in a cycle of size " << _seedMaxCyclePos);
  }

  _span          = gv_minspan;
  _seedMotif_int = new int[_span];

  for (int i = 0; i < _span; i++)
    _seedMotif_int[i] = 0;

  if (gv_signature_flag) {
    checksignature();
    setsignature();
  }

  if (_seedNbCyclePos > 0) {
    _seedCyclePos_int = new int[_seedNbCyclePos];
    for (int i = 0; i < _seedNbCyclePos; i++)
      _seedCyclePos_int[i] = i;
  } else {
    _seedCyclePos_int = NULL;
  }

  set_hmove();

  if (gv_nbruns == 0){
    next();
  } else {    // FIXME for gv_nbruns == 0 ? next is necessary ??
    random();
  }
}


/* Constructor
 */
seed::seed(string str_desc, const bool update_gv_span_weight) {
  // seed span
  {
    int pos = 0;
    _span   = 0;
    while(pos < (int)str_desc.length()) {
      convert(str_desc, &pos);
      _span++;
    }
    if (_span == 0) {
      _ERROR("Single Seed Parsing","empty string provided");
    }
  }

  _seedMotif_int = new int[_span];

  _seedCyclePos_int   = NULL;
  _seedNbCyclePos     = 0;
  _seedMaxCyclePos    = 1;

  /* seed shape */
  {
    int pos = 0;
    for (int i = 0; i < _span; i++)
      _seedMotif_int[i] = convert(str_desc, &pos);
  }

  set_hmove();

  /* reset global parameters to be consistent with the input seed */
  if (update_gv_span_weight) {
    if (gv_minspan > _span) {gv_minspan = _span; _WARNING("Single Seed Parsing","ajusting minspan to "<< _span);}
    if (gv_maxspan < _span) {gv_maxspan = _span; _WARNING("Single Seed Parsing","ajusting maxspan to "<< _span);}
    if (gv_weight_interval_flag) {
      double w =  weight();
      if (gv_minweight > w)  {gv_minweight = w; _WARNING("Single Seed Parsing","ajusting minweight to "<< w);}
      if (gv_maxweight < w)  {gv_maxweight = w; _WARNING("Single Seed Parsing","ajusting maxweight to "<< w);}
    }
  }
}


/* Destructor
 */
seed::~seed(){
  if  (_seedMotif_int){
    delete[] _seedMotif_int;
    _seedMotif_int = NULL;
  }
  if (_seedCyclePos_int) {
    delete[] _seedCyclePos_int;
    _seedCyclePos_int = NULL;
  }
}




/* Symetric check
 */
bool seed::symetric()  {
  for (int i = 0; i < _span/2; i++)
    if (_seedMotif_int[i] != _seedMotif_int[_span-1-i])
      return false;
  return true;
}


/* Equal
 */
bool seed::equal(seed * other) {
  if (_span != other->span())
    return false;

  for (int i = 0; i < _span; i++)
    if (_seedMotif_int[i] != other->_seedMotif_int[i])
      return false;
  return true;
}

/* ToString
 */
string seed::str() {
  ostringstream outs;

  for (int i = 0; i < _span; i++) {
    int c = this->_seedMotif_int[i];
    outs << rconvert(c);
  }
  if ( _seedCyclePos_int ) {
    outs << ":";

    /* convert to "seed begin position" */
    vector<int> cp = vector<int>(0);
    for (int p = 0; p < _seedNbCyclePos; p++)
      cp.push_back(((_seedCyclePos_int[p] - (_span%_seedMaxCyclePos) + _seedMaxCyclePos)%_seedMaxCyclePos) + 1);
    sort(cp.begin(),cp.end());

    for (int p = 0; p < _seedNbCyclePos; p++) {
      if (p > 0)
        outs << ",";
      outs << (cp[p]);
    }
    outs << "/" << (_seedMaxCyclePos);
  }
  return outs.str();
}


/* Weight of a given seed
 */
double seed::weight() {
  double weight = 0.0;
  for (int i = 0; i < _span; i++)
    weight += gv_bsel_weight[_seedMotif_int[i]];
  return weight;
}

/* Selectivity according to the weight (bernoulli model estimation)
 */
double seed::selectivityFromWeight() {
  return exp(weight()*log(gv_bsel_minprob));
}


/* Check if the current seed is acceptable
 */
bool seed::acceptable() {
  if (!gv_vectorized_flag &&
      (
       gv_bsel_weight[_seedMotif_int[0]]       <= (1e-32) ||
       gv_bsel_weight[_seedMotif_int[_span-1]] <= (1e-32)
       )
      ) {
    return false;
  }

  if (gv_symetric)
    if (!symetric())
      return false;

  double w = weight();

  if (gv_weight_interval_flag && (w < gv_minweight || w > gv_maxweight))
    return false;
  return true;
}


/* Check if the signature is realisable according to the span : stop the program otherwise
 */
void seed::checksignature() {
  int sumsignature = 0;
  for (int i = 0; i < gv_seed_alphabet_size; i++)
    sumsignature += gv_signature[i];
  if (sumsignature > gv_minspan ){
    _ERROR("Invalid Signature","size of the signature is " << sumsignature << " and does not fit inside minimal span " << gv_minspan);
  }
}

/* Set the very first seed signature
 */
void seed::setsignature() {
  int pos = 0;
  for (int b = gv_seed_alphabet_size-1; b >= 0 && pos < _span; b--){
    int lettersignature =  gv_signature[b];
    while(lettersignature > 0) { // fill with b
      _seedMotif_int[pos++] = b;
      lettersignature--;
    }
  }
  while (pos < _span) { // fill remaining part with jokers
    _seedMotif_int[pos++] = 0;
  }
}



/* Set all the elements of position >= first_pos that are <= b in the seed in the decreasing order
 * this function is used to enumerate signatures and inner used.
 * - first_pos gives the first position to be sorted
 * - b gives the maximal value that can be sorted (all others remain unchanged)
 */
void seed::reorder(int first_pos,int b) {
  int last_bm = -1;
  for (int i = first_pos; i < _span; i++) {
    if (_seedMotif_int[i] == b && last_bm >= 0) {
      //swap elements
      int e = _seedMotif_int[i];
      _seedMotif_int[i] = _seedMotif_int[last_bm];
      _seedMotif_int[last_bm] = e;
      i = last_bm - 1;
      last_bm = -1;
    } else {
      if (_seedMotif_int[i] < b && last_bm == -1)
        last_bm = i;
    }
  }
}

/* Generate the next seed
 */
int seed::next() {

  /* A) positions move first */
  if ( _seedCyclePos_int ) {
    int i = _seedNbCyclePos - 1;
    int j = _seedMaxCyclePos - 1;
    while (i>=0 && _seedCyclePos_int[i] == j) {
      i--;
      j--;
    }
    if (i < 0) {
      for (int j = 0; j < _seedNbCyclePos; j++)
        _seedCyclePos_int[j] = j;
      goto nextseed;
    } else {
      _seedCyclePos_int[i]++;
      for (int j = i+1; j < _seedNbCyclePos; j++) {
        _seedCyclePos_int[j] = _seedCyclePos_int[j-1] + 1;
      }
      return 1;
    }
  }

 nextseed:
  /* B) seed move second */
  if (gv_signature_flag || gv_signature_shuffle_from_m_pattern_flag) {
  start_sign:
    /* 1) non trivial signature enumeration */
    for (int b = 1; b < gv_seed_alphabet_size; b++ ){

      int first_b_from_right_to_left = -1;
      int last_bm_before_b           = -1;

      for (int j=_span-1; j >= 0; j--) {
        if ( _seedMotif_int[j] <  b ) {
          last_bm_before_b = j;
        } else if ( _seedMotif_int[j] == b && last_bm_before_b >= 0 ) {
          first_b_from_right_to_left = j;
          break;
        }
      }
      if (first_b_from_right_to_left >= 0) {
        //swap elements
        {
          int e = b;
          _seedMotif_int[first_b_from_right_to_left] =  _seedMotif_int[last_bm_before_b];
          _seedMotif_int[last_bm_before_b] = e;
        }
        //reorder elements b and <b after this swap position
        reorder(last_bm_before_b,b);

        //reorder the elements <b (only at position where such elements are)
        for (int bl = b-1; bl >= 0; bl-- )
          reorder(0,bl);
        while(!acceptable())
          goto start_sign;
        return 1;
      }
    }
    // end of the shuffle seed
    if (gv_signature_shuffle_from_m_pattern_flag) {
      //reorder elements b and <b after this swap position
      reorder(0,gv_seed_alphabet_size-1);
      return 0;
    }
    // span changing
    if (_span < gv_maxspan) {
      _span ++;
      delete[] _seedMotif_int;
      _seedMotif_int = new int[_span];
      setsignature();
      while(!acceptable())
        goto start_sign;
      return 1;
    } else {
      return 0;
    }


  } else {

    /* 2) complete enumeration */
  start:
    int j=_span-1;

    while(_seedMotif_int[j] == (gv_seed_alphabet_size-1)){
      j--;
      if (j < 0) {
        if (_span < gv_maxspan) {
          _span ++;
          delete[] _seedMotif_int;
          _seedMotif_int = new int[_span];
          goto reset;
        } else {
          return 0;
        }
      }
    }

    _seedMotif_int[j]++;
  reset:
    j++;
    for (int i = j; i < _span; i++) {
      _seedMotif_int[i] = 0;
    }

    while(!acceptable()){
      goto start;
    }
    return 1;
  }/* 2) */
}

/* Generate a random seed
 */
int seed::random(){

  /* A) select the seed shape */
  if (gv_signature_shuffle_from_m_pattern_flag) {
  swap_shuffle:
    for (int i = 0; i < _span; i++){
      int j = rand()%_span;
      int letter_tmp    = _seedMotif_int[i];
      _seedMotif_int[i] = _seedMotif_int[j];
      _seedMotif_int[j] = letter_tmp;
    }
    while(!acceptable())
      goto swap_shuffle; /* FIXME : can loop infinitely if a bad signature is given on -m */

  } else if (gv_signature_flag) {

    /* 1) signature random generation */
    // check signature
    checksignature();

    // select a new span
    int newspan = rand()%(gv_maxspan - gv_minspan + 1) + gv_minspan;
    if (newspan != _span) {
      _span = newspan;
      delete[] _seedMotif_int;
      _seedMotif_int = new int[_span];
    }


    // "fill and swap" randomization :

    // fill
    setsignature();
    // swap
  swap:
    for (int i = 0; i < _span; i++){
      int j = rand()%_span;
      int letter_tmp    = _seedMotif_int[i];
      _seedMotif_int[i] = _seedMotif_int[j];
      _seedMotif_int[j] = letter_tmp;
    }

    while(!acceptable())
      goto swap; /* FIXME : can loop infinitely if a bad signature is given */

  } else {

  start:

    /* 2) random  enumeration */
    int newspan = rand()%(gv_maxspan - gv_minspan + 1) + gv_minspan;
    if (newspan != _span) {
      _span = newspan;
      delete[] _seedMotif_int;
      _seedMotif_int = new int[_span];
    }

  refill:
    for (int i = 0; i < _span; i++) {
      _seedMotif_int[i] = rand() % gv_seed_alphabet_size;
    }



    /* heuristic weight trick to modify the seed to reach the reasonable weight */
    if (gv_weight_interval_flag) {

      double w = weight();

      if ((w < gv_minweight) || ( w > gv_maxweight)) {
        int pos_shuffle[_span];

        /* set a shuffle of seed positions to be considered */
        for (int i = 0; i < _span; i++)
          pos_shuffle[i] = i;
        for (int i = 1; i < _span; i++) {
          int p = rand()%(i);
          pos_shuffle[i] = pos_shuffle[p];
          pos_shuffle[p] = i;
        }

        /* increase weight loop */
        while (w < gv_minweight && (_seedMotif_int[pos_shuffle[_span-1]] < gv_seed_alphabet_size - 1)) {
          for (int i = 0; i < _span; i++) {
            int p = pos_shuffle[i];
            if (_seedMotif_int[p] < gv_seed_alphabet_size - 1) {
              int      bold = _seedMotif_int[p]++;
              double  delta = gv_bsel_weight[_seedMotif_int[p]] - gv_bsel_weight[bold];
              w  += delta;
              if (w >= gv_minweight)
                goto eoup;
            }
          }
        }

      eoup:

        /* set a shuffle of seed positions to be considered */
        for (int i = 0; i < _span; i++)
          pos_shuffle[i] = i;
        for (int i = 1; i < _span; i++) {
          int p = rand()%(i);
          pos_shuffle[i] = pos_shuffle[p];
          pos_shuffle[p] = i;
        }

        /* decrease weight loop */
        while (w > gv_maxweight && _seedMotif_int[pos_shuffle[_span-1]] > 0) {
          for (int i = 0; i < _span; i++) {
            int p = pos_shuffle[i];
            if (_seedMotif_int[p] > 0) {
              int      bold = _seedMotif_int[p]--;
              double  delta = gv_bsel_weight[_seedMotif_int[p]] - gv_bsel_weight[bold];
              w  += delta;

              if (w <= gv_minweight)
                goto eodown;
            }
          }
        }
      eodown:;
      }
    }

    while(!acceptable()){
      if (!(rand()%500)){
        goto start;
      } else {
        goto refill;
      }
    }
  }


  /* B) select the seed pos */
  if ( _seedCyclePos_int ) {
    for (int i = 0; i < _seedNbCyclePos; i++) {
    re:
      _seedCyclePos_int[i] = rand()%_seedMaxCyclePos;
      for (int j = 0; j < i; j++)
        if (_seedCyclePos_int[j] == _seedCyclePos_int[i]) /* FIXME : do this more efficiently ... */
          goto re;
    }
  }
  return 1;
}


/* Edit cycle positions
 */
int seed::next_cycleposedit() {
  if ( _seedCyclePos_int ) {

    if (_edit_cycle_pos_index >= _seedNbCyclePos)
      return 0;
    if (_edit_cycle_pos_index >= 0) {
    nxinc:
      _seedCyclePos_int[_edit_cycle_pos_index]++;
      for (int i = 0; i < _seedNbCyclePos; i++)
        if ((i != _edit_cycle_pos_index && _seedCyclePos_int[_edit_cycle_pos_index] == _seedCyclePos_int[i]) || _seedCyclePos_int[_edit_cycle_pos_index] == _edit_cycle_pos_old_value)
          goto nxinc;

      if (_seedCyclePos_int[_edit_cycle_pos_index] >= _seedMaxCyclePos) {
        _seedCyclePos_int[_edit_cycle_pos_index] = _edit_cycle_pos_old_value;
        _edit_cycle_pos_index++;
        if (_edit_cycle_pos_index >= _seedNbCyclePos)
          return 0;
        _edit_cycle_pos_old_value = _seedCyclePos_int[_edit_cycle_pos_index];
        _seedCyclePos_int[_edit_cycle_pos_index] = -1;
        goto nxinc;
      }
    } else {
      _edit_cycle_pos_index = 0;
      _edit_cycle_pos_old_value = _seedCyclePos_int[_edit_cycle_pos_index];
      _seedCyclePos_int[_edit_cycle_pos_index] = -1;
      goto nxinc;
    }
    return 1;
  }

  return 0;
}


/* Generate the next swap
 */
int seed::next_swap() {
 nextswapseed:
  VERB_FILTER(VERBOSITY_HIGH, MESSAGE__("[1/3] next_swap : reverse previous swap edit of (" << ((_swap_from+_swap_and_edit_shift)%_span) << "," << ((_swap_to+_swap_and_edit_shift)%_span) << ") " << (*this)););
  // reverse previous swap (if any)
  if (_edit_operation != 1 &&  _edit_operation != 2 && _swap_from < _span && _swap_to < _span) {
    int letter_tmp                                          = _seedMotif_int[(_swap_from+_swap_and_edit_shift)%_span];
    _seedMotif_int[(_swap_from+_swap_and_edit_shift)%_span] = _seedMotif_int[(_swap_to+_swap_and_edit_shift)%_span];
    _seedMotif_int[(_swap_to+_swap_and_edit_shift)%_span]   = letter_tmp;
  }
  VERB_FILTER(VERBOSITY_HIGH, MESSAGE__("[2/3] next_swap : current seed is " << (*this)););
  // move _swap_to and _swap_from positions
  do {
    _swap_to++;
    if (_swap_to >= _span) {
      _swap_from++;
      _swap_to = _swap_from + 1;
      if (_swap_to >= _span) {
      VERB_FILTER(VERBOSITY_HIGH, MESSAGE__("[3/3] next_swap : END of swap edit"););
        return 0;
      }
    }
  } while(_seedMotif_int[(_swap_from+_swap_and_edit_shift)%_span] == _seedMotif_int[(_swap_to+_swap_and_edit_shift)%_span]);
  // do the swap
  if (_edit_operation != 1 &&  _edit_operation != 2 && _swap_from < _span && _swap_to < _span) {
    int letter_tmp                                          = _seedMotif_int[(_swap_from+_swap_and_edit_shift)%_span];
    _seedMotif_int[(_swap_from+_swap_and_edit_shift)%_span] = _seedMotif_int[(_swap_to+_swap_and_edit_shift)%_span];
    _seedMotif_int[(_swap_to+_swap_and_edit_shift)%_span]   = letter_tmp;
  }
  if (!acceptable()) {
    VERB_FILTER(VERBOSITY_HIGH, MESSAGE__("[3/3] next_swap : new seed after swap edit of (" << ((_swap_from+_swap_and_edit_shift)%_span) << "," << ((_swap_to+_swap_and_edit_shift)%_span) << ") is not accepted ... renew " << (*this)););
    goto nextswapseed;
  }
  VERB_FILTER(VERBOSITY_HIGH, MESSAGE__("[3/3] next_swap : NEW seed after swap edit of (" << ((_swap_from+_swap_and_edit_shift)%_span) << "," << ((_swap_to+_swap_and_edit_shift)%_span) << ") is " << (*this)););
  return 1;
}


int seed::next_delete() {
  if (_span > gv_minspan) {
    // generate next delete
    while (_edit_from < _span-1) {
      _edit_from++;
      if (gv_bsel_weight[_seedMotif_int[_edit_from]] <= (1e-32)) {
        // >> this "if" added to avoid next_delete produces the same seed several times
        if (_edit_from%_span > 0 && _seedMotif_int[_edit_from] == _seedMotif_int[_edit_from-1])
          continue;
        // <<
        // (1/1) set
        _edit_from_symbol = _seedMotif_int[_edit_from%_span];
        for (int j = _edit_from+1; j < _span; j++)
          _seedMotif_int[j-1] =  _seedMotif_int[j];
        _span--;
        _edit_operation = 2;
        return 1;
      }
    }
  }
  _edit_from = -1;
  _edit_from_symbol = 0;
  return 0;
}

int seed::next_insert() {
  if (_span < gv_maxspan) {
    // generate next insert
    while (_edit_from_symbol  < gv_seed_alphabet_size) {
      if (gv_bsel_weight[_edit_from_symbol] <= (1e-32)) {
        while (_edit_from < _span) {
          _edit_from++;
          // >> this "if" added to avoid next_insert produces the same seed several times
          if(_edit_from > 0 && _seedMotif_int[_edit_from-1] == _edit_from_symbol)
            continue;
          // <<
          // (1/2) realloc
          {
            int * seedMotif_int_new = new int[_span+1];
            for(int i = 0; i < _span; i++)
              seedMotif_int_new[i] = _seedMotif_int[i];
            delete[] _seedMotif_int;
            _seedMotif_int = seedMotif_int_new;
          }
          // (2/2) and set
          for (int j = _span; j > _edit_from; j--)
            _seedMotif_int[j] = _seedMotif_int[j-1];
          _seedMotif_int[_edit_from] = _edit_from_symbol;
          _span++;
          _edit_operation = 1;
          return 1;
        }
      }
      _edit_from = -1;
      _edit_from_symbol++;
    }
  }
  _edit_from = -1;
  _edit_from_symbol = 0;
  return 0;
}


/* Generate the next edit
 */
int seed::next_edit() {
  VERB_FILTER(VERBOSITY_HIGH, MESSAGE__("[x/x] next_edit :  (" << _edit_operation << ") " << (*this)););

  // if not all edits have been tried
  if (_edit_operation <= 2) {

  next_edit_retry:

    // A) reverse previous operation
    {
      int previous_edit_operation = _edit_operation;
      if (previous_edit_operation == 1 && _edit_from >= 0) {
        _edit_from--;
        next_delete();
      } else {
        if (previous_edit_operation == 2 && _edit_from >= 0) {
          _edit_from--;
          next_insert();
        }
      }
      _edit_operation = previous_edit_operation;
    }

    // B) do the next operation
    if (gv_hmove_choice % 12 >= 6) {
      // B.1) insertion before deletion
      if (_edit_operation <= 1) {
        if (next_insert()) {
          if (!acceptable())
            goto next_edit_retry;
          return 1;
        }
      }
      if (next_delete()) {
        if (!acceptable())
          goto next_edit_retry;
        return 1;
      }
    } else {
      // B.2) deletion before insertion
      if (_edit_operation != 1) {
        if (next_delete()) {
          if (!acceptable())
            goto next_edit_retry;
          return 1;
        }
      }
      if (next_insert()) {
        if (!acceptable())
          goto next_edit_retry;
        return 1;
      }
    }
    // C) everything tried ... flag this even to be faster
    _edit_operation = 3;
  }
  return 0;
}




/* Generate the next hmove
 */
int seed::next_hmove() {
  if (gv_hmove_choice%3 == 0) {
    // A) edit seed positions
    if (next_cycleposedit())
      return 1;

    if (gv_hmove_choice%2 == 0) {
      // B) swap seed elements of non symetric seeds
      if (!gv_symetric && next_swap())
        return 1;
      // C) edit seed
      if (next_edit())
        return 1;
    } else {
      // B) edit seed
      if (next_edit())
        return 1;
      // C) swap seed elements of non symetric seeds
      if (!gv_symetric && next_swap())
        return 1;
    }
  } else {
    if (gv_hmove_choice%3 == 1) {
      // A) swap seed elements of non symetric seeds
      if (!gv_symetric && next_swap())
        return 1;
      if (gv_hmove_choice%2 == 0) {
        // B) edit seed positions
        if (next_cycleposedit())
          return 1;
        // C) edit seed
        if (next_edit())
          return 1;
      } else {
        // B) edit seed
        if (next_edit())
          return 1;
        // C) edit seed positions
        if (next_cycleposedit())
          return 1;
      }
    } else {
      // A) edit seed
      if (next_edit())
        return 1;
      if (gv_hmove_choice%2 == 0) {
        // B) edit seed positions
        if (next_cycleposedit())
          return 1;
        // C) swap seed elements of non symetric seeds
        if (!gv_symetric && next_swap())
          return 1;
      } else {
        // B) swap seed elements of non symetric seeds
        if (!gv_symetric && next_swap())
          return 1;
        // C) edit seed positions
        if (next_cycleposedit())
          return 1;
      }
    }
  }
  return 0;
}



/* Set the last hmove
 */
void seed::set_hmove() {

  _swap_from      = 0;
  _swap_to        = 0;
  _swap_and_edit_shift = rand()%_span;
  _edit_cycle_pos_index = -1;
  _edit_cycle_pos_old_value = 0;
  _edit_operation   = 0;
  _edit_from        = -1;
  _edit_from_symbol = 0;

  if (_seedCyclePos_int) {
    // shuffle positions without changing their values to have a little different order to optimize next time ...
    for (int i = 0; i < _seedNbCyclePos; i++) {
      int j = rand()%_seedNbCyclePos;
      int t                = _seedCyclePos_int[i];
      _seedCyclePos_int[i] = _seedCyclePos_int[j];
      _seedCyclePos_int[j] = t;
    }
  }
}


/* Reset all the  hmoves, come back to the original seed
 */
void seed::reset_hmove() {
  /* A) come back to the original seed */
  /* swap */
  if (!gv_symetric && _edit_operation != 1 &&  _edit_operation != 2 && _swap_from < _span && _swap_to < _span) {
    int letter_tmp                                          = _seedMotif_int[(_swap_from+_swap_and_edit_shift)%_span];
    _seedMotif_int[(_swap_from+_swap_and_edit_shift)%_span] = _seedMotif_int[(_swap_to+_swap_and_edit_shift)%_span];
    _seedMotif_int[(_swap_to+_swap_and_edit_shift)%_span]   = letter_tmp;
  }
  /* edit_pos */
  if (_edit_cycle_pos_index >= 0 && _edit_cycle_pos_index < _seedNbCyclePos) {
    _seedCyclePos_int[_edit_cycle_pos_index] = _edit_cycle_pos_old_value;
  }
  /* edit */
  if (_edit_operation == 1 && _edit_from >= 0) {
    _edit_from--;
    next_delete();
  } else {
    if (_edit_operation == 2 && _edit_from >= 0) {
      _edit_from--;
      next_insert();
    }
  }
  /* B) reset all the values after that */
  set_hmove();
}






  /** @brief set cycle properties of the seed
   *  @param cycle_pos is a vector of positions where the seed is allowed to match
   *  @param cycle_size is the maximal position where a position can be set
   */
void seed::setCyclePos(vector<int> cycle_pos, int cycle_size) {
    _seedCyclePos_int   = new int[cycle_pos.size()];
    for (int i = 0; i < (int)cycle_pos.size(); i++)
      _seedCyclePos_int[i] = cycle_pos[i];
    _seedNbCyclePos     = cycle_pos.size();
    _seedMaxCyclePos    = cycle_size;
  }


  /**
   * @brief convert a char into a number
   * @brief '0' -> 0, '1' -> 1 ... '9' -> 9, 'A' -> 10 ... 'Z' -> 35, 'a' -> 36 ...
   * @param s is the string where conversion is currently processed
   * @param pos is the first char inside s where the conversion has to be done
   * @return the converted value (integer that represents a seed letter) and increments i to the next symbol
   */
  int seed::convert(string & s, int * pos){

    const char * c = s.c_str();
    c += (*pos);

    if (gv_bsymbols_flag) {
      for (int i = 0; i < gv_seed_alphabet_size; i++)
        if (gv_bsymbols_array[i] == *c) {
          (*pos)++;
          return i;
        }
      _ERROR("Single Seed Parsing","invalid char \'" << c << "\' inside the seed descriptor");
    } else {

      if (*c >= '0' && *c <= '9') {
        if (*c - '0' >= gv_seed_alphabet_size) {
          _ERROR("Single Seed Parsing","out of bound char \'" << c << "\' inside the seed descriptor");
        }
        (*pos)++;
        return *c - '0';
      } else if (*c >= 'A' && *c <= 'Z') {
        if (*c - 'A' + 10 >= gv_seed_alphabet_size) {
          _ERROR("Single Seed Parsing","out of bound char \'" << c << "\' inside the seed descriptor");
        }
        (*pos)++;
        return *c - 'A' + 10;
      } else if (*c >= 'a' && *c <= 'z') {
        if (*c - 'a' + 36 >= gv_seed_alphabet_size) {
          _ERROR("Single Seed Parsing","out of bound char \'" << c << "\' inside the seed descriptor");
        }
        (*pos)++;
        return *c - 'a' + 36;
      } else if (*c == '#') {
        int code = 0;
        (*pos)++;
        c++;
        while (*c >= '0' && *c <= '9') {
          code = code * 10 + (*c - '0');
          (*pos)++;
          c++;
        }
        if (code >= 0 && code < gv_seed_alphabet_size) {
          return code;
        } else {
          _ERROR("Single Seed Parsing","out of bound code \'#" << code << "\' inside the seed descriptor");
        }
      } else {
          _ERROR("Single Seed Parsing","out of bound char \'" << c << "\' inside the seed descriptor");
          return 0;
      }
    }
    return 0;
  }

  /// convert a number into a char
string seed::rconvert(int i){
  if (gv_bsymbols_flag) {
    return string(1,gv_bsymbols_array[i]);
  } else {
      if (gv_seed_alphabet_size <= 62) {
        if (i < 10)
          return string(1,(char)i + '0');
        else
          if (i < 36)
            return string(1,(char)i - 10 + 'A');
          else
            return string(1,(char)i - 36 + 'a');
      } else {
        string s = "";
        while (i) {
          char c = '0' + (i%10);
          i /= 10;
          s = string(1,c) + s;
        }
        return string("#" + s);
      }
    }
}



#define MATCHES_AB(a,b)    (matchingmatrix[(a)][(b)])

/// returns first position where a seed hit occurs (-1 if no hit)
int seed::Hit(const vector<int> & alignment, const vector< vector <int> > & matchingmatrix) {
  if ( _seedCyclePos_int ) {
    /* create a copy of pos and sort it */
    vector<int> cp = vector<int>(0);
    for (int p = 0; p < _seedNbCyclePos; p++)
      cp.push_back(((_seedCyclePos_int[p] - (_span%_seedMaxCyclePos) + _seedMaxCyclePos)%_seedMaxCyclePos));
    sort(cp.begin(),cp.end());
    /* enumerate selected compatible positions */
    int u = 0;
    while (int v = cp[u % _seedNbCyclePos] + _seedMaxCyclePos * (u / _seedNbCyclePos) <= (int) alignment.size() - _span) {
      int i = 0;
      while (i < _span && MATCHES_AB(alignment[i+v],_seedMotif_int[i]))
        i++;
      if (i == _span)
        return v;
      u++;
    }
  } else {
    /* enumerate all compatible positions */
    for (int v = 0; v <= (int) alignment.size() - _span; v++) {
      int i = 0;
      while (i < _span && MATCHES_AB(alignment[i+v],_seedMotif_int[i]))
        i++;
      if (i == _span)
        return v;
    }
  }
  return -1;
}

/// returns number of positions where seed hits (0 if no hit)
int seed::mHits(const vector<int> & alignment, const vector< vector <int> > & matchingmatrix) {
  int count = 0;
  if ( _seedCyclePos_int ) {
    /* create a copy of pos and sort it */
    vector<int> cp = vector<int>(0);
    for (int p = 0; p < _seedNbCyclePos; p++)
      cp.push_back(((_seedCyclePos_int[p] - (_span%_seedMaxCyclePos) + _seedMaxCyclePos)%_seedMaxCyclePos));
    sort(cp.begin(),cp.end());
    /* enumerate selected compatible positions */
    int u = 0;
    while (int v = cp[u % _seedNbCyclePos] + _seedMaxCyclePos * (u / _seedNbCyclePos) <= (int) alignment.size() - _span) {
      int i = 0;
      while (i < _span && MATCHES_AB(alignment[i+v],_seedMotif_int[i]))
        i++;
      if (i == _span)
        count++;
      u++;
    }
  } else {
    /* enumerate all compatible positions */
    for (int v = 0; v <= (int) alignment.size() - _span; v++) {
      int i = 0;
      while (i < _span && MATCHES_AB(alignment[i+v],_seedMotif_int[i]))
        i++;
      if (i == _span)
        count++;
    }
  }
  return count;
}
