#ifndef __AUTOMATON_HH__
#define __AUTOMATON_HH__
// @{

// templates to avoid code duplication between "costs" / "double" / "count"

/** @brief Probabilistic/Cost/Count matrix product of the current automaton with a probabilistic automaton
 *  @param other is the second probabilistic automaton used for the product
 *  @param productSetFinalType indicates if product final states are the crossproduct of both automaton final states or only one of these
 *       @li PRODUCT_UNION_xxx           : automata "union" of final states (boolean values),
 *       @li PRODUCT_INTERSECTION_xxx    : automata "intersection" of final states (boolean values)
 *       @li PRODUCT_BUTNOT_xxx          : automata "this".final BUT NOT "other".final (boolean values)
 *       @li PRODUCT_NOTBUT_xxx          : automata NOT "this".final BUT "other".final
 *       @li PRODUCT_UNION_NO_LOOP_ADD   : automata "union" and "sum value" of final states, (integer value)
 *       @li PRODUCT_ADDHOC_xxx          : automata addhoc final function aff does this work
 *          with
 *       @li xxx : LOOP / NO_LOOP        : indicates if the final state is a single one absorbant (force it, otherwise keep it as in the "true" product)
 *  @param depth indicates the maximal depth that must be reached : extra states are non final selflooping states
 *         (by default, this value is greater than 2 Billions, but the given alignment length should be enought in most cases)
 *  @param aff function indicates the final value to be used, is is used when @param productSetFinalType = PRODUCT_ADDHOC_xxx
 *  @return a new matrix that only gives reachable states of the product.
 *  @see product
 *  @see matrices_step_product
 *  @warning use  PRODUCT_UNION_NO_LOOP_ADD  when manipulating integer values for final states
 *  @relates automaton
 */

template<typename T> inline matrix<T> * matrix_product(const automaton & other, const ProductSetFinalType productSetFinalType, const int depth, const AddHoc_Final_Func aff = NULL) const {

  VERB_FILTER(VERBOSITY_DEBUGGING, INFO__("== Matrix_Product == (productsetfinaltype:" << dec << productSetFinalType << ")"););

#ifdef ASSERTB
  if ((this->size() * other.size()) > (1<<28)) {
    _ERROR("matrix_product"," size of product automaton will \"certainly\" explode : better stop here ...");
  }
#endif

  matrix<T> * result = new matrix<T>();
  int i = result->addNewRow(TRUE);
  result->selfLoop(i,One<T>());

#ifdef USEMAPPRODUCT
  typedef less< pair<int,int> > lessp;
  typedef map< pair<int,int> , int, lessp > maptype;
  maptype statesNbIndex;
#define PRODINDEX(i) (i)
#else
  vector<int> statesNbIndex( this->size() * other.size(), 0);
#define PRODINDEX(i) ((i).first * other.size() + (i).second)
#endif

#ifdef USEQUEUEPRODUCT
  queue< pair<int,int> >  statesNbRemaining;
#else
  stack< pair<int,int> >  statesNbRemaining;
#endif

  // (0) final loop or not
  int final_loop = (productSetFinalType == PRODUCT_UNION_FINAL_LOOP || productSetFinalType == PRODUCT_INTERSECTION_FINAL_LOOP || productSetFinalType == PRODUCT_BUTNOT_FINAL_LOOP || productSetFinalType == PRODUCT_NOTBUT_FINAL_LOOP || productSetFinalType == PRODUCT_ADDHOC_FINAL_LOOP) ? TRUE : FALSE;

  // (1) start the product init state
  pair<int,int> indexN             = pair<int,int>(1,1);
  int stateNumber                  = result->addNewRow();

  statesNbIndex[PRODINDEX(indexN)] = stateNumber;
  statesNbRemaining.push(indexN);

#ifdef USEQUEUEPRODUCT
  // depth of the states being built
  int level_i                    = 0;
  int stateN_of_level_i          = stateNumber;
#endif

  // (2) take all the non-considered interesting cases pairs on both automatons
  while(!statesNbRemaining.empty()) {

    // current state remaining
#ifdef USEQUEUEPRODUCT
    pair<int,int> indexN  =  statesNbRemaining.front();
#else
    pair<int,int> indexN  =  statesNbRemaining.top();
#endif
    statesNbRemaining.pop();

    int stateNA = indexN.first;
    int stateNB = indexN.second;
    int stateN =  statesNbIndex[PRODINDEX(indexN)];

    VERB_FILTER(VERBOSITY_DEBUGGING, INFO__("$pop state:" << stateN););

    for(int a = 0; a < gv_align_alphabet_size; a++) {

      VERB_FILTER(VERBOSITY_DEBUGGING, INFO__("a = " << a););

      for(vector<transition>::const_iterator iterA = _states[stateNA]._next[a].begin(); iterA != _states[stateNA]._next[a].end(); iterA++) {
        for(vector<transition>::const_iterator iterB = other._states[stateNB]._next[a].begin(); iterB != other._states[stateNB]._next[a].end(); iterB++) {

          int stateAnext = iterA->_state;
          int stateBnext = iterB->_state;
          pair<int,int> indexNx = pair<int,int>(stateAnext,stateBnext);

          int stateNx    = 0;

          // final state
          int final_state = FALSE;


          switch (productSetFinalType) {

          case PRODUCT_UNION_FINAL_LOOP:
          case PRODUCT_UNION_NO_FINAL_LOOP:
            final_state =  ((this->_states[stateAnext]._final) || (other._states[stateBnext]._final)) ? TRUE : FALSE;
            break;
          case PRODUCT_UNION_NO_FINAL_LOOP_ADD:
            final_state =  ((this->_states[stateAnext]._final) + (other._states[stateBnext]._final));
            break;

          case PRODUCT_INTERSECTION_FINAL_LOOP:
          case PRODUCT_INTERSECTION_NO_FINAL_LOOP:
            final_state =  ((this->_states[stateAnext]._final) && (other._states[stateBnext]._final)) ? TRUE : FALSE;
            break;

          case PRODUCT_BUTNOT_FINAL_LOOP:
          case PRODUCT_BUTNOT_NO_FINAL_LOOP:
            final_state =  ((this->_states[stateAnext]._final) && (!(other._states[stateBnext]._final))) ? TRUE : FALSE;
            break;

          case PRODUCT_NOTBUT_FINAL_LOOP:
          case PRODUCT_NOTBUT_NO_FINAL_LOOP:
            final_state =  ((!(this->_states[stateAnext]._final)) && (other._states[stateBnext]._final)) ? TRUE : FALSE;
            break;

          case PRODUCT_ADDHOC_FINAL_LOOP:
          case PRODUCT_ADDHOC_NO_FINAL_LOOP:
            final_state =  aff(this->_states[stateAnext]._final,other._states[stateBnext]._final);
            break;

          }

          // add the "new" state, considering booleans "final_state" and "final_state"
          if (final_state && final_loop) {
            stateNx = 0;
          } else {
            if  (!statesNbIndex[PRODINDEX(indexNx)]) {

#ifdef USEQUEUEPRODUCT
              // compute level
              if (stateN > stateN_of_level_i) {
                stateN_of_level_i = (result->size() - 1);
                level_i++;
              }

              if (level_i <= depth ){
#endif
                // create a new state
                stateNx = result->addNewRow(final_state);
                statesNbIndex[PRODINDEX(indexNx)] = stateNx;
                statesNbRemaining.push(indexNx);

                VERB_FILTER(VERBOSITY_DEBUGGING, INFO__("$push state:" << dec << stateNx););

#ifdef USEQUEUEPRODUCT
              } else {
                // max level reached : goes to a "non final" loop state
                stateNx = result->addNewRow();
                result->selfLoop(stateNx,One<T>());
                }
#endif
            } else {
              stateNx = statesNbIndex[PRODINDEX(indexNx)];
            }
          }

          VERB_FILTER(VERBOSITY_DEBUGGING, INFO__("> add transition ( a:" << dec << a << ", q1:" << dec << stateN << ", q2:" << stateNx << " ) "););

          // add a transition on from stateN --a--> stateNx.
          T t = Transition<T>(a,iterB->_prob);
          result->addNewCell(stateN,stateNx,t);
        }// for (iterB)
      }// for (iterA)
    }// for (a)
  }//while stack nonempty

  // Free unused data needed to build the automaton
  statesNbIndex.clear();
  return result;
};

/** @brief Probabilistic/Cost matrix slice product of the current automaton with a probabilistic automaton
 *  @brief the slice product compute the "depth" first matrices
 *  @param other is the second probabilistic automaton used for the product
 *  @param productSetFinalType indicates if product final states are the crossproduct of both automaton final states or only one of these
 *       @li PRODUCT_UNION_xxx           : automata "union" of final states,
 *       @li PRODUCT_INTERSECTION_xxx    : automata "intersection" of final states
 *       @li PRODUCT_BUTNOT_xxx          : automata "this".final BUT NOT "other".final
 *       @li PRODUCT_NOTBUT_xxx          : automata NOT "this".final BUT "other".final
 *       @li PRODUCT_UNION_NO_LOOP_ADD   : automata "union" and "sum value" of final states, (integer value)
 *       @li PRODUCT_ADDHOC_xxx          : automata addhoc final function aff does this work
 *          with
 *       @li xxx : LOOP / NO_LOOP        : indicates if the final state is a single one absorbant (force it, otherwise keep it as in the "true" product)
 *  @param depth indicates the maximal depth that must be reached
 *  @param aff function indicates the final value to be used, is is used when @param productSetFinalType = PRODUCT_ADDHOC_xxx
 *  @return an array (size depth) of matrices such that each matrix gives the transition from any state to the next ...
 *  @see product
 *  @see matrix_product
 *  @relates automaton
 */

template<typename T> inline vector< matrix<T> * > * matrices_step_product(const automaton & other, const ProductSetFinalType productSetFinalType, const int depth, const AddHoc_Final_Func aff = NULL) const {

  VERB_FILTER(VERBOSITY_DEBUGGING, INFO__("== Matrix_Step_Product == (productsetfinaltype:" << dec << productSetFinalType << ")"););

#ifdef ASSERTB
  if ((this->size() * other.size()) > (1<<28)) {
    _ERROR("matrices_step_product"," size of product automaton will \"certainly\" explode : better stop here ...");
  }
#endif
  vector < matrix<T> * > * result = new vector< matrix<T> * > (depth+1,(matrix<T> *) NULL);

  double Bprob[2][other.size()];
  for (int i = 0; i < other.size(); i++)
    Bprob[0][i] = 0.0;
  Bprob[0][1] = 1.0;


  for(int i = 0; i <= depth; i++) {
    (*result)[i] = new matrix<T>();
    // state 0 on all the "matrices"
    int p = ((*result)[i])->addNewRow(TRUE);
    ((*result)[i])->selfLoop(p,One<T>());
    // state 1 (init)
    ((*result)[i])->addNewRow(FALSE);
  }

  // FIXME : only map is used
  typedef less< pair<int,int> > lessp;
  typedef map< pair<int,int> , int, lessp > maptype;
  maptype statesNbIndex[2];

  queue< pair<int,int> >  statesNbRemaining[2];

  // (0) final loop or not
  int final_loop = (productSetFinalType == PRODUCT_UNION_FINAL_LOOP || productSetFinalType == PRODUCT_INTERSECTION_FINAL_LOOP || productSetFinalType == PRODUCT_BUTNOT_FINAL_LOOP || productSetFinalType == PRODUCT_NOTBUT_FINAL_LOOP || productSetFinalType == PRODUCT_ADDHOC_FINAL_LOOP) ? TRUE : FALSE;

  // (1) build a matrix per level
  for(int level_i = 0; level_i < depth; level_i++) {

    // (1.a) build the remaing stats transitions from $level_i$ to $level_{i+1}$  (none when $level_i = 0$)
    while(!statesNbRemaining[level_i%2].empty()) {

      // current state remaining
      pair<int,int> indexN  =  statesNbRemaining[level_i%2].front();
      statesNbRemaining[level_i%2].pop();

      int stateNA = indexN.first;
      int stateNB = indexN.second;
      int stateN  = statesNbIndex[level_i%2][indexN];

      VERB_FILTER(VERBOSITY_DEBUGGING, INFO__("$pop state:" << stateN););

      for(int a = 0; a < gv_align_alphabet_size; a++) {

        for(vector<transition>::const_iterator iterA = _states[stateNA]._next[a].begin(); iterA != _states[stateNA]._next[a].end(); iterA++) {
          for(vector<transition>::const_iterator iterB = other._states[stateNB]._next[a].begin(); iterB != other._states[stateNB]._next[a].end(); iterB++) {

            int stateAnext = iterA->_state;
            int stateBnext = iterB->_state;
            pair<int,int> indexNx = pair<int,int>(stateAnext,stateBnext);

            int stateNx    = 0;

            // final state
            int final_state = FALSE;

            switch (productSetFinalType) {

            case PRODUCT_UNION_FINAL_LOOP:
            case PRODUCT_UNION_NO_FINAL_LOOP:
              final_state =  ((this->_states[stateAnext]._final) || (other._states[stateBnext]._final)) ? TRUE : FALSE;
              break;
            case  PRODUCT_UNION_NO_FINAL_LOOP_ADD:
              final_state =  ((this->_states[stateAnext]._final) + (other._states[stateBnext]._final));
              break;

            case PRODUCT_INTERSECTION_FINAL_LOOP:
            case PRODUCT_INTERSECTION_NO_FINAL_LOOP:
              final_state =  ((this->_states[stateAnext]._final) && (other._states[stateBnext]._final)) ? TRUE : FALSE;
              break;

            case PRODUCT_BUTNOT_FINAL_LOOP:
            case PRODUCT_BUTNOT_NO_FINAL_LOOP:
              final_state =  ((this->_states[stateAnext]._final) && (!(other._states[stateBnext]._final))) ? TRUE : FALSE;
              break;

            case PRODUCT_NOTBUT_FINAL_LOOP:
            case PRODUCT_NOTBUT_NO_FINAL_LOOP:
              final_state =  ((!(this->_states[stateAnext]._final)) && (other._states[stateBnext]._final)) ? TRUE : FALSE;
              break;

            case PRODUCT_ADDHOC_FINAL_LOOP:
            case PRODUCT_ADDHOC_NO_FINAL_LOOP:
              final_state =  aff(this->_states[stateAnext]._final,other._states[stateBnext]._final);
              break;
            }

            // add the "new" state, considering booleans "final_state" and "final_state"
            if (final_state && final_loop) {
              stateNx = 0;
            } else {
              if  (!statesNbIndex[(level_i+1)%2][indexNx]) {
                stateNx                               = (*result)[level_i+1]->addNewRow(final_state);
                statesNbIndex[(level_i+1)%2][indexNx] = stateNx;
                statesNbRemaining[(level_i+1)%2].push(indexNx);

                VERB_FILTER(VERBOSITY_DEBUGGING, INFO__("$push state:" << dec << stateNx););

              } else {
                stateNx = statesNbIndex[(level_i+1)%2][indexNx];
              }
            }

            VERB_FILTER(VERBOSITY_DEBUGGING, INFO__("> add transition ( a:" << dec << a << ", q1:" << dec << stateN << ", q2:" << stateNx << " ) "););

            // add a transition on from stateN --a--> stateNx.
            T t = Transition<T>(a,iterB->_prob);
            (*result)[level_i]->addNewCell(stateN,stateNx,t);
          }// for (iterB)
        }// for (iterA)
      }// for (a)
    }//while stack nonempty


    // clean old map
    statesNbIndex[level_i%2].clear();


    // (1.b) build the marginal states transitions and distribution
    int stateNA = (_init_states.size() ? _init_states[level_i%(_init_states.size())] : 1);
    // start the product init state for the B where P(B) > 0
    for (int stateNB = 0; stateNB < other.size(); stateNB++) {
      if (Bprob[level_i%2][stateNB] > 0 ) {

        for(int a = 0; a < gv_align_alphabet_size; a++) {

          for(vector<transition>::const_iterator iterA = _states[stateNA]._next[a].begin(); iterA != _states[stateNA]._next[a].end(); iterA++) {
            for(vector<transition>::const_iterator iterB = other._states[stateNB]._next[a].begin(); iterB != other._states[stateNB]._next[a].end(); iterB++) {

              int stateAnext = iterA->_state;
              int stateBnext = iterB->_state;
              pair<int,int> indexNx = pair<int,int>(stateAnext,stateBnext);
              int stateNx           = statesNbIndex[(level_i+1)%2][indexNx];
              if  (!stateNx) {
                stateNx  = (*result)[level_i+1]->addNewRow(this->_states[stateAnext]._final);
                statesNbIndex[(level_i+1)%2][indexNx] = stateNx;
                statesNbRemaining[(level_i+1)%2].push(indexNx);
              }
              T t = Transition<T>(a,Bprob[level_i%2][stateNB]*iterB->_prob);
              (*result)[level_i]->addNewCell(1,stateNx,t);
            }
          }
        }
      }
    }

    // (1.c) update the marginal probabilities for B
    for (int stateNB = 0; stateNB < other.size(); stateNB++) {
      Bprob[(level_i+1)%2][stateNB] = 0.0;
    }
    for (int stateNB = 0; stateNB < other.size(); stateNB++) {
      for(int a = 0; a < gv_align_alphabet_size; a++) {
        for(vector<transition>::const_iterator iterB = other._states[stateNB]._next[a].begin(); iterB != other._states[stateNB]._next[a].end(); iterB++) {
          Bprob[(level_i+1)%2][iterB->_state] +=  Bprob[level_i%2][stateNB] * iterB->_prob;
        }
      }
    }


  } // for(level_i)

  // Free unused data needed to build the automaton
  statesNbIndex[(depth+1)%2].clear();
  //statesNbRemaining[level_i%2].clear(); FIXME

  return result;
};

// @}

#endif
