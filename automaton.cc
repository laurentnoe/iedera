#include "automaton.h"
#include <algorithm>

//#define ASSERTB
//#define BUILD

#ifdef ASSERTB
#include <assert.h>
#endif

#define DETERMINISTIC(s,a) (this->_states[s]._next[a].size() == 1)
#define MATCHES_AB(a,b)    (matchingmatrix[(a)][(b)])
#define SCORES_AB(a,b)     (scoringmatrix[(a)][(b)])


// output method for the current automaton
ostream& operator<<(ostream& os, const automaton& automaton) {
  os << automaton._states.size()  << endl;
  // display each state
  for (int i = 0; i < (int)automaton._states.size(); i++){
    os << dec << i << "\t" << automaton._states[i]._final << endl;
    // display for each state forward transitions
    for (int a = 0; a < gv_align_alphabet_size; a++) {
      os << "\t" << dec << a << "\t" << dec << automaton._states[i]._next[a].size()  << endl;
      for (vector<transition>::const_iterator iter = automaton._states[i]._next[a].begin(); iter != automaton._states[i]._next[a].end(); iter++) {
        os << "\t\t"  << dec <<  (iter->_state) << "\t" << (iter->_prob) << endl;
      }
    }
  }
  return os;
}


istream& operator>>(istream& is, automaton& automaton) {

  // clear previous automaton
  for (int i = 0; i < (int)automaton._states.size(); i++){
    automaton._states[i].clear();
  }
  automaton._states.clear();

  // read automaton size
  int size = 0;
  is >> size;
#ifdef DEBUGREADING
  cerr << size << endl;
#endif
  if (size <= 0) {
    cerr << "> when reading automaton size" << endl;
    _ERROR("operator>>","incorrect size "<< size);
  }
  // add states first
  for (int i = 0; i < size; i++){
    automaton.addNewState();
  }

  // display each state
  for (int statefrom = 0; statefrom < size; statefrom++){
    int state = 0, final = 0;
    is >>  state >>  final;
#ifdef DEBUGREADING
    cerr << state << "\t" << final << endl;
#endif

    if (state < 0 || state >= size) {
      cerr << "> when reading automaton state " << state << endl;
      _ERROR("operator>>","incorrect state " << state  << " (not in [0.."<<(size-1)<<"])");
    }
    if (final != 0 && final != 1) {
      cerr << "> when reading automaton state " << state << endl;
      _ERROR("operator>>","incorrect state final flag " << final  << " (not in [0..1])");
    }

    automaton._states[statefrom]._final = final;

    // display for each state forward transitions
    double probability_sum = 0.0;
    for (int a = 0; a < gv_align_alphabet_size; a++){
      int x = 0, alistsize = 0;
      is >> x >> alistsize;
#ifdef DEBUGREADING
      cerr << "\t" << x << "\t" <<  alistsize << endl;
#endif

      if (x != a) {
        cerr << "> when reading automaton state " << state << endl;
        _ERROR("operator>>","incorrect alphabet letter " << x << " (should be ordered and thus be " << a << ")");
      }
      if (alistsize < 0 ) {
        cerr << "> when reading automaton state " << state << endl;
        _ERROR("operator>>","incorrect transition list size "<< alistsize << " for letter " << x << " (should be >= 0)");
      }

      for (int j = 0; j < alistsize; j++){
        int    stateto = 0;
        double probability  = 0.0;
        is >>  stateto >> probability;
#ifdef DEBUGREADING
        cerr << "\t\t" << stateto << "\t" << probability << endl;
#endif
        probability_sum += probability;

        if (stateto < 0 || stateto >= size) {
          cerr << "> when reading automaton state " << state << endl;
          _ERROR("operator>>","incorrect transition to state "<< stateto << " for letter " << x << " (this state does not exist)");
        }
        automaton.addNewTransition(a,statefrom,stateto,probability);
      }
    }
    if ( probability_sum > 1.0001 || probability_sum < 0.9999 ){
      cerr << "> when reading state \"" << state << "\"  letters " << endl;
      _ERROR("operator>>","incorrect probability sum "<< probability_sum);
    }
  }
  return is;
}


// "dot" format output (use "graphviz\dotty" or "graphviz\dot file.dot -T ps -o file.eps" to visualise it)

void automaton::dot(ostream& os) const {

  os << "digraph A {"  << endl;
  os << "\t fontsize =\"8\"" << endl;

  // display each state
  for (int i = 0; i < (int)_states.size(); i++) {
    os << "\t \"node" << i << "\" [ label = \"" << i;
    if (_states[i]._final) {
      os << "(" << (_states[i]._final) << ")";
    }
    os << "\",shape = " << ((_states[i]._final)?"doublecircle":"circle") << " ];"<<endl;
  }

  // display forward transitions
  for (int i = 0; i < (int)_states.size(); i++) {
    for (int a = 0; a < gv_align_alphabet_size; a++) {
      for (vector<transition>::const_iterator iter = _states[i]._next[a].begin(); iter != _states[i]._next[a].end(); iter++) {
        os << "\t \"node" << i <<"\"  -> \"node" << (iter->_state) << "\" [label = \"" << a << //" (" << (iter->_prob) << ")" <<
          "\"];" << endl;
      }
    }
  }
  os << "}" << endl;
}

// gap-system Fr package (use "LoadPackage("FR");" inside gap )
// http://www.gap-system.org/Manuals/pkg/fr-2.1.1/doc/chap0.html

void automaton::gapFR(ostream& os) const {
  os << "MealyElement(";
  // (1/4) domain
  os << "Domain([";
  for (int a = 0; a < gv_align_alphabet_size; a++) {
    if (a > 0)
      os << ",";
    os << (a+1);
  }
  os << "]),";
  // (2/4) transitions
  os << "[";
  for (int i = 0; i < (int)_states.size(); i++) {
    if (i > 0)
        os << ",";
    os << "[";
    for (int a = 0; a < gv_align_alphabet_size; a++) {
      if (a > 0)
        os << ",";
      for (vector<transition>::const_iterator iter = _states[i]._next[a].begin(); iter != _states[i]._next[a].end(); iter++) {
        if (iter != _states[i]._next[a].begin()) {
          _ERROR("gapFR","non determinitic automaton");
        }
        os << (iter->_state+1);
      }
    }
    os << "]";
  }
  os << "],";
  // (3/4) output {ouput final states on previous transitions to have Mealy}
  os << "[";
  for (int i = 0; i < (int)_states.size(); i++) {
    if (i > 0)
      os << ",";
    os << "[";
    for (int a = 0; a < gv_align_alphabet_size; a++) {
      if (a > 0)
        os << ",";
      for (vector<transition>::const_iterator iter = _states[i]._next[a].begin(); iter != _states[i]._next[a].end(); iter++) {
        if (iter != _states[i]._next[a].begin()) {
          _ERROR("gapFR","non determinitic automaton");
        }
        os << (_states[iter->_state]._final+1);
      }
    }
    os << "]";
  }
  os << "],";
  // (4/4) init
  os << "2";
  os << ");" << endl;
}


// gap-system Automata package (use "LoadPackage("Automata");" inside gap )
// http://www.gap-system.org/Manuals/pkg/automata/doc/chap0.html

void automaton::gapAutomata(ostream& os) const {
  bool deterministic = true;
  // (0/6) check in deterministic or not
  for (int i = 0; i < (int)_states.size(); i++) {
    for (int a = 0; a < gv_align_alphabet_size; a++) {
      if (!DETERMINISTIC(i,a)) {
        deterministic = false;
         goto end_check_det;
      }
    }
  }
 end_check_det:
  // (1/6) det/nondet
  os << "Automata(";
  if (deterministic)
    os << "\"det\"," << endl;
  else
    os << "\"nondet\"," << endl;
  // (2/6) number of states
  os << (_states.size()) << ",";
  // (3/6) alphabet size
  os << gv_align_alphabet_size << "," << endl;
  // (4/6) transitions
  os << "[";
  for (int i = 0; i < (int)_states.size(); i++) {
    if (i > 0)
      os << ",";
    os << "[";
    for (int a = 0; a < gv_align_alphabet_size; a++) {
      if (a > 0)
        os << ",";
      // zero or one transition : zero gives ",," wheras one give ",<int>,"
      if (_states[i]._next[a].size() <= 1) {
        if (_states[i]._next[a].begin() != _states[i]._next[a].end())
          os << (_states[i]._next[a].begin()->_state+1);
      } else {
        // more than one transition gives ",[<int>,<int>,...],
        os << "[";
        for (vector<transition>::const_iterator iter = _states[i]._next[a].begin(); iter != _states[i]._next[a].end(); iter++) {
          if (iter != _states[i]._next[a].begin())
            os << ",";
          os << (iter->_state+1);
        }
        os << "]";
      }
    }
    os << "]";
  }
  os << "]," << endl;
  // (5/6) init states
  os << "[2]," << endl;
  // (6/6) final states
  bool outfinal = false;
  for (int i = 0; i < (int)_states.size(); i++) {
    os << "[";
    if (_states[i]._final) { // FIXME NO WAY TO FIX SEVERAL SETS OF FINAL STATES ?
      if (outfinal)
        os << ",";
      outfinal = true;
      os << (i+1);
    }
    os << "]" << endl;
  }
  os << ");" << endl;
}


// This method create a new state and return number of states.
inline int automaton::addNewState(const int final) {
  _states.push_back(state(final));
  return _states.size() - 1;
}


// This method add a transition between two states on letter "a"
void automaton::addNewTransition(const int a, const int startingState, const int endingState, const double prob) {

#ifdef ASSERTB
  if (a < 0 || a >= gv_align_alphabet_size) {
    _ERROR("addNewTransition","transition letter out of range");
  }

  if (startingState >= (int)_states.size() || endingState >= (int)_states.size()) {
    _ERROR("addNewTransition","state required does not exist");
  }
#endif

  // forward linking
#ifdef ASSERTB
  // do not add twice the same state in the "a" backward list
  for (vector<transition>::const_iterator iter = _states[startingState]._next[a].begin(); iter != _states[startingState]._next[a].end(); iter++) {
    if ( iter->_state == endingState) {
      cerr << "> when linking state q2:" << dec << endingState << " with state q1:" << dec << startingState << " on letter a:" << dec << a << endl;
      _ERROR("addNewTransition"," state q2 has already a transition on \"a\" ");
    }
  }
#endif
  _states[startingState]._next[a].push_back(transition(endingState,prob));
  ///@todo{FIXME: push_front ?}
}


// This method modifies the transition probability (a,startingState,endingState)
int automaton::changeTransitionProb(const int a, const int startingState, const int endingState, const double prob) {

#ifdef ASSERTB
  if (a < 0 || a >= gv_align_alphabet_size) {
    _ERROR(":changeTransitionProb","transition letter out of range");
  }

  if (startingState >= (int)_states.size() || endingState >= (int)_states.size()) {
    _ERROR(":changeTransitionProb","state required does not exist");
  }
#endif


  // forward link pr
  for (vector<transition>::iterator iter = _states[startingState]._next[a].begin(); iter != _states[startingState]._next[a].end(); iter++) {
    if (iter->_state == endingState) {
      iter->_prob = prob;
      break;
    }
  }
  return 0;
}



// Give the current probability of transition
double automaton::Prob(const int a, const int startingState, const int endingState) const {

  // forward link pr
  for (vector<transition>::const_iterator iter = _states[startingState]._next[a].begin(); iter != _states[startingState]._next[a].end(); iter++) {
    if ( iter->_state == endingState) {
      return iter->_prob;
    }
  }
  return 0.0;
}



// This method selfloops the given state on all letters
void automaton::selfLoop(const int stateNb) {
  for (int a = 0; a < gv_align_alphabet_size; a++)
    addNewTransition(a,stateNb,stateNb,(+1e+0)/gv_align_alphabet_size);
}




// Checks if there is at least on transition startingState --(a)--> ...
bool automaton::hasTransition( int a, int startingState) const {
#ifdef ASSERTB
  assert(startingState >= 0);
  assert(startingState < (int)_states.size());
  assert(a >= 0);
  assert(a < gv_align_alphabet_size);
#endif
  return _states[startingState]._next[a].size() > 0;
}


// Automaton product
automaton * automaton::product(const automaton & other, const ProductSetFinalType productSetFinalType, const ProductProbabilityType thisOrOtherIsProbabilist, const int depth, const AddHoc_Final_Func aff, const int shift) const {

  VERB_FILTER(VERBOSITY_DEBUGGING, INFO__("== Product == (productsetfinaltype:" << dec << productSetFinalType << ")"););

#ifdef ASSERTB
  if ((this->size() * other.size()) > (1<<28)) {
    _ERROR("product"," size of product automaton will \"certainly\" explode : better stop here ...");
  }
#endif

  automaton * result = new automaton();
  int StateFinal = result->addNewState(TRUE);
  result->selfLoop(StateFinal);

#ifdef USEMAPPRODUCT
  typedef less< pair<int,int> > lessp;
  typedef map< pair<int,int>, int, lessp > maptype;
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
  pair<int,int> indexInit             = pair<int,int>(1,1);
  int stateInit                    = result->addNewState();

  statesNbIndex[PRODINDEX(indexInit)] = stateInit;
  statesNbRemaining.push(indexInit);

  if (
#ifndef NOMATRIX
      gv_subalignment_flag && ((this->_init_states.size() > 0) ||
#endif
                               (other._init_states.size() > 0)
#ifndef NOMATRIX
                               )
#endif
      ) {
    result->_init_states.push_back(stateInit);
  }

#ifdef USEQUEUEPRODUCT
  // depth of the states being built
  int level_i                    = 0;
  int stateN_of_level_i          = stateInit;
#endif

  // (2) take all the non-considered interesting cases pairs on both automatons
  while (!statesNbRemaining.empty()) {

    // current state remaining
#ifdef USEQUEUEPRODUCT
    pair<int,int> indexN  =  statesNbRemaining.front();
#else
    pair<int,int> indexN  =  statesNbRemaining.top();
    #warning "the automaton product function is not compiled with a Queue, but with a Stack :  \"depth\" parameter cannot be used on this implementation"
#endif

    statesNbRemaining.pop();

    int stateNA = indexN.first;
    int stateNB = indexN.second;
    int stateN =  statesNbIndex[PRODINDEX(indexN)];

#ifdef USEQUEUEPRODUCT
    // compute level
    if (stateN > stateN_of_level_i) {
      stateN_of_level_i = (result->size() - 1);
      level_i++;
    }
#endif

    VERB_FILTER(VERBOSITY_DEBUGGING, INFO__("$pop state:" << stateN););

    // compute normalizing factor if both are probabilist automata
    double normalize_this_other = 0.0;
    if (thisOrOtherIsProbabilist == PRODUCT_BOTH_ARE_PROBABILIST) {
      for (int a = 0; a < gv_align_alphabet_size; a++) {
        for (vector<transition>::const_iterator iterA = _states[stateNA]._next[a].begin(); iterA != _states[stateNA]._next[a].end(); iterA++) {
          for (vector<transition>::const_iterator iterB = other._states[stateNB]._next[a].begin(); iterB != other._states[stateNB]._next[a].end(); iterB++) {
            normalize_this_other += iterA->_prob * iterB->_prob;
          }
        }
      }
    }

    for (int a = 0; a < gv_align_alphabet_size; a++) {

      VERB_FILTER(VERBOSITY_DEBUGGING, INFO__("a = " << a););

      for (vector<transition>::const_iterator iterA = _states[stateNA]._next[a].begin(); iterA != _states[stateNA]._next[a].end(); iterA++) {
        for (vector<transition>::const_iterator iterB = other._states[stateNB]._next[a].begin(); iterB != other._states[stateNB]._next[a].end(); iterB++) {

          int stateAnext = iterA->_state;
          int stateBnext = iterB->_state;
          pair<int,int> indexNx  = pair<int,int>(stateAnext, stateBnext);

#ifdef USEQUEUEPRODUCT

          // shifter rules
          if (shift) {
            if (shift > 0) {
              if (level_i < shift) {
                stateBnext = 1;
                indexNx    = pair<int,int>(stateAnext, stateBnext);
              }
            } else {
              if (level_i < -shift) {
                stateAnext = 1;
                indexNx    = pair<int,int>(stateAnext, stateBnext);
              }
            }
          }
#else
#warning "the automaton product function is not compiled with Queue, but with a Stack : \"shift\"  parameters cannot be  used on this implementation"
#endif






          int stateNx    = 0;

          // final state
          int final_state = 0;

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
            stateNx = StateFinal;
          } else {
            if  (!statesNbIndex[PRODINDEX(indexNx)]) {

#ifdef USEQUEUEPRODUCT
              if (level_i <= depth) {
#endif
                // create a new state
                stateNx = result->addNewState(final_state);
                statesNbIndex[PRODINDEX(indexNx)] = stateNx;
                statesNbRemaining.push(indexNx);

                if (
#ifndef NOMATRIX
                    gv_subalignment_flag && ((this->_init_states.size() > 0) ||
#endif
                                             (other._init_states.size() > 0)
#ifndef NOMATRIX
                                             )
#endif
                    &&
                    (((this->_init_states.size() > 0) && (stateAnext == this->_init_states[result->_init_states.size()%(this->_init_states.size())])) || ((this->_init_states.size() == 0) && (stateAnext == 1)))
                    &&
                    (((other._init_states.size() > 0) && (stateBnext == other._init_states[result->_init_states.size()%(other._init_states.size())])) || ((other._init_states.size() == 0) && (stateBnext == 1)))
                    )
                  {
                    result->_init_states.push_back(stateNx);
                  }

                VERB_FILTER(VERBOSITY_DEBUGGING, INFO__("$push state:" << dec << stateNx););

#ifdef USEQUEUEPRODUCT

              } else {
                // max level reached : goes to a "non final" loop state
                stateNx = result->addNewState(final_state);
                result->selfLoop(stateNx);
              }
#endif
            } else {
              stateNx = statesNbIndex[PRODINDEX(indexNx)];
            }
          }

          VERB_FILTER(VERBOSITY_DEBUGGING, INFO__("> add transition ( a:" << dec << a << ", q1:" << dec << stateN << ", q2:" << stateNx << " ) "););

          // add a transition on from stateN --a--> stateNx.
          switch (thisOrOtherIsProbabilist) {
          case PRODUCT_NONE_IS_PROBABILIST:
            result->addNewTransition(a,stateN,stateNx);
            break;
          case PRODUCT_THIS_IS_PROBABILIST:
            result->addNewTransition(a,stateN,stateNx,iterA->_prob);
            break;
          case PRODUCT_OTHER_IS_PROBABILIST:
            result->addNewTransition(a,stateN,stateNx,iterB->_prob);
            break;
          case PRODUCT_BOTH_ARE_PROBABILIST:
            result->addNewTransition(a,stateN,stateNx,iterA->_prob * iterB->_prob / normalize_this_other);
            break;
          }
        }// for (listB)
      }// for (listA)
    }// for (a)
  }//while stack nonempty

  // Free unused data needed to build the automaton
  statesNbIndex.clear();
  return result;
}




// This method computes the m hits automaton of the previous one (see product, it is quite similar)
automaton * automaton::mhit(unsigned int m, int depth) const {

  VERB_FILTER(VERBOSITY_DEBUGGING, INFO__("== MHit == (m:" << dec << m << ")"););

  automaton * result = new automaton();
  result->addNewState(TRUE);
  result->selfLoop(0);

  // pair <int,int>
  //   - first int is the state on the original automaton "this"
  //   - second is the number of finals "added to go to this state"


#ifdef USEMAPMHIT
  typedef map< pair<int,int>, int> maptype;
  maptype statesNbIndex;
#define MHITINDEX(i) (i)
#else
  vector<int> statesNbIndex( this->size() * m, 0);
#define MHITINDEX(i) ((i).first * m + (i).second)
#endif

#ifdef USEQUEUEMHIT
  queue< pair<int,unsigned int> >  statesNbRemaining;
#else
  stack< pair<int,unsigned int> >  statesNbRemaining;
#endif

  // (1) start the mhits init state
  pair<int,unsigned int> indexN    = pair<int,unsigned int>(1,0);
  int stateNumber                  = result->addNewState();
  statesNbIndex[MHITINDEX(indexN)] = stateNumber;
  statesNbRemaining.push(indexN);

#ifndef NOMATRIX
  if (gv_subalignment_flag && (this->_init_states.size() > 0)) {
    result->_init_states.push_back(1);
  }
#endif

#ifdef USEQUEUEMHIT
  // depth of the states beeing built
  int level_i                    = 0;
  int stateN_of_level_i          = stateNumber;
#endif

  // (2) take the non considered pair
  while (!statesNbRemaining.empty()) {

    // current state remaining
#ifdef USEQUEUEMHIT
    pair<int,unsigned int> indexN  =  statesNbRemaining.front();
#else
    pair<int,unsigned int> indexN  =  statesNbRemaining.top();
#endif
    statesNbRemaining.pop();

    int thisStateN = indexN.first;
    int     mhitN  = indexN.second;
    int     stateN = statesNbIndex[MHITINDEX(indexN)];

    VERB_FILTER(VERBOSITY_DEBUGGING, INFO__("$pop state:" << stateN););

    for (int a = 0; a < gv_align_alphabet_size; a++) {

      VERB_FILTER(VERBOSITY_DEBUGGING, INFO__("a = " << a););

      for (vector<transition>::const_iterator iter = _states[thisStateN]._next[a].begin(); iter != _states[thisStateN]._next[a].end(); iter++) {

        int thisStateNextN  = iter->_state;
        unsigned int mhitNx = mhitN + this->_states[thisStateNextN]._final;

        int stateNx         = 0;

        // add the "new" state, considering booleans "final_state" and "final_state"
        if (mhitNx >= m) {
          stateNx = 0;
        } else {
          pair<int, unsigned int> indexNx = pair<int, unsigned int>(thisStateNextN, mhitNx);

          if  (!statesNbIndex[MHITINDEX(indexNx)]) {

#ifdef USEQUEUEMHIT
            // compute level
            if (stateN > stateN_of_level_i) {
              stateN_of_level_i = (result->size() - 1);
              level_i++;
            }

            if (level_i <= depth ){
#endif
              // create a new state
              stateNx = result->addNewState();
              statesNbIndex[MHITINDEX(indexNx)] = stateNx;
              statesNbRemaining.push(indexNx);

#ifndef NOMATRIX
              if ( ///@todo{FIXME : not sure that this is correct ... need to be checked twice}
                  gv_subalignment_flag  && (this->_init_states.size() > 0)
                  &&
                  (thisStateNextN == this->_init_states[result->_init_states.size()%(this->_init_states.size())])
                  &&
                  (mhitNx == 0)
                   )
                {
                  result->_init_states.push_back(stateNx);
                }
#endif

              VERB_FILTER(VERBOSITY_DEBUGGING, INFO__("$push state:" << dec << stateNx););

#ifdef USEQUEUEMHIT
            } else {
              // max level reached : goes to a "non final" loop state
              stateNx = result->addNewState();
              result->selfLoop(stateNx);
            }
#endif
          } else {
            stateNx = statesNbIndex[MHITINDEX(indexNx)];
          }
        }

        VERB_FILTER(VERBOSITY_DEBUGGING, INFO__("> add transition ( a:" << dec << a << ", q1:" << dec << stateN << ", q2:" << stateNx << " ) "););

        // add a transition on from stateN --a--> stateNx.
        result->addNewTransition(a,stateN,stateNx,iter->_prob);
      }// for (list)
    }// for (a)
  }//while stack nonempty

  // Free unused data needed to build the automaton
  statesNbIndex.clear();
  return result;
}




/*
 * This method builds a simple Bernoulli
 * probabilistic automaton given a probability
 * vector for each alignment alphabet letter.
 */


int automaton::Automaton_Bernoulli(const vector <double> & p /*[gv_align_alphabet_size]*/) {
#ifdef ASSERTB
  assert((int)p.size() == gv_align_alphabet_size);
#endif
  addNewState();
  int InitState_I = addNewState();
  for (int a = 0; a < gv_align_alphabet_size; a++) {
    addNewTransition(a, InitState_I,  InitState_I,p[a]);
  }
  return 0;
}


/*
 * This methodes builds a Markov
 * probabilistic automaton given a probability
 * vector for each word of size k on the alignment alphabet.
 */

int automaton::Automaton_Markov(const vector<double> & p, const int k ){

  int ApowK     = 1;
  int ApowKplus = gv_align_alphabet_size;


  // [A] build the trie
  int nextbase  = 1;

  // empty ending state
  addNewState(TRUE);
  selfLoop(0);

  int InitState_I = addNewState();    // starting at state 1


  for (int d = 1; d <= k; d++, ApowK *= gv_align_alphabet_size, ApowKplus *= gv_align_alphabet_size) {

    nextbase += ApowK;

    // [a] compute weights
    // (a.1) denumerators
    vector<double> probasum  = vector<double>(ApowK,    0.0);
    int i = 0;
    for (int buklet = 0; buklet < ApowK; buklet++) {
      for (int ibuklet = 0; ibuklet < (int)p.size(); ibuklet += ApowK) {
        probasum[buklet] += p[i++];
      }
    }
    // (a.2) numerators
    vector<double> proba     = vector<double>(ApowKplus,0.0);
    int j = 0;
    for (int buklet = 0; buklet < ApowKplus; buklet++) {
      for (int ibuklet = 0; ibuklet < (int)p.size(); ibuklet += ApowKplus) {
        proba[buklet]    += p[j++];
      }
    }

    // [b] set states an transitions
    for (int buklet = 0; buklet < ApowKplus; buklet++) {
      int state = addNewState();
      int base  = state -  InitState_I - 1;
      addNewTransition(base%gv_align_alphabet_size,
                       base/gv_align_alphabet_size +  InitState_I,
                       state,
                       proba[buklet]/probasum[buklet/gv_align_alphabet_size]);
    } // buklet
  } // depth


  // [B] add final transitions
  for (int b = 0; b < (int)(ApowK); b ++) {
    int bplus = b * gv_align_alphabet_size;
    double psum = 0.0;
    for (int a = 0; a < gv_align_alphabet_size; a++) {
      psum += p[bplus + a];
    }

    for (int a = 0; a < gv_align_alphabet_size; a++) {
      addNewTransition(a ,
                       nextbase +  (b        )%ApowK, // previous word
                       nextbase +  (bplus + a)%ApowK, // next word
                       p[bplus + a]/psum);
    }
  }

  return 0;
}



/*
 * this method builds an Homogeneous alignment model automaton given
 *
 * state 0 : accepting state
 * state 1 : starting state
 * state 2 : rejecting state
 *
 *
 *   /\/\...
 *  /\/\/...
 *
 */

int automaton::Automaton_Homogeneous(const vector<int> & scores,
                                     const int length) {
  bool AcceptingStateReachable = false;
  int  Maxscore = 0;

  for (int i = 0; i < gv_align_alphabet_size; i++) {
    Maxscore = MAX(scores[i],Maxscore);
  }

  if (Maxscore == 0) {
    _ERROR("Automaton_Homogeneous","not positive score is provided");
  }

  // set the vector statesOfScores[2][currentscore][lastmaxscore]
  vector< vector< vector<int> > > statesOfScore = vector< vector< vector <int> > > ( 2, vector< vector <int>  > (Maxscore*(length-1)+1,vector<int>(0)));
  for (int currentscore = 0;  currentscore < Maxscore*(length-1)+1;  currentscore++) {
    statesOfScore[0][currentscore] = vector<int>(Maxscore*(length-1)+1-currentscore,0);
    statesOfScore[1][currentscore] = vector<int>(Maxscore*(length-1)+1-currentscore,0);
  }

  int stateFinalAccepting        = addNewState(TRUE);
  selfLoop(stateFinalAccepting);
  int stateInit                  = addNewState();
  int stateReject                = addNewState();
  selfLoop(stateReject);

  // insert the first state in the table
  statesOfScore[0][0][0] = stateInit;

  for (int l = 0; l < length; l++ ) {
    vector< vector<int> > & statesFrom =  statesOfScore[l%2];
    vector< vector<int> > & statesTo   =  statesOfScore[(l+1)%2];
    for (int currentscore = 0; currentscore <= Maxscore*l; currentscore++ ) {
      for (int lastreachedscore = 0; lastreachedscore <= Maxscore*l - currentscore; lastreachedscore++ ) {
        int stateFrom = statesFrom[currentscore][lastreachedscore];
        if ( stateFrom > 0 ) {
          for (int a = 0; a < gv_align_alphabet_size; a++) {
            int scoreTo = currentscore+scores[a];
            if (scoreTo < 0) {
              addNewTransition(a,stateFrom,stateReject);
            } else {
              if (l == (length - 1)) {
                if (scores[a] >= lastreachedscore) {
                  addNewTransition(a,stateFrom,stateFinalAccepting);
                  AcceptingStateReachable = true;
                } else {
                  addNewTransition(a,stateFrom,stateReject);
                }
              } else {
                int lastreachedscoreTo = MAX(lastreachedscore-scores[a],0);
                int stateTo = statesTo[scoreTo][lastreachedscoreTo];
                if (stateTo == 0) {
                  stateTo = statesTo[scoreTo][lastreachedscoreTo] = addNewState();
                }
                addNewTransition(a,stateFrom,stateTo);
              }
            }
          }
        }
        statesFrom[currentscore][lastreachedscore] = 0;
      }
    }
  }
  statesOfScore.clear();
  if (!AcceptingStateReachable) {
    _ERROR("Automaton_Homogeneous","cannot reach the given score according to -u and -x parameters");
  }
  return 0;
}


/*
 * Build a counting automaton of alignments symbols >= min_a  (according to the final state "0" that has a "1" as its "final" field)
 */
int automaton::Automaton_CountAlphabetSymbols(const int min_a) {
  addNewState(1);
  addNewState(0);
  for (int i = 0; i < 2; i++)
    for (int a = 0; a < gv_align_alphabet_size; a++)
      addNewTransition(a, i, (a >=  min_a) ? 0 : 1);
  return 0;
}


/*
 * Build an automaton that accepts only sequences with > max_cost
 */


int automaton::Automaton_Lossless(const vector<int> & costs,
                                  const int max_cost) {
  vector<int> stateOfCost = vector<int> (max_cost+1,0); /* from 0 to max_cost */
  queue<int> costNbRemaining;

  // state 0 is final and accepting
  int stateReject = addNewState(TRUE);
  selfLoop(stateReject);
  // state 1 is init state (and of cost 0)
  stateOfCost[0] = addNewState();
  costNbRemaining.push(0);
  while (!costNbRemaining.empty()) {
    int c_from = costNbRemaining.front();
    int s_from = stateOfCost[c_from];
    costNbRemaining.pop();
    for (int a = 0; a < gv_align_alphabet_size; a++) {
      int s_to = stateReject;
      int c_to = c_from + costs[a];
      if (c_to <= max_cost) {
        if (stateOfCost[c_to]) {
          s_to = stateOfCost[c_to];
        } else {
          s_to = addNewState();
          stateOfCost[c_to] = s_to;
          costNbRemaining.push(c_to);
        }
      }
      addNewTransition(a,s_from,s_to);
    }
  }
  return 0;
}



/*
 * This method builds a cyclic automaton
 *
 */


int automaton::Automaton_Cycle(const int cycle,
                               const int * const final_list,
                               const int final_nb) {

  int final_i = 0;

  int Finalstate_I = addNewState(TRUE);
  selfLoop(Finalstate_I);

  // create a sorted list of final states (do not modify the current one used by "seed")
  int * final_list_sorted = new int[final_nb];
  for (int i = 0; i < final_nb; i++)
    final_list_sorted[i] = final_list[i];
  sort(final_list_sorted, final_list_sorted + final_nb);

  // create states
  for (int i = 0; i < cycle; i++) {
    int u;
    if  (final_i < (int)final_nb && final_list_sorted[final_i] == i) {
      u = addNewState(TRUE); // final
      final_i++;
    } else {
      u = addNewState(FALSE); // non final
    }

    // set init_states when subalignment is activated
#ifndef NOMATRIX
    if (gv_subalignment_flag) {
      _init_states.push_back(u);
    }
#endif
  }

  // link states
  for (int i = 0; i < cycle; i++) {
    for (int a = 0; a < gv_align_alphabet_size; a++) {
      addNewTransition(
                       a,
                       1 + (i)%cycle,
                       1 + (i+1)%cycle
                       );
    }
  }

  delete[] final_list_sorted;

  return 0;
}



/*
 * This method generates a new minimized automaton
 * using Hopcroft minimization
 *
 */

typedef pair<int,int> ABpair;
#define REVERSESTATE(c,a,b) (((c) == (a))?(b):(((c) == (b))?(a):(c)))

automaton * automaton::Hopcroft() const {

  VERB_FILTER(VERBOSITY_DEBUGGING, INFO__("== Hopcroft == "););

  int nbStates = this->size();

  vector<int>                   stclass(nbStates);  // class for a given state
  vector<list<int> >            block(nbStates);    // list of states for a given class
  vector<list<int>::iterator >  location(nbStates); // pointer to each stat inside the given list
  vector<int>                   card(nbStates,0);   // cardinality of a given class
  vector<vector<int> >          SET(nbStates,vector<int>(gv_align_alphabet_size,0)); // stack membership
  vector<vector<vector<int> > > PREV(nbStates,vector< vector<int> >(gv_align_alphabet_size, vector<int>(0))); // previous state(s) on 'a'

  int nbClasses = 0;

  // (0) some preprocessing on final labels to get "final_to_class" function
  int max_final = 0;
  int min_final = INT_INFINITY;
  for (int i = 0; i < nbStates; i++) {
    int f = _states[i]._final;
    if (f) {
      min_final = MIN(min_final,f);
      max_final = MAX(max_final,f);
    }
  }

  vector< int > final_to_class(max_final+1,-1);
  {
    vector< int >  count(max_final+1,0); // cardinality of a given final label
    for (int i = 0; i < nbStates; i++)
      count[_states[i]._final]++;
    for (int f = 0; f <= max_final; f++)
      if (count[f])
        final_to_class[f] = nbClasses++;
    // invert "0" and "min_final"
    // (to keep some states where "final = 0" in the class 1 ->
    // -> see why after when rebuiding the minimized automaton [*] )
    final_to_class[0] = min_final;
    final_to_class[1] = 0;
  }


  // (1) initialization
  for (int i = 0; i < nbStates; i++) {
    int _class_ = final_to_class[_states[i]._final];
    if (!card[_class_])
      block[_class_] = list<int>();
    stclass[i] = _class_;
    block[_class_].push_front(i);
    location[i] = block[_class_].begin();
    card[_class_]++;
    // PREV
    for (int a = 0; a < gv_align_alphabet_size; a++)
      for (vector<transition>::const_iterator
             iter  = _states[i]._next[a].begin();
           iter != _states[i]._next[a].end();
           iter ++ ) {
        PREV[iter->_state][a].push_back(i);
      }
  }

  // initial stack
  stack<ABpair> Set;
  // find the max set and put the other(s) in the refining Set
  {
    int max = 0;
    int max_class = 0;
    for (int _class = 0; _class < nbClasses; _class++) {
      if (max < card[_class]) {
        max       = card[_class];
        max_class = _class;
      }
    }
    for (int _class = 0; _class < nbClasses; _class++) {
      if (_class != max_class) {
        for (int a = 0; a < gv_align_alphabet_size; a++) {
          Set.push(ABpair(a,_class));
          SET[_class][a] = 1;
        }
      }
    }
  }

  // used to count (P_a_inv)Intersert(B), states membership to P_a_inv, and B membership to intersection
  vector<int> card_P_a_inv_Intersect_B(nbStates,0);
  vector<int> Bsplit(nbStates,0);

  // (2) main hopcroft loop
  while (!Set.empty()) {
    //(P,a) <- First(S)
    ABpair aP = Set.top(); Set.pop();
    int a = aP.first;
    int P = aP.second;

    SET[P][a] = 0;

    // (2.a) IsRefined precomputation step
    list<int> & Plist               = block[P];
    list<int>   Plist_a_inv         = list<int>();
    // for each state in (P)
    for (list<int>::iterator iter_Plist  = Plist.begin();
         iter_Plist != Plist.end();
         iter_Plist ++ ) {

      int Pstate = *iter_Plist;

      // for each predecessor *--a-->(P)
      for (vector<int>::iterator iter_Pstate_list_a_inv  = PREV[Pstate][a].begin();
           iter_Pstate_list_a_inv != PREV[Pstate][a].end();
           iter_Pstate_list_a_inv ++ ) {

        int Pstate_a_inv = *iter_Pstate_list_a_inv;
        int b = stclass[Pstate_a_inv];
        Plist_a_inv.push_back(Pstate_a_inv);
        card_P_a_inv_Intersect_B[b]++;
      }// transitions
    }// iter_Plist



    int oldNbClasses  = nbClasses;
    vector<int> Bprev = vector<int>();

    // (2.b) foreach a-1(P) state
    for (list<int>::iterator i  = Plist_a_inv.begin();
         i != Plist_a_inv.end();
         i++){
      int Pstate_a_inv = *i;
      int b            = stclass[Pstate_a_inv];
      // if b is refined by a-1(P)

      if ( card_P_a_inv_Intersect_B[b] > 0 && card[b] - card_P_a_inv_Intersect_B[b] > 0) {


        int bp;
        if ( Bsplit[b] != 0 ) {
          // a previous state has been already moved from b to b'=Bsplit[b]
          bp = Bsplit[b];
        } else {
          // first state to be splited from b
          block[nbClasses]  = list<int>();
          bp = Bsplit[b]    = nbClasses;
          Bprev.push_back(b);
          nbClasses++;
        }

        block[b].erase(location[Pstate_a_inv]);
        block[bp].push_front(Pstate_a_inv);
        location[Pstate_a_inv] = block[bp].begin();
        card[bp]++;
        stclass[Pstate_a_inv] = bp;

      } else {

        card_P_a_inv_Intersect_B[b] = 0;

      }
    } // for each a-1(P)
    Plist_a_inv.clear();

    // (2.c) foreach Bp state
    // Update(S,B,B',B'')
    int i = 0;
    for (int bp = oldNbClasses; bp < nbClasses; bp++) {
      int b = Bprev[i++];
      // reset tables
      Bsplit[b] = 0;
      card[b]  -= card[bp];
      card_P_a_inv_Intersect_B[b] = 0;

      for (int a = 0; a < gv_align_alphabet_size; a++) {

        if (SET[b][a]) {
          SET[bp][a] = 1;
          Set.push(ABpair(a,bp));
        } else {
          if (card[b] > card[bp]) {
            SET[bp][a] = 1;
            Set.push(ABpair(a,bp));
          } else {
            SET[b][a] = 1;
            Set.push(ABpair(a,b));
          }
        }
      }
    }
    Bprev.clear();

  }//while (!empty)


  // Clear Some Hopcroft data
  Bsplit.clear();
  card_P_a_inv_Intersect_B.clear();
  location.clear();
  card.clear();
  SET.clear();


  // create the minimized automaton

  // [*] see why class "1" is "non-final" before
  automaton * result = new automaton();
  for (int c = 0; c < nbClasses; c++) {
    result->addNewState(_states[*(block[c].begin())]._final);

    if (block[c].size() > 1) {
       VERB_FILTER(VERBOSITY_DEBUGGING, INFO__(
          "merging : ";
          for(list<int>::iterator i = block[c].begin(); i != block[c].end(); i++) {
             cerr << (*i) << " ";
          }
       ););
    }
  }

  // put the transitions with the first state (class "stclass[1]") on pos "1"
  for (int c = 0; c < nbClasses; c++) {
    for (int a = 0; a < gv_align_alphabet_size; a++) {
#ifdef ASSERTB
      assert(DETERMINISTIC(*(block[c].begin()),a));
#endif
      result->addNewTransition(a,
                               REVERSESTATE(c,stclass[1],1),
                               REVERSESTATE(stclass[(_states[*(block[c].begin())]._next[a].begin()->_state)],stclass[1],1)
                               );
    }
  }

#ifndef NOMATRIX
  if (gv_subalignment_flag && (_init_states.size() > 0)) {
    for (int i = 0;i<(int)_init_states.size();i++)
      result->_init_states.push_back(REVERSESTATE(stclass[_init_states[i]],stclass[1],1));
  }
#endif

  stclass.clear();
  block.clear();
  return result;
}


// check isomorphism

bool automaton::isIsomorphTo(const automaton & other) const {


  // differents sizes
  if (this->size() != other.size()) {
    return false;
  }


  // check the mapping
  vector<int> map_this_to_other = vector<int>(this->size(),-1);
  stack< pair<int,int> >  statesNbRemaining;

  statesNbRemaining.push(pair<int,int>(1,1));

  while (!statesNbRemaining.empty()) {

    pair<int,int> states = statesNbRemaining.top();
    statesNbRemaining.pop();

    if (map_this_to_other[states.first]<0) {

      // not already found before
      map_this_to_other[states.first] = states.second;

      for (int a = 0; a < gv_align_alphabet_size; a++) {
        // push neighboors
        pair<int,int> states_next = pair<int,int>(this->_states[states.first]._next[a].begin()->_state,other._states[states.second]._next[a].begin()->_state);
        statesNbRemaining.push(states_next);
      }

    } else {
      // check otherwise the consistency
      if ( map_this_to_other[states.first] != states.second ) {
        return false;
      }
    }
  } // while not empty
  return true;
}

// build a very simple linear automaton from the seed with a reject bag states

int automaton::Automaton_SeedLinearMatching(const seed& s,
                                            const vector< vector <int> > & matchingmatrix) {

  // Get the seed span
  int   motif_span = s.span();
  int * motif      = s.table();

  // create a first state [0] and put it as the final one
  int Finalstate_I = addNewState(TRUE);
  selfLoop(Finalstate_I);

  // create a second state [1]  : it will be the initial one
  int Initstate_I = addNewState();

  // bag reject set [2]
  int RejectBagstate_I  = addNewState();
  selfLoop(RejectBagstate_I);

  // build a linear automaton with fail as a reject state ...
  int Prevstate_I = Initstate_I;
  for (int i = 0; i < motif_span; i++) {
    int b = motif[i];

    int Nextstate_I = Finalstate_I;
    if (i < motif_span - 1)
      Nextstate_I = addNewState();

    for (int a = 0; a < gv_align_alphabet_size; a++){
      if (MATCHES_AB(a,b)){
        addNewTransition(a,Prevstate_I,Nextstate_I);
      } else {
        addNewTransition(a,Prevstate_I,RejectBagstate_I);
      }
    }
    Prevstate_I = Nextstate_I;
  }

  return 0;
}

/** @class SeedPrefixesMatchingSet_old
 *  @brief simple "inner class" to keep states @f$ <X,t> @f$ and @f$ k @f$, 
 *         during SeedPrefixMatching_old method
 *  @see Automaton_SeedPrefixesMatching_old
 */

class SeedPrefixesMatchingSet_old {

public:
  /// set of prefixes matching
  int X;
  /// lenght of the last run of '1'
  short t;
  /// length of the maximal prefix in the set @f$ X @f$ : @f$ k = max(X) @f$
  short k;
  /** @brief build a SeedPrefixesMatchingSet_old
   */
  SeedPrefixesMatchingSet_old(int X = 0, short t = 0, short k = 0) : X(X), t(t), k(k) {};
};

#define SEEDPREFIXESMATCHING_OLD_X(s)          (statesSeedPrefixesMatchingSet_old[(s)].X)
#define SEEDPREFIXESMATCHING_OLD_T(s)          (statesSeedPrefixesMatchingSet_old[(s)].t)
#define SEEDPREFIXESMATCHING_OLD_K(s)          (statesSeedPrefixesMatchingSet_old[(s)].k)

int automaton::Automaton_SeedPrefixesMatching_old(const seed& s,
                                                  const vector< vector <int> > & matchingmatrix,
                                                  const bool nomerge) {

  VERB_FILTER(VERBOSITY_DEBUGGING, INFO__("== SeedPrefixesMatching_old == (nomerge:" << dec << nomerge << ")"););

  // Get the seed span
  int   motif_span = s.span();
  int * motif      = s.table();

#ifdef ASSERTB
  if (motif_span <= 0) {
    _ERROR("SeedPrefixesMatching_old","null or negative span");
  }

  for (int i = 0; i < motif_span; i++)
    if ( motif[i] < 0 || motif[i] >= gv_seed_alphabet_size)
      _ERROR("SeedPrefixesMatching_old","incorrect seed element");
#endif

VERB_FILTER(VERBOSITY_DEBUGGING, INFO__(
  " motif_span :" << motif_span << ", motif:";
   for (int i = 0; i < motif_span; i++)
     cerr << motif[i] << " ";
   cerr << endl;
););

  // Compute xset_bitsize and X
  vector<int> L(motif_span+1,0);
  int xset_bitsize  = 0; // bitsize of the set
  int xset_size     = 0; // size of the set
  L[0] = -1;
  for (int i = 0; i<motif_span; i++) {
    if (motif[i] < (gv_seed_alphabet_size - (gv_matching_symbol_flag?1:0))) {
      xset_bitsize++;
      L[xset_bitsize] = i;
      VERB_FILTER(VERBOSITY_DEBUGGING, INFO__(" L[" << xset_bitsize << "] : " <<  L[xset_bitsize] << endl;););
    }
  }
  xset_size = 1 << xset_bitsize;

  VERB_FILTER(VERBOSITY_DEBUGGING, INFO__(" xset_size : " << xset_size << ",\t" << " xset_bitsize : " << xset_bitsize << endl;););

  // Compute TCODE[t] to get an index by use of IndexNb = TCODE[t] + X coding system
  vector<int> TCODE(motif_span+1+(nomerge?1:0),0);
  TCODE[0] = 0;
  for (int maxxsetsize = xset_size, l = xset_bitsize, t = 1; t <= motif_span+(nomerge?1:0); t++) {
    if (L[l] == motif_span+(nomerge?1:0)-t) {
      maxxsetsize >>= 1;
      l--;
    }
    TCODE[t] = TCODE[t-1] + maxxsetsize;
    VERB_FILTER(VERBOSITY_DEBUGGING, INFO__(" TCODE[" << t << "] : " << TCODE[t] << endl;););
  }

  // 1) Precompute E,F functions
  vector< vector<int> > FX(motif_span+(nomerge?1:0),vector<int>(gv_align_alphabet_size,0));
  vector< vector<int> > Fk(motif_span+(nomerge?1:0),vector<int>(gv_align_alphabet_size,0));
  vector< vector< vector<int> > > EX(motif_span+(nomerge?1:0), vector< vector<int> >((xset_bitsize+1), vector<int>(gv_align_alphabet_size,0)));
  vector< vector< vector<int> > > Ek(motif_span+(nomerge?1:0), vector< vector<int> >((xset_bitsize+1), vector<int>(gv_align_alphabet_size,0)));

  // F[t][a]
  for (int t = 0; t<motif_span+(nomerge?1:0); t++) {
    for (int i = 1; i<= xset_bitsize && L[i] <= t; i++) {
      int r = 1 << (i-1);
      int b = motif[L[i]];

      for (int a = 0; a < gv_align_alphabet_size; a++){
        if (MATCHES_AB(a,b)){
          FX[t][a] = (t>0) ? (FX[t-1][a]) | r : r;
          Fk[t][a] = i;
        }
      }
    }
  }


  // E[t][a]
  for (int i = 1; i<=xset_bitsize; i++) {
    int L_i = L[i];
    int b   = motif[L_i];

    for (int a = 0; a < gv_align_alphabet_size; a++) {
      // does b match the mismatch a
      if (MATCHES_AB(a,b)) {
        // if X[ib]-t-1 is a joker position, do a mask to have its bit position.
        for (int k = 1; k < i; k++) {
          int t = L[i] - L[k] - 1;
#ifdef ASSERTB
          assert(t>=0);
          assert(t<motif_span);
          assert(k>0);
          assert(k<=xset_bitsize);
#endif
          EX[t][k][a] = 1<<(i-1);
          Ek[t][k][a] = i;
        }
      }
    }
  }


  // 2) automaton build

  // fast state index (designed to retrieve a state given <X,t> code
  vector<int> statesNbIndex(TCODE[motif_span+(nomerge?1:0)]+1, 0);
  // queue/stack used to store non preprocessed states <X,t> code
#ifdef USEQUEUEAUTOMATON
  queue<int>  statesNbRemaining;
#else
  stack<int>  statesNbRemaining;
#endif
  // keep each state information (X,t,k) inside this table
  vector<SeedPrefixesMatchingSet_old> statesSeedPrefixesMatchingSet_old(0);

  // create a first state [0] and put it as the final one
  int Finalstate_I = addNewState(TRUE);
  statesSeedPrefixesMatchingSet_old.push_back(SeedPrefixesMatchingSet_old(0,0,0));
  selfLoop(Finalstate_I);

  // create a second state [1]  : it will be the initial one
  int Initstate_I   = addNewState();
  statesSeedPrefixesMatchingSet_old.push_back(SeedPrefixesMatchingSet_old(0,0,0));
  statesNbIndex[0]  = Initstate_I;
  statesNbRemaining.push(Initstate_I);


  while (!statesNbRemaining.empty()) {

    // current state remaining
#ifdef USEQUEUEAUTOMATON
    int Xstate_I = statesNbRemaining.front();
#else
    int Xstate_I = statesNbRemaining.top();
#endif
    statesNbRemaining.pop();

#ifdef ASSERTB
    assert(Xstate_I >= 0);
    assert(Xstate_I <  (int)_states.size());
#endif
    int Xstate_X = SEEDPREFIXESMATCHING_OLD_X(Xstate_I);
    int Xstate_T = SEEDPREFIXESMATCHING_OLD_T(Xstate_I);
    int Xstate_K = SEEDPREFIXESMATCHING_OLD_K(Xstate_I);

    VERB_FILTER(VERBOSITY_DEBUGGING, INFO__("$pop  I:" << Xstate_I << " < X:" << Xstate_X << ", t:" << Xstate_T <<", k:"<< Xstate_K <<" >"));

#ifdef ASSERTB
    assert(Xstate_X >= 0);
    assert(Xstate_X < xset_size);
    assert(Xstate_T >= 0);
    assert(Xstate_T < motif_span + (nomerge?1:0));
    assert(Xstate_K >= 0);
    assert(Xstate_K <= xset_bitsize);
    assert(TCODE[Xstate_T] + Xstate_X < TCODE[motif_span+(nomerge?1:0)]);
#endif



    int XPstate_Ix = 0, XPstate_I = 0, XPstate_X = 0;

    if (Xstate_K > 0) {
      XPstate_X   = Xstate_X ^ (1 << (Xstate_K-1));    // get X'
      XPstate_Ix  = TCODE[Xstate_T] + XPstate_X;   // get (X',t) index
      XPstate_I   = statesNbIndex[XPstate_Ix];
    }

    for (int a = 0; a < gv_align_alphabet_size; a++) {

      VERB_FILTER(VERBOSITY_DEBUGGING, INFO__("a = " << a););

      // next state according to letter "a"
      int Ystate_X = Xstate_X;
      int Ystate_T = Xstate_T;
      int Ystate_K = Xstate_K;

      if ( a == (gv_align_alphabet_size - (gv_matching_symbol_flag?1:0))) {
        Ystate_T++;
      } else {

        if ( Xstate_K > 0 ) {
#ifdef ASSERTB
          assert(DETERMINISTIC(XPstate_I,a));
#endif
          int YPstate_I = _states[XPstate_I]._next[a].begin()->_state;
          // YP result
          int YPstate_X = SEEDPREFIXESMATCHING_OLD_X(YPstate_I);
          int YPstate_K = SEEDPREFIXESMATCHING_OLD_K(YPstate_I);
          //  Y result
#ifdef ASSERTB
          assert(Xstate_T >= 0);
          assert(Xstate_T < motif_span+(nomerge?1:0));
          assert(Xstate_K > 0);
          assert(Xstate_K <= xset_bitsize);
#endif

          Ystate_X  =      YPstate_X | EX[Xstate_T][Xstate_K][a];
          Ystate_K  = MAX( YPstate_K , Ek[Xstate_T][Xstate_K][a] );
#ifdef BUILD
          cout << "**V "<< endl;
          cout << "YPstate_X:" <<  YPstate_X << ",  YPstate_K:" << YPstate_K << endl;
          cout << "Ystate_X:"  <<  Ystate_X  << ",  Ystate_K:"  << Ystate_K << endl;

#endif
        } else {
          Ystate_X  = FX[Xstate_T][a];
          Ystate_K  = Fk[Xstate_T][a];

#ifdef BUILD
          cout << "*U "<< endl;
          cout << "Ystate_X:"  <<  Ystate_X << ",  Ystate_K:" << Ystate_K<< endl;

#endif

        } // if ( Xstate_K > 0 )
        Ystate_T  = 0;

      }// if (a == gv_alphabet_size - 1)

      // Y result Index
      int Ystate_Ix = TCODE[Ystate_T] + Ystate_X;
      int Ystate_I  = statesNbIndex[Ystate_Ix];
      int final     = ((Ystate_T + L[Ystate_K]) >= (motif_span - 1)) ? TRUE : FALSE;


      //Y : either final state (->0) or "after final state" thus a previous one
      if (Ystate_T + L[Ystate_K] >= motif_span - 1 + (nomerge?1:0)) {
        if (nomerge){ // final states are not merged and  their transitions are consided as "normal states" :
          // the only difference is their maximal prefix matching that cannot be extended more that "span"
          if (Ystate_K > 0)
            Ystate_X ^= 1<<(Ystate_K-1);
          else
            Ystate_T--;
          Ystate_Ix = TCODE[Ystate_T] + Ystate_X;
          Ystate_I  = statesNbIndex[Ystate_Ix];
        }else{ // final states merged and absorbant (->0)
          Ystate_I = 0;
        }
      } else {
        // else check if the state has not been already created
        if (!Ystate_I) {
          // create it and push it inside the list of states which need update.
          Ystate_I = addNewState(final);
          statesSeedPrefixesMatchingSet_old.push_back(SeedPrefixesMatchingSet_old(Ystate_X,Ystate_T,Ystate_K));
          statesNbIndex[Ystate_Ix] =  Ystate_I;
          statesNbRemaining.push(Ystate_I);

          VERB_FILTER(VERBOSITY_DEBUGGING, INFO__("$push I:"<< Ystate_I << " < X:" << hex << Ystate_X << ", t:" << dec << Ystate_T << ", k:" << Ystate_K << " >"););

        }
      }

      VERB_FILTER(VERBOSITY_DEBUGGING, INFO__("> add transition ( a:" << dec << a << ", q1:" << dec << Xstate_I << ", q2:" << Ystate_I << " )"););

      addNewTransition(a,Xstate_I,Ystate_I);
#ifdef ASSERTB
      assert(DETERMINISTIC(Xstate_I,a));
#endif
    }// for all letter in ALPHABET.
  }// while states remain.

  // Free unused data needed to build the automaton
  L.clear();
  FX.clear();
  Fk.clear();
  EX.clear();
  Ek.clear();
  statesNbIndex.clear();
  statesSeedPrefixesMatchingSet_old.clear();
  return 0;
}






/** @class SeedPrefixesMatchingSet
 *  @brief simple "inner class" to keep states @f$ <X,t> @f$ and @f$ k @f$, together with @f$ Xp @f$ to keep previous track, 
 *         during SeedPrefixMatching method
 *  @see Automaton_SeedPrefixesMatching
 */

class SeedPrefixesMatchingSet {

public:
  /// set of prefixes matching
  int X;
  /// lenght of the last run of '1'
  short t;
  /// length of the maximal prefix in the set @f$ X @f$ : @f$ k = max(X) @f$
  short k;
  /// gives the state @f$< X' = X / max(X), t > @f$ (this state does exists and is created before @f$ <X,t> @f$)
  int Xp;
  /** @brief gives the last state @f$ <X'',t> @f$ created that verifies :
   *   @li @f$ X''/max(X'') = X @f$
   *   @li @f$ max(X'') @f$ is greatest among all created @f$ X'' @f$ verifying the previous condition
   */
  int RevMaxXp;
  /** @brief build a SeedPrefixesMatchingSet
   */
  SeedPrefixesMatchingSet(int X = 0,short t = 0,short k = 0,int Xp = 0,int RevMaxXp = 0) : X(X), t(t), k(k), Xp(Xp), RevMaxXp(RevMaxXp) {};
};

#define SEEDPREFIXESMATCHING_X(s)          (statesSeedPrefixesMatchingSet[(s)].X)
#define SEEDPREFIXESMATCHING_T(s)          (statesSeedPrefixesMatchingSet[(s)].t)
#define SEEDPREFIXESMATCHING_K(s)          (statesSeedPrefixesMatchingSet[(s)].k)
#define SEEDPREFIXESMATCHING_XP(s)         (statesSeedPrefixesMatchingSet[(s)].Xp)
#define SEEDPREFIXESMATCHING_REVMAXXP(s)   (statesSeedPrefixesMatchingSet[(s)].RevMaxXp)


int automaton::Automaton_SeedPrefixesMatching(const seed& s,
                                              const vector< vector <int> > & matchingmatrix,
                                              const bool nomerge) {

  VERB_FILTER(VERBOSITY_DEBUGGING, INFO__("== SeedPrefixesMatching == (nomerge:" << dec << nomerge << ")"););

  // Get the seed span
  int   motif_span = s.span();
  int * motif      = s.table();

#ifdef ASSERTB
  if (motif_span <= 0) {
    _ERROR("SeedPrefixesMatching","null or negative span");
  }

  for (int i = 0; i < motif_span; i++)
    if ( motif[i] < 0 || motif[i] >= gv_seed_alphabet_size)
      _ERROR("SeedPrefixesMatching","incorrect seed element");
#endif

#ifdef BUILD
  cout << " motif_span :" << motif_span << ", motif:";
  for (int i = 0; i < motif_span; i++ )
    cout << motif[i] << " ";
  cout << endl;
#endif


  // Compute xset_bitsize and X
  vector<int> L(motif_span+1,0);
  int xset_bitsize = 0; // bitsize of the set
  L[0] = -1;
  for (int i = 0; i < motif_span; i++) {
    if (motif[i] < (gv_seed_alphabet_size - (gv_matching_symbol_flag?1:0))) {
      xset_bitsize++;
      L[xset_bitsize] = i;
#ifdef BUILD
      cout << " L[" << xset_bitsize << "] : " <<  L[xset_bitsize] << endl;
#endif
    }
  }

#ifdef BUILD
  {
    int xset_size = 1 << xset_bitsize;
    cout << " xset_size :"     << xset_size << ",\t" <<
            " xset_bitsize : " << xset_bitsize << endl;
  }
#endif


  // 1) Precompute E,F functions
  vector< vector<int> > FX(motif_span+(nomerge?1:0),vector<int>(gv_align_alphabet_size,0));
  vector< vector<int> > Fk(motif_span+(nomerge?1:0),vector<int>(gv_align_alphabet_size,0));
  vector< vector< vector<int> > > EX(motif_span+(nomerge?1:0), vector< vector<int> >((xset_bitsize+1), vector<int>(gv_align_alphabet_size,0)));
  vector< vector< vector<int> > > Ek(motif_span+(nomerge?1:0), vector< vector<int> >((xset_bitsize+1), vector<int>(gv_align_alphabet_size,0)));

  // F[t][a]
  for (int t = 0; t < motif_span+(nomerge?1:0); t++) {
    for (int i = 1; i<= xset_bitsize && L[i] <= t; i++) {
      int r = 1<<(i-1);
      int b = motif[L[i]];

      for (int a = 0; a < gv_align_alphabet_size; a++){
        if (MATCHES_AB(a,b)){
          FX[t][a] = (t>0) ? (FX[t-1][a]) | r : r;
          Fk[t][a] = i;
        }
      }
    }
  }


  // E[t][a]
  for (int i = 1; i<=xset_bitsize; i++) {
    int L_i = L[i];
    int b   = motif[L_i];

    for (int a = 0; a < gv_align_alphabet_size; a++) {
      // does b match the mismatch a
      if (MATCHES_AB(a,b)) {
        // if X[ib]-t-1 is a joker position, do a mask to have its bit position.
        for (int k = 1; k < i; k++) {
          int t = L[i] - L[k] - 1;
#ifdef ASSERTB
          assert(t>=0);
          assert(t<motif_span);
          assert(k>0);
          assert(k<=xset_bitsize);
#endif
          EX[t][k][a] = 1<<(i-1);
          Ek[t][k][a] = i;
        }
      }
    }
  }


  // 2) automaton build
  // queue/stack used to store non preprocessed states <X,t> code
  queue<int>  statesNbRemaining;

  // keep each state information (X,t,k) inside this table
  vector<SeedPrefixesMatchingSet> statesSeedPrefixesMatchingSet(0);

  // create a first state [0] and put it as the final one
  int Finalstate_I = addNewState(TRUE);
  statesSeedPrefixesMatchingSet.push_back(SeedPrefixesMatchingSet(0,0,0,0,0));
  selfLoop(Finalstate_I);

  // create a second state [1]  : it will be the initial one
  int Initstate_I   = addNewState();
  statesSeedPrefixesMatchingSet.push_back(SeedPrefixesMatchingSet(0,0,0,1,1));


  // create first level states and push them in the queue
  int State_R = 0;
  for (int a = 0; a < gv_align_alphabet_size; a++) {
    bool matches = MATCHES_AB(a,motif[0]);
    int  final   = ((motif_span == 1) && matches) ? TRUE : FALSE;

    if ( a == (gv_align_alphabet_size - (gv_matching_symbol_flag?1:0))) {
      if (final && !nomerge) {
        addNewTransition(a,Initstate_I,Finalstate_I);
      } else {
        int State_I = addNewState(final);
        statesSeedPrefixesMatchingSet.push_back(SeedPrefixesMatchingSet(0,1,0,Initstate_I,-1));
        statesNbRemaining.push(State_I);
        addNewTransition(a,Initstate_I,State_I);
      }
    } else {
      if (matches) {
        if (final && !nomerge) {
          State_R = Finalstate_I;
        } else {
          if (!State_R) {
            State_R = addNewState(final);
            statesSeedPrefixesMatchingSet.push_back(SeedPrefixesMatchingSet(1,0,1,Initstate_I,-1));
            statesNbRemaining.push(State_R);
          }
        }
        addNewTransition(a,Initstate_I,State_R);
      } else {
        addNewTransition(a,Initstate_I,Initstate_I);
      }
    }
  }






  // main loop
  while (!statesNbRemaining.empty()) {

    // current state Xstate
    int Xstate_I = statesNbRemaining.front();
    statesNbRemaining.pop();

#ifdef ASSERTB
    assert(Xstate_I >= 0);
    assert(Xstate_I <  (int)_states.size());
#endif

    int Xstate_X    = SEEDPREFIXESMATCHING_X (Xstate_I);
    int Xstate_T    = SEEDPREFIXESMATCHING_T (Xstate_I);
    int Xstate_K    = SEEDPREFIXESMATCHING_K (Xstate_I);

#ifdef ASSERTB
    assert(Xstate_X >= 0);
    assert(Xstate_X < (1<<xset_bitsize));
    assert(Xstate_T >= 0);
    assert(Xstate_T < motif_span + (nomerge?1:0));
    assert(Xstate_K >= 0);
    assert(Xstate_K <= xset_bitsize);
#endif

    // XPstate : defined as "Xstate" where the maximal matching prefix fails ...
    int XPstate_I   = SEEDPREFIXESMATCHING_XP(Xstate_I);

#ifdef ASSERTB
    assert(XPstate_I >= 0);
    assert(XPstate_I < (int)_states.size());
#endif

    VERB_FILTER(VERBOSITY_DEBUGGING, INFO__("$pop  I:" << Xstate_I << " < X:" << Xstate_X << ", t:" << Xstate_T <<", k:"<< Xstate_K << " > "););

    // for each letter (a) : compute the transition (Xstate)--a-->Ystate [and determine Ystate]
    for (int a = 0; a < gv_align_alphabet_size; a++) {

      VERB_FILTER(VERBOSITY_DEBUGGING, INFO__("a = " << a););

      // set at true if Ystate has not been already created before
      bool Y_not_there = false;

      // YPstate : defined as (XPstate)--a-->
      int YPstate_I   = _states[XPstate_I]._next[a].begin()->_state;

#ifdef ASSERTB
      assert(YPstate_I >= 0);
      assert(YPstate_I < (int)_states.size());
#endif

      int YPstate_T   = SEEDPREFIXESMATCHING_T(YPstate_I);
      int YPstate_X   = SEEDPREFIXESMATCHING_X(YPstate_I);
      int YPstate_K   = SEEDPREFIXESMATCHING_K(YPstate_I);

#ifdef ASSERTB
      assert(YPstate_X >= 0);
      assert(YPstate_X < (1<<xset_bitsize));
      assert(YPstate_T >= 0);
      assert(YPstate_T < motif_span + (nomerge?1:0));
      assert(YPstate_K >= 0);
      assert(YPstate_K <= xset_bitsize);
#endif

      // YPRstate ...value is either -1 if undefined,
      // otherwise >= 0 to indicate the last state created that
      // verifies YPRstate / maximalprefixmatching{YPRstate} = YPstate.
      int YPRstate_I  = SEEDPREFIXESMATCHING_REVMAXXP(YPstate_I);

#ifdef ASSERTB
      assert(YPRstate_I >= (-1));
      assert(YPRstate_I < (int)_states.size());
#endif

      // Ystate_I : next state according to letter "a"
      // by default values are given here
      // and will be processed
      int Ystate_I    = Xstate_I;
      int Ystate_X    = Xstate_X;
      int Ystate_T    = Xstate_T;
      int Ystate_K    = Xstate_K;

      // (a == match)
      if ( a == (gv_align_alphabet_size - (gv_matching_symbol_flag?1:0))) {
        // UPDATE STATE VALUES
        Ystate_T++;         // (a) add a 1 to the run of 1 and
        // NOMERGE : CHECK FINAL STATE
        if (nomerge && L[Ystate_K] + Ystate_T >= motif_span) {
          if (Ystate_K > 0) {
            Ystate_I = YPstate_I;
            Ystate_T = YPstate_T;
            Ystate_X = YPstate_X;
            Ystate_K = YPstate_K;
          } else {
            Ystate_I = Xstate_I;
            Ystate_T = Xstate_T;
            Ystate_X = Xstate_X;
            Ystate_K = Xstate_K;
          }
        } else {
          //UPDATE STATE INDEX
          Y_not_there = true;
        }
      } else {
        // (a != match)
        // either (K > 0) or (K == 0)
        if ( Xstate_K > 0 ) {

#ifdef ASSERTB
          assert(DETERMINISTIC(XPstate_I,a));
#endif
          // (b) add a possible prefix when compared with XP ??
          //   - either EX[Xstate_T][Xstate_K][a] == 0
          //     and  Y = YP : thus  Y already exists
          //   - otherwise EX = 1<<Ek and thus the possible Y may exist ( it can be YPR )

          if (EX[Xstate_T][Xstate_K][a]) {

#ifdef BUILD
            cerr << "EX != 0 : maximal prefix does match " << endl;
#endif
            // UPDATE STATE VALUES
            Ystate_X  =      YPstate_X | EX[Xstate_T][Xstate_K][a];
            Ystate_K  = MAX( YPstate_K , Ek[Xstate_T][Xstate_K][a]);

            // NOMERGE : CHECK FINAL STATE
            if (nomerge && L[Ystate_K] >= motif_span) {
#ifdef ASSERTB
              assert(YPstate_X ==  (Ystate_X^(1 << (Ystate_K-1))));
#endif
              goto Ex0;
            }

            // UPDATE STATE INDEX
            if (YPRstate_I >=0 && SEEDPREFIXESMATCHING_X(YPRstate_I) == Ystate_X  /* COMMENT:NOT NEEDED && SEEDPREFIXESMATCHING_T(YPRstate_I) == 0 */ ) {
#ifdef ASSERTB
              assert(SEEDPREFIXESMATCHING_T(YPRstate_I) == 0);
#endif
              Ystate_I = YPRstate_I; // state already exists (has been processed before)
            } else {
              Y_not_there = true;     // state must be created for the first time
            }
          } else {
          Ex0:
#ifdef BUILD
            cerr << "Ex == 0 : maximal prefix doesnt match : Y is thus YP" << endl;
#endif

            // UPDATE STATE VALUES
            Ystate_X  =      YPstate_X;
            Ystate_K  =      YPstate_K;
            // NOMERGE : CHECK FINAL STATE
            // ** not needed **
            // UPDATE STATE INDEX
            Ystate_I  =      YPstate_I;
          }

#ifdef BUILD
          cerr << "**V " << endl;
          cerr << "YPstate_X:" <<  YPstate_X << ",  YPstate_K:" << YPstate_K << endl;
          cerr << "Ystate_X:"  <<  Ystate_X  << ",  Ystate_K:"  << Ystate_K  << endl;
#endif
        } else { // Xstate is empty ... and we have read 1^t.a  since ...

          // (c) initial case when reading 1^t . a
          //     three subcases can be deduced here :
          //     - either there is one prefix matching on the full 1^t.a span ? in this case does such state exists ?
          //     - either the maximal prefix does not
          //
#ifdef BUILD
          cerr << " FX on 1^"<< Xstate_T << "." << a << endl;
#endif

          // UPDATE STATE VALUES
          Ystate_X   = FX[Xstate_T][a];
          Ystate_K   = Fk[Xstate_T][a];

          // NOMERGE : CHECK FINAL STATE
          if (nomerge && L[Ystate_K] >= motif_span) {
#ifdef ASSERTB
            assert(YPstate_X ==  (Ystate_X^(1 << (Ystate_K-1))));
#endif
            Ystate_I = YPstate_I;
          } else {

            // UPDATE INDEX
            if (YPRstate_I >= 0 && SEEDPREFIXESMATCHING_X(YPRstate_I) == Ystate_X && SEEDPREFIXESMATCHING_T(YPRstate_I) == 0) { // if FX does gives a new prefix when compared with YP
              Ystate_I  = YPRstate_I;
              //cerr << "*";
            } else {
              if (SEEDPREFIXESMATCHING_X(YPstate_I) == Ystate_X && SEEDPREFIXESMATCHING_T(YPstate_I) == 0) {
                Ystate_I  = YPstate_I;
                //cerr << "+";
              } else {
                Y_not_there = true;
              }
            }
          }
#ifdef BUILD
          cerr << "*U "<< endl;
          cerr << "Ystate_X:"  <<  Ystate_X << ",  Ystate_K:" << Ystate_K<< endl;
#endif
        }// if ( Xstate_K > 0 )

        Ystate_T  = 0;

      }// if (a == gv_alphabet_size - 1)


      int final = ((Ystate_T + L[Ystate_K]) >= (motif_span - 1)) ? TRUE : FALSE;

      if (!nomerge  && final) {
        Ystate_I = 0;
      } else {
        // else check if the state has not been already created
        if (Y_not_there) {
          // create it and push it inside the list of states which need update.
          Ystate_I = addNewState(final);
          SEEDPREFIXESMATCHING_REVMAXXP(YPstate_I) = Ystate_I;
          statesSeedPrefixesMatchingSet.push_back(SeedPrefixesMatchingSet(Ystate_X,Ystate_T,Ystate_K,YPstate_I,0));
          statesNbRemaining.push(Ystate_I);

          VERB_FILTER(VERBOSITY_DEBUGGING, INFO__("$push I:"<< Ystate_I << " < X:" << hex << Ystate_X << ", t:" << dec << Ystate_T << ", k:" << Ystate_K << " >"););

        }
      }

      VERB_FILTER(VERBOSITY_DEBUGGING, INFO__("> add transition ( a:" << dec << a << ", q1:" << dec << Xstate_I << ", q2:" << Ystate_I << " )"););

      addNewTransition(a,Xstate_I,Ystate_I);
#ifdef ASSERTB
      assert(DETERMINISTIC(Xstate_I,a));
#endif
    }// for all letter in ALPHABET.
  }// while states remain.

  // Free unused data needed to build the automaton
  L.clear();
  FX.clear();
  Fk.clear();
  EX.clear();
  Ek.clear();
  statesSeedPrefixesMatchingSet.clear();
  return 0;
}

/** @class SeedPrefixesMatchingDMSet
 *  @brief simple "inner class" to keep, for a set of seeds of size @f$ n @f$, the states @f$ <X[1..n],t> @f$ and their maximal prefix @f$ k[1..n] @f$, together with the coverage inside a vector @f$ p @f$, 
 *         during SeedPrefixMatching_CoverageDM method
 *  @see Automaton_SeedPrefixesMatching_CoverageDM
 */

class SeedPrefixesMatchingDMSet {

public:
  /// set of prefixes matching for each seed
  vector<int> X;
  /// lenght of the last run of '1'
  short t;
  /// length of the maximal prefix in the set @f$ X @f$ : @f$ k = max(X) @f$
  vector<short> k;
  /** set of seed matching symbols as a suffix of the maximal ongoing seed prefix
   *  this vector can be "padded right" by run of '0' (uncovered),
   *  but "padded left" run of '0' must be removed.
   */
  vector<int> p;

  /// final state value (here can be 0 or any coverage value > 0)
  int final;
  /** @brief build a SeedPrefixesMatchingSet
   */
  SeedPrefixesMatchingDMSet(vector<int> X,  vector<short> k, vector<int> p, short t, int final) : t(t), final(final) {
    this->X = vector<int>(X);
    this->k = vector<short>(k);
    this->p = vector<int>(p);
  };
  /** @brief clear a SeedPrefixesMatchingSet
   */
  ~SeedPrefixesMatchingDMSet() {
    this->X.clear();
    this->k.clear();
    this->p.clear();
  };


  /** @brief compare two  SeedPrefixesMatchingDMSet (used for the map<SeedPrefixesMatchingDMSet>)
   *  @param lhs is the left SeedPrefixesMatchingDMSet to be compared
   *  @param rhs is the right SeedPrefixesMatchingDMSet to be compared
   *  @return true iif (lhs < rhs)
   */

  friend bool operator<(const SeedPrefixesMatchingDMSet& lhs, const SeedPrefixesMatchingDMSet& rhs) {
    /* compare T size */
    if (lhs.t < rhs.t)
      return true;
    else
      if (lhs.t > rhs.t)
        return false;

    /* then compare "k" size first */
    for (unsigned x = 0 ; x < lhs.k.size() ; x++) {
      if (lhs.k[x] < rhs.k[x]) {
        return true;
      } else {
        if (lhs.k[x] > rhs.k[x]) {
          return false;
        }
      }
    }

    /* then go to "X" details */
    for (unsigned x = 0 ; x < lhs.k.size() ; x++) {
      if (lhs.X[x] < rhs.X[x]) {
        return true;
      } else {
        if (lhs.X[x] > rhs.X[x]) {
          return false;
        }
      }
    }
    /* then compare "P" size */
    if (lhs.p.size() < rhs.p.size())
      return true;
    else
      if (lhs.p.size() > rhs.p.size())
        return false;

    /* then compare "P" elements */
    for (unsigned i = 0 ; i < lhs.p.size() ; i++) {
      if (lhs.p[i] < rhs.p[i]) {
        return true;
      } else {
        if (lhs.p[i] > rhs.p[i]) {
          return false;
        }
      }
    }
    /* then finally : compare final */
    if (lhs.final < rhs.final)
      return true;
    else
      if (lhs.final > rhs.final)
        return false;

    return false;
  };
};

#define SEEDPREFIXESMATCHINGDM_S(s)          (statesSeedPrefixesMatchingDMSet[(s)])
#define SEEDPREFIXESMATCHINGDM_X(s)          (statesSeedPrefixesMatchingDMSet[(s)].X)
#define SEEDPREFIXESMATCHINGDM_T(s)          (statesSeedPrefixesMatchingDMSet[(s)].t)
#define SEEDPREFIXESMATCHINGDM_K(s)          (statesSeedPrefixesMatchingDMSet[(s)].k)
#define SEEDPREFIXESMATCHINGDM_P(s)          (statesSeedPrefixesMatchingDMSet[(s)].p)
#define SEEDPREFIXESMATCHINGDM_FINAL(s)      (statesSeedPrefixesMatchingDMSet[(s)].final)



#define RESIZE() {                                                            \
   int _max_pref_ = 0;                                                        \
   for (int x = 0; x < nb_motifs; x++) {                                      \
     int _pref_ = MIN(L[x][Ystate_K[x]]+Ystate_T,motifs_span[x]);             \
         _max_pref_ = MAX(_max_pref_,_pref_);                                 \
   }                                                                          \
   int _delta_ = (int)Ystate_P.size() - _max_pref_;                           \
   if (_delta_ > 0) {                                                         \
     Ystate_P.erase (                                                         \
       Ystate_P.begin(),                                                      \
       Ystate_P.begin() + _delta_ - 1                                         \
     );                                                                       \
   }                                                                          \
}


/*
 * Do a "first" Normalization step (can be removed but this reduce the number of states)
 * - it checks if the coverage value at some position is already higher than the maximal
 *   one than can be reached by the ongoing hit or future hits : as such (since it has
 *   already be counted before) it can be erased.
 */

#define NORMALIZE() {                                                                               \
   int _P_size = Ystate_P.size();                                                                   \
   vector<int> P_reachable(_P_size,0);                                                              \
   for (int x = 0; x < nb_motifs; x++) {                                                            \
     /* (1.a) lower "t" part */                                                                     \
     int _t_max_ = MIN(Ystate_T, motifs_span[x]);                                                   \
     for (int _t_ = 1; _t_ <= _t_max_; _t_++) {                                                     \
       for (int _u_ = MAX(0, _t_ - _P_size) ; _u_ < _t_; _u_++) {                                   \
         P_reachable[_P_size - _t_ + _u_] =                                                         \
           MAX(weight_seed_alphabet[motifs[x][_u_]],                                                \
               P_reachable[_P_size - _t_ + _u_]);                                                   \
       }                                                                                            \
     }                                                                                              \
     /* (1.b) upped part without "t" */                                                             \
     for (int _k_ = 0; _k_ < Ystate_K[x] && (L[x][_k_+1] + Ystate_T <= motifs_span[x]); _k_++) {    \
       if (Ystate_X[x] & (1 << _k_)) {                                                              \
         int _r_ = L[x][_k_+1] + Ystate_T + 1;                                                      \
         for (int _u_ = MAX(0, _r_ - _P_size); _u_ < _r_; _u_++) {                                  \
           P_reachable[_P_size - _r_ + _u_] =                                                       \
             MAX(weight_seed_alphabet[motifs[x][_u_]],                                              \
                 P_reachable[_P_size - _r_ + _u_]);                                                 \
         }                                                                                          \
       }                                                                                            \
     }                                                                                              \
   }                                                                                                \
   /* (2) update and simplify */                                                                    \
   for (int _i_=0; _i_ < _P_size; _i_++) {                                                          \
     if (P_reachable[_i_] < Ystate_P[_i_]) {                                                        \
       Ystate_P[_i_] = P_reachable[_i_];                                                            \
     }                                                                                              \
   }                                                                                                \
}



/* Do a "second" Normalization step (can be removed but this reduce the number of states)
 * - it check if the coverage value at some position (before the "t" run, because no need
 *   to look inside it), can be "pushed right" to shorten the coverage string
 *
 *        v     v
 *        |     |
 *   #--#-@---#-@----..
 *    #--#-@---#-@---..
 *      #-#-@---#-@--..
 *        |     |#-@-..
 * if P =.2.....0..
 *        |     |
 * this can be replaced (and shortened)
 *        |     |
 * by P =.0.....2..
 *
 */

#define NORMALIZE2() {                                                                                    \
   int _P_size_ = Ystate_P.size();                                                                        \
   /* (1) upped part without "t" ("t" part not needed here since non equivalent) */                       \
   for (int _u_1_ = 0; _u_1_ < _P_size_-Ystate_T; _u_1_++) {                                              \
     if (Ystate_P[_u_1_]) {                                                                               \
       for (int _u_2_ = _u_1_+1; _u_2_ < _P_size_-Ystate_T; _u_2_++) {                                    \
         if (Ystate_P[_u_2_] < Ystate_P[_u_1_]) {                                                         \
           /* do positions _u_1_ and _u_2_ have equivalent seed elements under each seed prefix? */       \
           for (int x = 0; x < nb_motifs; x++) {                                                          \
             for (int _k_ = 0; _k_ < Ystate_K[x] && (L[x][_k_+1] + Ystate_T <= motifs_span[x]); _k_++) {  \
               if (Ystate_X[x] & (1 << _k_)) {                                                            \
                 /* length of the ongoing prefix */                                                       \
                 int _r_ = L[x][_k_+1] + Ystate_T + 1;                                                    \
                 /* motifs pos (if correct) */                                                            \
                 int _i_1_ = _r_ - _P_size_ + _u_1_;                                                      \
                 int _i_2_ = _r_ - _P_size_ + _u_2_;                                                      \
                 /* motifs elements (jocker if outside) */                                                \
                 int _s_1_ = ((_i_1_<0) || (_i_1_>=_r_)) ? 0 : motifs[x][_i_1_];                          \
                 int _s_2_ = ((_i_2_<0) || (_i_2_>=_r_)) ? 0 : motifs[x][_i_2_];                          \
                 if (weight_seed_alphabet[_s_1_] != weight_seed_alphabet[_s_2_]) { goto non_equiv; }      \
               }                                                                                          \
             }                                                                                            \
           }                                                                                              \
           /* equiv : exchange _u_2_ and _u_1_ elements */                                                \
           {                                                                                              \
             int u = Ystate_P[_u_2_];                                                                     \
             Ystate_P[_u_2_] = Ystate_P[_u_1_];                                                           \
             Ystate_P[_u_1_] = u;                                                                         \
           }                                                                                              \
         non_equiv:;                                                                                      \
         }                                                                                                \
       }                                                                                                  \
     }                                                                                                    \
   }                                                                                                      \
 }


/*
 */



#define REMOVE_LEFT_ZEROS() {             \
  unsigned _k_index_ = 0;                 \
  while (_k_index_ < Ystate_P.size()) {   \
    if (Ystate_P[_k_index_] > 0)          \
      break;                              \
    _k_index_++;                          \
  };                                      \
  if (_k_index_)                          \
    Ystate_P.erase (                      \
      Ystate_P.begin(),                   \
      Ystate_P.begin()+_k_index_          \
    );                                    \
}



#define FINAL_COUNT_AND_SET(Ystate_P,Ystate_final,x) {                           \
  int _i_ = motifs_span[x]-1, _j_ = Ystate_P.size()-1;                           \
  for (; (_i_ >= 0) && (_j_ >= 0) ; _i_--, _j_--) {                              \
    int _u_ = weight_seed_alphabet[motifs[x][_i_]] - Ystate_P[_j_];              \
    if (_u_ > 0) {                                                               \
      Ystate_final += _u_;                                                       \
      Ystate_P[_j_] += _u_;                                                      \
    }                                                                            \
  }                                                                              \
  if (_i_ >= 0 && _j_ < 0) { /* when the vector was left shrink  */              \
    while (_i_ >= 0) {                                                           \
      int _u_ = weight_seed_alphabet[motifs[x][_i_]];                            \
      Ystate_P.insert(Ystate_P.begin(),_u_);                                     \
      Ystate_final += _u_;                                                       \
      _i_--;                                                                     \
    }                                                                            \
  }                                                                              \
}



int automaton::Automaton_SeedPrefixesMatching_CoverageDM(const vector<seed *> & s,
                                                         const vector<int> weight_seed_alphabet,
                                                         const vector< vector <int> > & matchingmatrix) {

  VERB_FILTER(VERBOSITY_DEBUGGING, INFO__("== SeedPrefixesMatching_CoverageDM == "););

  const int  nb_motifs = s.size();
  vector<const int *>  motifs;
  vector<int>     motifs_span;
  int max_motifs_span = 0;

  // Get spans and patterns for all seeds
  for (int x = 0; x < nb_motifs; x++) {
    motifs_span.push_back(s[x]->span());
    motifs.push_back(s[x]->table());
    max_motifs_span = MAX(max_motifs_span,motifs_span[x]);
#ifdef BUILD
    cerr << "Seed " << x << endl << (*s[x]) << endl;
#endif
  }

#ifdef ASSERTB
  for (int x = 0; x < nb_motifs; x++) {
    if (motifs_span[x] <= 0) {
      _ERROR("SeedPrefixesMatching","null or negative span");
    }

    for (int i = 0; i < motifs_span[x]; i++) {
      if ( motifs[x][i] < 0 || motifs[x][i] >= gv_seed_alphabet_size) {
        _ERROR("SeedPrefixesMatching","incorrect seed element");
      }
    }
  }
#endif

#ifdef BUILD
  for (int x = 0; x < nb_motifs; x++) {
    cerr << " motif_span :" << motifs_span[x] << ", motif:";
    for (int i = 0; i < motifs_span[x]; i++ )
      cerr << motifs[x][i] << " ";
    cerr << endl;
  }
#endif


  // Compute xset_bitsize and L for all seeds
  vector < vector<int> > L(nb_motifs);
  vector <int> xset_bitsize(nb_motifs,0); // bitsize of the set
  for (int x = 0; x < nb_motifs; x++) {
    L[x] = vector<int>(motifs_span[x]+1,0);
    L[x][0] = -1;
#ifdef BUILD
    cerr << " L[0] : " << L[x][0] << endl;
#endif
    for (int i = 0; i < motifs_span[x]; i++) {
      if (motifs[x][i] < (gv_seed_alphabet_size - (gv_matching_symbol_flag?1:0))) {
        xset_bitsize[x]++;
        L[x][xset_bitsize[x]] = i;
#ifdef BUILD
        cerr << " L[" << xset_bitsize[x] << "] : " << L[x][xset_bitsize[x]] << endl;
#endif
      }
    }
#ifdef BUILD
    {
      int xset_size = 1 << xset_bitsize[x];
      cerr << " xset_size : "     << xset_size << ",\t" <<
              " xset_bitsize : " << xset_bitsize[x] << endl;
    }
#endif
  }

  // 1) Precompute E,F functions
  vector< vector < vector<int> > > FX(nb_motifs);
  vector< vector < vector<int> > > Fk(nb_motifs);
  vector< vector < vector < vector<int> > > > EX(nb_motifs);
  vector< vector < vector < vector<int> > > > Ek(nb_motifs);


  for (int x = 0; x < nb_motifs; x++) {
    // for the seed x=$0, gives the full value 'X' after a run of 1^t.a   ;   where t=$1 and a=$2
    FX[x] = vector < vector<int> >(motifs_span[x]+1,vector<int>(gv_align_alphabet_size,0));
    Fk[x] = vector < vector<int> >(motifs_span[x]+1,vector<int>(gv_align_alphabet_size,0));
    // F[x][t][a]
    for (int t = 0; t < motifs_span[x]+1; t++) {
      for (int i = 1; i <= xset_bitsize[x] && L[x][i] <= t; i++) {
        int r = 1<<(i-1);
        int b = motifs[x][L[x][i]];

        for (int a = 0; a < gv_align_alphabet_size; a++){
          if (MATCHES_AB(a,b)){
            FX[x][t][a] = (t>0) ? (FX[x][t-1][a]) | r : r;
            Fk[x][t][a] = i;
          }
        }
      }
    }
    // for the seed x=$0, update the (non '#') matching state k of 'X' after a run 1^t.a   ;   where t=$1, k=$2, a=$3
    EX[x] = vector< vector < vector<int> > >(motifs_span[x]+1, vector< vector<int> >((xset_bitsize[x]+1), vector<int>(gv_align_alphabet_size,0)));
    Ek[x] = vector< vector < vector<int> > >(motifs_span[x]+1, vector< vector<int> >((xset_bitsize[x]+1), vector<int>(gv_align_alphabet_size,0)));
    // E[x][t][k][a]
    for (int i = 1; i <= xset_bitsize[x]; i++) {
      int L_i = L[x][i];
      int b   = motifs[x][L_i];

      for (int a = 0; a < gv_align_alphabet_size; a++) {
        // does b match the mismatch a
        if (MATCHES_AB(a,b)) {
          // if X[ib]-t-1 is a joker position, do a mask to have its bit position.
          for (int k = 1; k < i; k++) {
            int t = L[x][i] - L[x][k] - 1;
#ifdef ASSERTB
            assert(t>=0);
            assert(t<motifs_span[x]);
            assert(k>0);
            assert(k<=xset_bitsize[x]);
#endif
            EX[x][t][k][a] = 1<<(i-1);
            Ek[x][t][k][a] = i;
          }
        }
      }
    }
  }

  // 2) automaton build
  // queue/stack used to store non preprocessed states <X,t> code
  queue<int> statesNbRemaining;

  // keep each state information (X,t,k) inside this table
  vector<SeedPrefixesMatchingDMSet> statesSeedPrefixesMatchingDMSet;

  // keep the index for each state in this map
  map<SeedPrefixesMatchingDMSet, int> statesNbIndex;


  // create a first state [0] and put it as the final one
  int Finalstate_I = addNewState(TRUE);
  statesSeedPrefixesMatchingDMSet.push_back(SeedPrefixesMatchingDMSet(vector<int>(nb_motifs,0),vector<short>(nb_motifs,0),vector<int>(),0,1));
  selfLoop(Finalstate_I);

  // create a second state [1]  : it will be the initial one
  int Initstate_I   = addNewState();
  SeedPrefixesMatchingDMSet s_i = SeedPrefixesMatchingDMSet(vector<int>(nb_motifs,0),vector<short>(nb_motifs,0),vector<int>(),0,0);
  statesSeedPrefixesMatchingDMSet.push_back(s_i);

  // push it to the Stack and the Map
  statesNbRemaining.push(Initstate_I);
  statesNbIndex[s_i] = Initstate_I;

  while (!statesNbRemaining.empty()) {

    // (A) Get current state Xstate
    int Xstate_I = statesNbRemaining.front();
    statesNbRemaining.pop();

#ifdef ASSERTB
    assert(Xstate_I >= 0);
    assert(Xstate_I <  (int)_states.size());
#endif

    vector<int>   Xstate_X = SEEDPREFIXESMATCHINGDM_X(Xstate_I);
    short         Xstate_T = SEEDPREFIXESMATCHINGDM_T(Xstate_I);
    vector<short> Xstate_K = SEEDPREFIXESMATCHINGDM_K(Xstate_I);
    vector<int>   Xstate_P = SEEDPREFIXESMATCHINGDM_P(Xstate_I);
    int       Xstate_final = SEEDPREFIXESMATCHINGDM_FINAL(Xstate_I);
#ifdef ASSERTB
    for (int x = 0; x < nb_motifs; x++) {
      assert(Xstate_X[x] >= 0);
      assert(Xstate_X[x] < (1<<xset_bitsize[x]));
      assert(Xstate_T >= 0);
      assert(Xstate_T <= max_motifs_span);
      assert(Xstate_K[x] >= 0);
      assert(Xstate_K[x] <= xset_bitsize[x]);
    }
#endif



VERB_FILTER(VERBOSITY_DEBUGGING, INFO__(
    "$pop I:" << Xstate_I << endl;
    cerr << "<";
    for (int x = 0; x < nb_motifs; x++) {
      cerr << " X[" << x << "]:" << hex << Xstate_X[x] << ", k[" << x << "]:" << dec << Xstate_K[x] << " ; ";
    }
    cerr << " { t:" << Xstate_T << " } ; ";
    cerr << " { P:";
    for (unsigned u = 0; u < Xstate_P.size(); u++) {
      cerr << Xstate_P[u];
    }
    cerr << " }  ; final:" << Xstate_final << " >";
););

    // for each letter (a) : compute the transition (Xstate)--a-->Ystate [and determine Ystate]
    for (int a = 0; a < gv_align_alphabet_size; a++) {

      /* check if the a == allways_matching symbol apply ? */
      bool a_is_one = a == (gv_align_alphabet_size - (gv_matching_symbol_flag?1:0));

      VERB_FILTER(VERBOSITY_DEBUGGING, INFO__("a = " << a););

      // (B) try to compute Ystate_s (use a "copy" to get the correct size for all allocators + the size "t")
      SeedPrefixesMatchingDMSet Ystate_s = SeedPrefixesMatchingDMSet(vector<int>(Xstate_X),vector<short>(Xstate_K),vector<int>(Xstate_P),Xstate_T,0);
      vector<int>             & Ystate_X     = Ystate_s.X;
      short                   & Ystate_T     = Ystate_s.t;
      vector<short>           & Ystate_K     = Ystate_s.k;
      vector<int>             & Ystate_P     = Ystate_s.p;
      int                     & Ystate_final = Ystate_s.final;

      if ( a_is_one ) { // (a == match)
        // (B.1) matches (update the final)
        Ystate_T++;
#ifdef BUILD
        cerr << "\t\t" << "after T increase : " << Xstate_T << " -> " << Ystate_T << endl;
        for (int x = 0; x < nb_motifs; x++) {
          cerr << "\t\t" << "Ystate_X[" << x << "]:"  << Ystate_X[x]  << ",  Ystate_K[" << x << "]:"  << Ystate_K[x]  << endl;
        }
#endif
      } else { // (a != match)
        for (int x = 0; x < nb_motifs; x++) {
          // init Y_X[x] Y_K[x] both to 0
          Ystate_X[x] = 0;
          Ystate_K[x] = 0;
          // (B.2.1) if prefixes from X_state already matches
          if (Xstate_K[x] > 0) {
            // recompute the full set (local memory here)
            for (int i = 1 ; i <= Xstate_K[x] ; i++) {
              // compute for each maching prefix "i" followed by "1^t.a" ones
              if ((Xstate_X[x]) & 1<<(i-1)) {
                Ystate_X[x] =      Ystate_X[x] | EX[x][MIN(Xstate_T,motifs_span[x])][i][a];
                Ystate_K[x] = MAX( Ystate_K[x] , Ek[x][MIN(Xstate_T,motifs_span[x])][i][a]);
#ifdef BUILD
                cerr << "\t\t" << "after EX {Pref:" << i << "}.1^" << Xstate_T << "." << a << endl;
                cerr << "\t\t" << "Ystate_X[" << x << "]:"  << Ystate_X[x]  << ",  Ystate_K[" << x << "]:" << Ystate_K[x] << endl;
#endif
#ifdef ASSERTB
                assert(Ystate_X[x] < (1<<xset_bitsize[x]));
                assert(Ystate_K[x] <= xset_bitsize[x]);
#endif
              }
            }
          }
          // (B.2.2) merge the result with the "1^t.a"
          Ystate_X[x] =      Ystate_X[x] | FX[x][MIN(Xstate_T,motifs_span[x])][a];
          Ystate_K[x] = MAX( Ystate_K[x] , Fk[x][MIN(Xstate_T,motifs_span[x])][a]);
#ifdef BUILD
          cerr << "\t\t" << "after FX on 1^"<< Xstate_T << "." << a << endl;
          cerr << "\t\t" << "Ystate_X[" << x << "]:" << Ystate_X[x] << ",  Ystate_K[" << x << "]:" << Ystate_K[x] << endl;
#endif
#ifdef ASSERTB
          assert(Ystate_X[x] < (1<<xset_bitsize[x]));
          assert(Ystate_K[x] <= xset_bitsize[x]);
#endif
        }
        Ystate_T = 0;
      }// if (a_is_one)

      // NOW : add a new element to the right of the P vector
      Ystate_P.push_back(0);
#ifdef BUILD
      cerr << "\t\t" << "after push_back on Ystate_P " << endl;
      cerr << "\t\t { P:";
      for (unsigned u = 0; u < Ystate_P.size(); u++) {
        cerr << "\t\t" << Ystate_P[u];
      }
      cerr << "\t\t } >" << endl;
#endif

      // (C) Limiter for all cases encountered / then final for the "exact" ones
      // (C.1) Limit run when oveflow the max_motif_span
      Ystate_final = 0;
      if (Ystate_T > max_motifs_span) {
        Ystate_T = max_motifs_span;
        // (C.2) Final computation and mask fixing
        // BEFORE :
        // final += nb_motifs;
        // NOW :
#ifdef BUILD
        cerr << "\t\t" << "after long 1^t ; { t:" << Ystate_T << "}" << endl;
#endif
      for (int x = 0; x < nb_motifs; x++) {
        FINAL_COUNT_AND_SET(Ystate_P,Ystate_final,x);
      }
#ifdef BUILD
      cerr << "\t\t" << "after final count on long 1^t " << endl;
      cerr << "\t\t { P:";
      for (unsigned u = 0; u < Ystate_P.size(); u++) {
        cerr << "\t\t" << Ystate_P[u];
      }
      cerr << "\t\t } ; final:" << Ystate_final << " >" << endl;
#endif
      } else {
        for (int x = 0; x < nb_motifs; x++) {
          // if a seed "overmatch" ... and is thus too old ...
          if (L[x][Ystate_K[x]] + Ystate_T >= motifs_span[x]) {
            // and if the "overmatch" include a bit "1" from X : we can do some cleaning in X
            if (Ystate_K[x] > 0) {
              // remove the current maximal bit "1" in X, and find the previous "1" (if any)
              do {
                Ystate_X[x] &= ~(1<<(Ystate_K[x]-1));
                Ystate_K[x]--;
              } while ((Ystate_K[x] > 0) && !(Ystate_X[x] & (1<<(Ystate_K[x]-1))));
            }
          }
          // otherwise this is a "T too long" problem for this seed : its not a problem at the end

          // (C.2) Final computation and mask fixing
          // BEFORE :
          //final += ((Ystate_T + L[x][Ystate_K[x]]) >= (motifs_span[x]-1)) ? 1 : 0;
          // NOW :
          if ((Ystate_T + L[x][Ystate_K[x]]) >= (motifs_span[x]-1)) {
            FINAL_COUNT_AND_SET(Ystate_P,Ystate_final,x);
#ifdef BUILD
            cerr << "\t\t" << "after final count on [x:" << x << "] " << endl;
            cerr << "\t\t { P:";
            for (unsigned u = 0; u < Ystate_P.size(); u++) {
              cerr << "\t\t" << Ystate_P[u];
            }
            cerr << "\t\t } ; final:" << Ystate_final << " >" << endl;
#endif
          }
        }
      }
      //NOW : resize, normalize, then remove_left the P vector

      // resize
      RESIZE();
#ifdef BUILD
      cerr << "\t\t" << "after resize" << endl;
      cerr << "\t\t { P:";
      for (unsigned u = 0; u < Ystate_P.size(); u++) {
        cerr << "\t\t" << Ystate_P[u];
      }
      cerr << "\t\t } ; final:" << Ystate_final << " >" << endl;
#endif

      // normalize : theses two steps are optional !! but help reducing the automaton size ...

      // decrease/remove unecessary coverage values ...
      NORMALIZE();
#ifdef BUILD
      cerr << "\t\t" << "after normalize" << endl;
      cerr << "\t\t { P:";
      for (unsigned u = 0; u < Ystate_P.size(); u++) {
        cerr << "\t\t" << Ystate_P[u];
      }
      cerr << "\t\t } ; final:" << Ystate_final << " >" << endl;
#endif

      // search for equivalent positions : and push them right / sort them
      NORMALIZE2();
#ifdef BUILD
      cerr << "\t\t" << "after normalize2" << endl;
      cerr << "\t\t { P:";
      for (unsigned u = 0; u < Ystate_P.size(); u++) {
        cerr << "\t\t" << Ystate_P[u];
      }
      cerr << "\t\t } ; final:" << Ystate_final << " >" << endl;
#endif

      // remove left ending zeros
      REMOVE_LEFT_ZEROS();
#ifdef BUILD
      cerr << "\t\t" << "after remove left" << endl;
      cerr << "\t\t { P:";
      for (unsigned u = 0; u < Ystate_P.size(); u++) {
        cerr << "\t\t" << Ystate_P[u];
      }
      cerr << "\t\t } ; final:" << Ystate_final << " >" << endl;
#endif

      // (D) Insert if Ystate does not exists
      int Ystate_I = 0;
      if (statesNbIndex.find(Ystate_s) == statesNbIndex.end()) {
        // create it and push it inside the list of states which need update.
        Ystate_I = addNewState(Ystate_final);
        statesNbIndex[Ystate_s] = Ystate_I;
        statesSeedPrefixesMatchingDMSet.push_back(SeedPrefixesMatchingDMSet(vector<int>(Ystate_X),vector<short>(Ystate_K),vector<int>(Ystate_P),Ystate_T,Ystate_final));
        statesNbRemaining.push(Ystate_I);

VERB_FILTER(VERBOSITY_DEBUGGING, INFO__(
        "$push I:"<< Ystate_I << endl;
        cerr << "<";
        for (int x = 0; x < nb_motifs; x++) {
          cerr << " Y[" << x << "]:" << hex << Ystate_X[x] << ", k[" << x << "]:" << dec << Ystate_K[x] << " ; ";
        }
        cerr << " { t:" <<  Ystate_T << " } ; ";
        cerr << " { P:";
        for (unsigned u = 0; u < Ystate_P.size(); u++) {
          cerr << Ystate_P[u];
        }
        cerr << " } ; final:" << Ystate_final << " >";
););

      } else {
        Ystate_I = statesNbIndex[Ystate_s];
        VERB_FILTER(VERBOSITY_DEBUGGING, INFO__("$exist I:" << Ystate_I););
      }

      VERB_FILTER(VERBOSITY_DEBUGGING, INFO__("> add transition ( a:" << dec << a << ", q1:" << dec << Xstate_I << ", q2:" << Ystate_I << " )"););

      addNewTransition(a,Xstate_I,Ystate_I);
#ifdef ASSERTB
      assert(DETERMINISTIC(Xstate_I,a));
#endif
    } // for all letter in ALPHABET.
  } // while states remain.

  // Free unused data needed to build the automaton
  L.clear();
  FX.clear();
  Fk.clear();
  EX.clear();
  Ek.clear();
  statesNbIndex.clear();
  statesSeedPrefixesMatchingDMSet.clear();
  return 0;
}





/** @class SeedPrefixesMatchingCostSet
 *  @brief simple "inner class" to keep states @f$ <X,t> @f$ and  @f$ k + Cost @f$
 *         during SeedPrefixMatching method
 *  @see Automaton_SeedPrefixesMatchingCost
 */

class SeedPrefixesMatchingCostSet {

public:
  /// set of prefixes matching
  int X;
  /// lenght of the last run of '1'
  short t;
  /// length of the maximal prefix in the set @f$ X @f$ : @f$ k = max(X) @f$
  short k;
  /// gives the state @f$< X' = X / max(X), t > @f$ (this state does exists and is created before @f$ <X,t> @f$)
  int Xp;
  /** @brief gives the last state @f$ <X'',t> @f$ created that verifies :
   *   @li @f$ X''/max(X'') = X @f$
   *   @li @f$ max(X'') @f$ is greatest among all created @f$ X'' @f$ verifying the previous condition
   */
  int RevMaxXp;
  /// gives the Cost at this state
  int Cost;
  /** @brief build a SeedPrefixesMatchingCostSet
   */
  SeedPrefixesMatchingCostSet(int X = 0, short t = 0, short k = 0, int Xp = 0, int RevMaxXp = 0, int Cost = 0) : X(X), t(t), k(k), Xp(Xp), RevMaxXp(RevMaxXp), Cost(Cost) {};
};

#define SEEDPREFIXESMATCHINGCOST_X(s)          (statesSeedPrefixesMatchingCostSet[(s)].X)
#define SEEDPREFIXESMATCHINGCOST_T(s)          (statesSeedPrefixesMatchingCostSet[(s)].t)
#define SEEDPREFIXESMATCHINGCOST_K(s)          (statesSeedPrefixesMatchingCostSet[(s)].k)
#define SEEDPREFIXESMATCHINGCOST_XP(s)         (statesSeedPrefixesMatchingCostSet[(s)].Xp)
#define SEEDPREFIXESMATCHINGCOST_REVMAXXP(s)   (statesSeedPrefixesMatchingCostSet[(s)].RevMaxXp)
#define SEEDPREFIXESMATCHINGCOST_COST(s)       (statesSeedPrefixesMatchingCostSet[(s)].Cost)



int automaton::Automaton_SeedPrefixesMatchingCost(const seed& s,
                                                  const vector< vector <int> > & matchingmatrix,
                                                  const bool nomerge,
                                                  const vector<int> & costs,
                                                  const int cost_threshold) {

  VERB_FILTER(VERBOSITY_DEBUGGING, INFO__("== SeedPrefixesMatchingCost == (nomerge:" << dec << nomerge << ")"););

  // Get the seed span
  int   motif_span = s.span();
  int * motif      = s.table();

#ifdef ASSERTB
  if (motif_span <= 0) {
    _ERROR("SeedPrefixesMatchingCost","null or negative span");
  }

  for (int i = 0; i < motif_span; i ++)
    if ( motif[i] < 0 || motif[i] >= gv_seed_alphabet_size)
      _ERROR("SeedPrefixesMatchingCost","incorrect seed element");
#endif

#ifdef BUILD
  cerr << " motif_span :" << motif_span << ", motif:";
  for (int i = 0; i < motif_span; i++ )
    cerr << motif[i] << " ";
  cerr << endl;
#endif


  // Compute xset_bitsize and X
  vector<int> L(motif_span+1,0);
  int xset_bitsize = 0; // bitsize of the set
  L[0] = -1;
  for (int i = 0; i < motif_span; i++) {
    if (motif[i] < (gv_seed_alphabet_size - (gv_matching_symbol_flag?1:0))) {
      xset_bitsize++;
      L[xset_bitsize] = i;
#ifdef BUILD
      cerr << " L[" << xset_bitsize << "] : " << L[xset_bitsize] << endl;
#endif
    }
  }

#ifdef BUILD
  {
    int xset_size = 1 << xset_bitsize;
    cerr << " xset_size :"     << xset_size << ",\t" <<
      " xset_bitsize : " << xset_bitsize << endl;
  }
#endif


  // 1) Precompute E,F functions
  vector< vector<int> > FX(motif_span+(nomerge?1:0),vector<int>(gv_align_alphabet_size,0));
  vector< vector<int> > Fk(motif_span+(nomerge?1:0),vector<int>(gv_align_alphabet_size,0));
  vector< vector< vector<int> > > EX(motif_span+(nomerge?1:0), vector< vector<int> >((xset_bitsize+1), vector<int>(gv_align_alphabet_size,0)));
  vector< vector< vector<int> > > Ek(motif_span+(nomerge?1:0), vector< vector<int> >((xset_bitsize+1), vector<int>(gv_align_alphabet_size,0)));

  // F[t][a]
  for (int t = 0; t < motif_span+(nomerge?1:0); t++) {
    for (int i = 1; i<= xset_bitsize && L[i] <= t; i++) {
      int r = 1<<(i-1);
      int b = motif[L[i]];

      for (int a = 0; a < gv_align_alphabet_size; a++){
        if (MATCHES_AB(a,b)){
          FX[t][a] = (t>0) ? (FX[t-1][a]) | r : r;
          Fk[t][a] = i;
        }
      }
    }
  }


  // E[t][a]
  for (int i = 1; i<=xset_bitsize; i++) {
    int L_i = L[i];
    int b   = motif[L_i];

    for (int a = 0; a < gv_align_alphabet_size; a++) {
      // does b match the mismatch a
      if (MATCHES_AB(a,b)) {
        // if X[ib]-t-1 is a joker position, do a mask to have its bit position.
        for (int k = 1; k < i; k++) {
          int t = L[i] - L[k] - 1;
#ifdef ASSERTB
          assert(t>=0);
          assert(t<motif_span);
          assert(k>0);
          assert(k<=xset_bitsize);
#endif
          EX[t][k][a] = 1<<(i-1);
          Ek[t][k][a] = i;
          //
        }
      }
    }
  }


  // 2) automaton build
  // queue/stack used to store non preprocessed states <X,t> code
  queue<int>  statesNbRemaining;
  // keep each state information (X,t,k) inside this table
  vector<SeedPrefixesMatchingCostSet> statesSeedPrefixesMatchingCostSet(0);

  // create a first state [0] and put it as the final one
  int Finalstate_I = addNewState(TRUE);
  statesSeedPrefixesMatchingCostSet.push_back(SeedPrefixesMatchingCostSet(0,0,0,0,0,INT_INFINITY));
  selfLoop(Finalstate_I);

  // create a second state [1]  : it will be the initial one
  int Initstate_I = addNewState();
  statesSeedPrefixesMatchingCostSet.push_back(SeedPrefixesMatchingCostSet(0,0,0,1,1,0));

  // bag reject set [2]
  int RejectBagstate_I  = addNewState();
  statesSeedPrefixesMatchingCostSet.push_back(SeedPrefixesMatchingCostSet(0,0,0,0,0,INT_INFINITY));
  selfLoop(RejectBagstate_I);

  // create first level states and push them in the queue
  int State_R = 0;
  for (int a = 0; a < gv_align_alphabet_size; a++) {
    bool matches = MATCHES_AB(a,motif[0]);
    int  final   =  ((motif_span == 1) && matches) ? TRUE : FALSE;
    int  cost    = SEEDPREFIXESMATCHINGCOST_COST(Initstate_I) + costs[a];
    if (cost > cost_threshold) {
      addNewTransition(a,Initstate_I,RejectBagstate_I);
    } else {
      if ( a == (gv_align_alphabet_size - (gv_matching_symbol_flag?1:0))) {
        if (final && !nomerge) {
          addNewTransition(a,Initstate_I,Finalstate_I);
        } else {
          int State_I = addNewState(final);
          statesSeedPrefixesMatchingCostSet.push_back(SeedPrefixesMatchingCostSet(0,1,0,Initstate_I,-1,cost));
          statesNbRemaining.push(State_I);
          addNewTransition(a,Initstate_I,State_I);
        }
      } else {
        if (matches) {
          if (final && !nomerge) {
            State_R = Finalstate_I;
          } else {
            if (!State_R) {
              State_R = addNewState(final);
              statesSeedPrefixesMatchingCostSet.push_back(SeedPrefixesMatchingCostSet(1,0,1,Initstate_I,-1,cost));
              statesNbRemaining.push(State_R);
            }
          }
          addNewTransition(a,Initstate_I,State_R);
        } else {
          addNewTransition(a,Initstate_I,Initstate_I);
        }
      }
    }
  }



  // main loop
  while (!statesNbRemaining.empty()) {

    // current state Xstate
    int Xstate_I = statesNbRemaining.front();
    statesNbRemaining.pop();

#ifdef ASSERTB
    assert(Xstate_I >= 0);
    assert(Xstate_I <  (int)_states.size());
#endif

    int Xstate_X      = SEEDPREFIXESMATCHINGCOST_X (Xstate_I);
    int Xstate_T      = SEEDPREFIXESMATCHINGCOST_T (Xstate_I);
    int Xstate_K      = SEEDPREFIXESMATCHINGCOST_K (Xstate_I);
    int Xstate_cost   = SEEDPREFIXESMATCHINGCOST_COST (Xstate_I);

#ifdef ASSERTB
    assert(Xstate_X >= 0);
    assert(Xstate_X < (1<<xset_bitsize));
    assert(Xstate_T >= 0);
    assert(Xstate_T < motif_span + (nomerge?1:0));
    assert(Xstate_K >= 0);
    assert(Xstate_K <= xset_bitsize);
    assert(Xstate_cost <= cost_threshold);
#endif

    // XPstate : defined as "Xstate" where the maximal matching prefix fails ...
    int XPstate_I   = SEEDPREFIXESMATCHINGCOST_XP(Xstate_I);

#ifdef ASSERTB
    assert(XPstate_I >= 0);
    assert(XPstate_I < (int)_states.size());
#endif

    VERB_FILTER(VERBOSITY_DEBUGGING, INFO__("$pop  I:" << Xstate_I << " < X:" << Xstate_X << ", t:" << Xstate_T <<", k:"<< Xstate_K <<" >"););

    // for each letter (a) : compute the transition (Xstate)--a-->Ystate [and determine Ystate]
    for (int a = 0; a < gv_align_alphabet_size; a++) {
      int Ystate_cost   = Xstate_cost + costs[a];

      VERB_FILTER(VERBOSITY_DEBUGGING, INFO__("a = " << a););

      if (Ystate_cost > cost_threshold) {
        addNewTransition(a,Xstate_I,RejectBagstate_I);
      } else {
        // set at true if Ystate has not been already created before
        bool Y_not_there = false;

        // YPstate : defined as (XPstate)--a-->
        int YPstate_I   = _states[XPstate_I]._next[a].begin()->_state;

#ifdef ASSERTB
        assert(YPstate_I >= 0);
        assert(YPstate_I < (int)_states.size());
#endif

        int YPstate_T   = SEEDPREFIXESMATCHINGCOST_T(YPstate_I);
        int YPstate_X   = SEEDPREFIXESMATCHINGCOST_X(YPstate_I);
        int YPstate_K   = SEEDPREFIXESMATCHINGCOST_K(YPstate_I);

#ifdef ASSERTB
        assert(YPstate_X >= 0);
        assert(YPstate_X < (1<<xset_bitsize));
        assert(YPstate_T >= 0);
        assert(YPstate_T < motif_span + (nomerge?1:0));
        assert(YPstate_K >= 0);
        assert(YPstate_K <= xset_bitsize);
#endif

        // YPRstate ...value is either -1 if undefined,
        // otherwise >= 0 to indicate the last state created that
        // verifies YPRstate / maximalprefixmatching{YPRstate} = YPstate.
        int YPRstate_I  = SEEDPREFIXESMATCHINGCOST_REVMAXXP(YPstate_I);

#ifdef ASSERTB
        assert(YPRstate_I >= (-1));
        assert(YPRstate_I < (int)_states.size());
#endif

        // Ystate_I : next state according to letter "a"
        // by default values are given here
        // and will be processed
        int Ystate_I    = Xstate_I;
        int Ystate_X    = Xstate_X;
        int Ystate_T    = Xstate_T;
        int Ystate_K    = Xstate_K;

        // (a == match)
        if ( a == (gv_align_alphabet_size - (gv_matching_symbol_flag?1:0))) {
          // UPDATE STATE VALUES
          Ystate_T++;         // (a) add a 1 to the run of 1 and
          // NOMERGE : CHECK FINAL STATE
          if (nomerge && L[Ystate_K] + Ystate_T >= motif_span) {
            if (Ystate_K > 0) {
              Ystate_I = YPstate_I;
              Ystate_T = YPstate_T;
              Ystate_X = YPstate_X;
              Ystate_K = YPstate_K;
            } else {
              Ystate_I = Xstate_I;
              Ystate_T = Xstate_T;
              Ystate_X = Xstate_X;
              Ystate_K = Xstate_K;
            }
          } else {
            //UPDATE STATE INDEX
            Y_not_there = true;
          }
        } else {
          // (a != match)
          // either (K > 0) or (K == 0)
          if ( Xstate_K > 0 ) {

#ifdef ASSERTB
            assert(DETERMINISTIC(XPstate_I,a));
#endif
            // (b) add a possible prefix when compared with XP ??
            //   - either EX[Xstate_T][Xstate_K][a] == 0
            //     and  Y = YP : thus  Y already exists
            //   - otherwise EX = 1<<Ek and thus the possible Y may exist ( it can be YPR )

            if (EX[Xstate_T][Xstate_K][a]) {

#ifdef BUILD
              cerr << "EX != 0 : maximal prefix does match " << endl;
#endif
              // UPDATE STATE VALUES
              Ystate_X  =      YPstate_X | EX[Xstate_T][Xstate_K][a];
              Ystate_K  = MAX( YPstate_K , Ek[Xstate_T][Xstate_K][a]);

              // NOMERGE : CHECK FINAL STATE
              if (nomerge && L[Ystate_K] >= motif_span) {
#ifdef ASSERTB
                assert(YPstate_X ==  (Ystate_X^(1 << (Ystate_K-1))));
#endif
                goto Ex0;
              }

              // UPDATE STATE INDEX
              if (YPRstate_I >=0 && SEEDPREFIXESMATCHINGCOST_X(YPRstate_I) == Ystate_X /* COMMENT:NOT NEEDED && SEEDPREFIXESMATCHING_T(YPRstate_I) == 0 */ ) {
#ifdef ASSERTB
                assert(SEEDPREFIXESMATCHINGCOST_T(YPRstate_I) == 0);
#endif
                Ystate_I = YPRstate_I; // state already exists (has been processed before)
              } else {
                Y_not_there = true;     // state must be created for the first time
              }
            } else {
            Ex0:
#ifdef BUILD
              cerr << "Ex == 0 : maximal prefix doesnt match : Y is thus YP" << endl;
#endif

              // UPDATE STATE VALUES
              Ystate_X  =      YPstate_X;
              Ystate_K  =      YPstate_K;
              // NOMERGE : CHECK FINAL STATE
              // ** not needed **
              // UPDATE STATE INDEX
              Ystate_I  =      YPstate_I;
            }

#ifdef BUILD
            cerr << "**V " << endl;
            cerr << "YPstate_X:" << YPstate_X << ",  YPstate_K:" << YPstate_K << endl;
            cerr << "Ystate_X:"  << Ystate_X  << ",  Ystate_K:"  << Ystate_K  << endl;
#endif
          } else { // Xstate is empty ... and we have read 1^t.a  since ...

            // (c) initial case when reading 1^t . a
            //     three subcases can be deduced here :
            //     - either there is one prefix matching on the full 1^t.a span ? in this case does such state exists ?
            //     - either the maximal prefix does not
            //
#ifdef BUILD
            cerr << " FX on 1^"<< Xstate_T << "." << a << endl;
#endif

            // UPDATE STATE VALUES
            Ystate_X   = FX[Xstate_T][a];
            Ystate_K   = Fk[Xstate_T][a];

            // NOMERGE : CHECK FINAL STATE
            if (nomerge && L[Ystate_K] >= motif_span) {
#ifdef ASSERTB
              assert(YPstate_X ==  (Ystate_X^(1 << (Ystate_K-1))));
#endif
              Ystate_I = YPstate_I;
            } else {

              // UPDATE INDEX
              if (YPRstate_I >= 0 && SEEDPREFIXESMATCHINGCOST_X(YPRstate_I) == Ystate_X && SEEDPREFIXESMATCHINGCOST_T(YPRstate_I) == 0) { // if FX does gives a new prefix when compared with YP
                Ystate_I  = YPRstate_I;
                //cerr << "*";
              } else {
                if (SEEDPREFIXESMATCHINGCOST_X(YPstate_I) == Ystate_X && SEEDPREFIXESMATCHINGCOST_T(YPstate_I) == 0) {
                  Ystate_I  = YPstate_I;
                  //cerr << "+";
                } else {
                  Y_not_there = true;
                }
              }
            }
#ifdef BUILD
            cerr << "*U "<< endl;
            cerr << "Ystate_X:"  << Ystate_X << ",  Ystate_K:" << Ystate_K << endl;
#endif
          }// if ( Xstate_K > 0 )

          Ystate_T  = 0;

        }// if (a == gv_alphabet_size - 1)

        int final = (Ystate_T + L[Ystate_K] >= motif_span - 1) ? TRUE : FALSE;

        if (!nomerge  && final) {

          Ystate_I = 0;

        } else {

          // else check if the state has not been already created
          if (Y_not_there) {
            // create it and push it inside the list of states which need update.
            Ystate_I = addNewState(final);
            SEEDPREFIXESMATCHINGCOST_REVMAXXP(YPstate_I) = Ystate_I;
            statesSeedPrefixesMatchingCostSet.push_back(SeedPrefixesMatchingCostSet(Ystate_X,Ystate_T,Ystate_K,YPstate_I,0,Ystate_cost));
            statesNbRemaining.push(Ystate_I);

            VERB_FILTER(VERBOSITY_DEBUGGING, INFO__("$push I:"<< Ystate_I << " < X:" << hex << Ystate_X << ", t:" << dec << Ystate_T << ", k:" << Ystate_K << " >"););

          }
        }

        VERB_FILTER(VERBOSITY_DEBUGGING, INFO__("> add transition ( a:" << dec << a << ", q1:" << dec << Xstate_I << ", q2:" << Ystate_I << " )"););

        addNewTransition(a,Xstate_I,Ystate_I);

        statesSeedPrefixesMatchingCostSet[Ystate_I].Cost = MIN(Ystate_cost,statesSeedPrefixesMatchingCostSet[Ystate_I].Cost);


#ifdef ASSERTB
        assert(DETERMINISTIC(Xstate_I,a));
#endif
      }// YstateCost <= cost_threshold
    }// for all letter in ALPHABET.
  }// while states remain.

  // Free unused data needed to build the automaton
  L.clear();
  FX.clear();
  Fk.clear();
  EX.clear();
  Ek.clear();

  statesSeedPrefixesMatchingCostSet.clear();
  return 0;
}













/** @class SeedBuhlerSet
 *  @brief simple "inner class" to keep states properties
 *         during SeedBuhler (Aho Corasick) method
 *  @see Automaton_SeedBuhler
 */

class SeedBuhlerSet {

public:
  /// Aho-Corasick Fail function for this state
  int   Fail;
  /// Number of transition that need to be completed for this state
  short Rest;
  /// Level (distance from the root) in the AC-tree for this state
  short Level;
  /** @brief build a SeedBuhlerSet
   */
  SeedBuhlerSet(int Fail = 0, short Rest = 0, short Level = 0) : Fail(Fail), Rest(Rest), Level(Level) {};
};

#define SEEDBUHLER_FAIL(s)  (statesSeedBuhlerSet[(s)].Fail)
#define SEEDBUHLER_REST(s)  (statesSeedBuhlerSet[(s)].Rest)
#define SEEDBUHLER_LEVEL(s) (statesSeedBuhlerSet[(s)].Level)




int automaton::Automaton_SeedBuhler(const seed& s,
                                    const vector< vector <int> > & matchingmatrix,
                                    const bool nomerge) {

  // Get the seed span
  int   motif_span  = s.span();
  int * motif       = s.table();

#ifdef ASSERTB
  if (motif_span <= 0) {
    _ERROR("Automaton_SeedBuhler","negative span");
  }
#endif

  // (1) build the trie
  queue<int>  statesNbRemaining;
  // keep each state information (Fail,Rest,Level) inside this table
  vector<SeedBuhlerSet> statesSeedBuhlerSet(0);

  // create a first state [0] and put it as the final one
  int Finalstate_I = addNewState(TRUE);
  statesSeedBuhlerSet.push_back(SeedBuhlerSet(0,0,0));
  selfLoop(Finalstate_I);

  // create a second state [1] : it will be the initial one
  int Initstate_I = addNewState();
  statesSeedBuhlerSet.push_back(SeedBuhlerSet(0,0,0));
  statesNbRemaining.push(Initstate_I);

  while (!statesNbRemaining.empty()) {

    // current state remaining
    int stateN =  statesNbRemaining.front();
    statesNbRemaining.pop();
    int i        =  SEEDBUHLER_LEVEL(stateN);
    SEEDBUHLER_REST(stateN) = gv_align_alphabet_size;
    int final    = (i >= motif_span - 1) ? TRUE : FALSE;
    int b = motif[i];

    for (int a = 0; a < gv_align_alphabet_size; a++) {
      if (MATCHES_AB(a,b)){
        int StateX   = 0;
        if (i < motif_span - 1 || nomerge){
          StateX = addNewState(final); // add new state
          statesSeedBuhlerSet.push_back(SeedBuhlerSet(0,0,i+1));
          statesNbRemaining.push(StateX);
        }
        addNewTransition(a,stateN,StateX);
        SEEDBUHLER_REST(stateN)--;
      }//matches
    }//for a
  }//while states remaining


  // (2) add failure links

  // init first states
  int  r  = 1;
  SEEDBUHLER_FAIL(r) = r;
  for (int a = 0; a < SEEDBUHLER_REST(r); a++){
    addNewTransition(a,r,SEEDBUHLER_FAIL(r));
#ifdef ASSERTB
    assert(DETERMINISTIC(r,a));
#endif
  }

  for (int a = SEEDBUHLER_REST(r); a < gv_align_alphabet_size; a++) {
#ifdef ASSERTB
    assert(DETERMINISTIC(r,a));
#endif
    int p = _states[r]._next[a].begin()->_state;
    SEEDBUHLER_FAIL(p) = r;
  }

  // init last state
  SEEDBUHLER_FAIL(0) = 0;
  SEEDBUHLER_REST(0) = 0;

  // intermediate states:  BORDER function from Aho-Corasick*/
  for (int r = 2; r < size(); r++){
    for (int a = 0; a < SEEDBUHLER_REST(r); a++){
#ifdef ASSERTB
      assert(DETERMINISTIC(SEEDBUHLER_FAIL(r),a));
#endif
      //if (!(_states[r]._final)) /*FIXME || nomerge ??*/
      addNewTransition(a,r,_states[SEEDBUHLER_FAIL(r)]._next[a].begin()->_state);
#ifdef ASSERTB
      assert(DETERMINISTIC(r,a));
#endif
    }

    for (int a = SEEDBUHLER_REST(r); a < gv_align_alphabet_size; a++){
#ifdef ASSERTB
      assert(DETERMINISTIC(r,a));
#endif
      int p = _states[r]._next[a].begin()->_state;
      int s = SEEDBUHLER_FAIL(r);
#ifdef ASSERTB
      assert(DETERMINISTIC(s,a));
#endif
      SEEDBUHLER_FAIL(p) = _states[s]._next[a].begin()->_state;
    }
  }
  statesSeedBuhlerSet.clear();
  return 0;
}




/** @class SeedScoreSet
 *  @brief simple "inner class" to keep states properties
 *         during SeedScore method
 *  @see Automaton_SeedScore
 */

class SeedScoreSet {

public:
  /// Aho-Corasick Fail function for this state
  int   Fail;
  /// Score reached at this state
  short Score;
  /// Level (distance from the root) in the AC-tree for this state
  short Level;
  /** @brief build a SeedScoreSet
   */
  SeedScoreSet(int Fail = 0, short Score = 0, short Level = 0) : Fail(Fail), Score(Score), Level(Level) {};
};

#define SEEDSCORE_FAIL(s)  (statesSeedScoreSet[(s)].Fail)
#define SEEDSCORE_SCORE(s) (statesSeedScoreSet[(s)].Score)
#define SEEDSCORE_LEVEL(s) (statesSeedScoreSet[(s)].Level)


int automaton::Automaton_SeedScore(const seed & s,
                                   const vector< vector <int> > & matchingmatrix,
                                   const vector< vector<int> > & scoringmatrix,
                                   const int scoringthreehold,
                                   const bool nomerge) {

  // Get the seed span
  int   motif_span  = s.span();
  int * motif       = s.table();

#ifdef ASSERTB
  if (motif_span <= 0) {
    _ERROR("Automaton_SeedScore","negative span");
  }
#endif

  // (0) compute the suffix best score (sbs) :
  // it gives, for a suffix of length i the best suffix score
  // reached : thus is prefix_score[i] + bbs[i] < scoringthreehold
  // then it is not possible to reach the require scoringthreehold.

  vector<int> sbs(motif_span,0);
  sbs[0] = 0;

  // for all suffixes length
  for (int i = 1; i < motif_span; i++){
    int b  = motif[motif_span - i];
    int bs = -INT_INFINITY;
    // FIXME : this part can be precomputed out of the function
    for (int a = 0; a < gv_align_alphabet_size; a++ ) {
      if (MATCHES_AB(a,b)){
        bs = MAX(bs,SCORES_AB(a,b));
      }
    }
    sbs[i] = sbs[i-1] + bs;
#ifdef BUILD
    fprintf(stderr," sbs[%d]=%d\n",i,sbs[i]);
#endif
  }



  // (1) build the trie
  queue<int>  statesNbRemaining;

  // keep each state information (Fail,Rest,Level) inside this table
  vector<SeedScoreSet> statesSeedScoreSet(0);

  // create a first state [0] and put it as the final one
  int Finalstate_I = addNewState(TRUE);
  statesSeedScoreSet.push_back(SeedScoreSet(0,0,0));
  selfLoop(Finalstate_I);

  // create a second state [1]  : it will be the initial one
  int Initstate_I = addNewState();
  statesSeedScoreSet.push_back(SeedScoreSet(0,0,0));
  statesNbRemaining.push(Initstate_I);

  while (!statesNbRemaining.empty()) {

    // current state remaining
    int stateN =  statesNbRemaining.front();
    statesNbRemaining.pop();
    int i        = SEEDSCORE_LEVEL(stateN);
    int final    = (i >= motif_span - 1) ? TRUE : FALSE;

    int b = motif[i];
    for (int a = 0; a < gv_align_alphabet_size; a++) {
      if (MATCHES_AB(a,b)){
        if (SEEDSCORE_SCORE(stateN) + SCORES_AB(a,b) + sbs[motif_span - (i+1)] >= scoringthreehold){
          int StateX = 0;
          if ((!nomerge &&  i < motif_span - 1) || (nomerge && i < motif_span)) {
            StateX = addNewState(final);
            statesSeedScoreSet.push_back(SeedScoreSet(0,SEEDSCORE_SCORE(stateN) + SCORES_AB(a,b),i+1));
            if (!final) {
              statesNbRemaining.push(StateX);
            }
          }
          addNewTransition(a,stateN,StateX);
#ifdef ASSERTB
          assert(DETERMINISTIC(stateN,a));
#endif
        }
      }//prefix penalty
    }//for a
  }//while states remaining

  sbs.clear();

  // (2) add failure links

  // init first states
  int  r  = 1;
  SEEDSCORE_FAIL(r) = r;
  for (int a = 0; a < gv_align_alphabet_size; a++ ) {
    if (!hasTransition(a,r)) {
      addNewTransition(a,r,SEEDSCORE_FAIL(r));
    } else {
      int p = _states[r]._next[a].begin()->_state;
      SEEDSCORE_FAIL(p) = r;
    }
#ifdef ASSERTB
    assert(DETERMINISTIC(r,a));
#endif
  }

  // init last state
  SEEDSCORE_FAIL(0) = 0;


  // intermediate states:  BORDER function from Aho-Corasick*/
  for (int r = 2; r < size(); r++){
    for (int a = 0; a < gv_align_alphabet_size; a++){
      if (!hasTransition(a,r)) {
#ifdef ASSERTB
        assert(DETERMINISTIC(SEEDSCORE_FAIL(r),a));
#endif
        //if (!(_states[r]._final)) //FIXME /*FIXME || nomerge ??*/
        addNewTransition(a,r,_states[SEEDSCORE_FAIL(r)]._next[a].begin()->_state);
#ifdef ASSERTB
        assert(DETERMINISTIC(r,a));
#endif

      } else {

#ifdef ASSERTB
        assert(DETERMINISTIC(r,a));
#endif
        int p = _states[r]._next[a].begin()->_state;
        int s = SEEDSCORE_FAIL(r);
#ifdef ASSERTB
        assert(DETERMINISTIC(s,a));
#endif
        SEEDSCORE_FAIL(p) = _states[s]._next[a].begin()->_state;
      }
    }
  }
  statesSeedScoreSet.clear();
  return 0;
}










/** @class SeedScoreCostSet
 *  @brief simple "inner class" to keep states properties,
 *         during SeedScoreCost method
 *  @see Automaton_SeedScore
 */

class SeedScoreCostSet {

public:
  /// Aho-Corasick Fail function for this state
  int   Fail;
  /// Score reached at this state
  short Score;
  /// Level (distance from the root) in the AC-tree for this state
  short Level;
  /// gives the Cost at this state
  int   Cost;
  /** @brief build a SeedScoreCostSet
   */
  SeedScoreCostSet(int Fail = 0, short Score = 0, short Level = 0,int Cost = 0) : Fail(Fail), Score(Score), Level(Level), Cost(Cost) {};
};

#define SEEDSCORECOST_FAIL(s)  (statesSeedScoreCostSet[(s)].Fail)
#define SEEDSCORECOST_SCORE(s) (statesSeedScoreCostSet[(s)].Score)
#define SEEDSCORECOST_LEVEL(s) (statesSeedScoreCostSet[(s)].Level)
#define SEEDSCORECOST_COST(s)  (statesSeedScoreCostSet[(s)].Cost)

int automaton::Automaton_SeedScoreCost(const seed & s,
                                       const vector< vector <int> > & matchingmatrix,
                                       const vector< vector<int> > & scoringmatrix,
                                       const int scoringthreehold, const bool nomerge,
                                       const vector<int> & costs,
                                       const int cost_threshold) {

  // Get the seed span
  int   motif_span  = s.span();
  int * motif       = s.table();

#ifdef ASSERTB
  if (motif_span <= 0) {
    _ERROR("Automaton_SeedScoreCost","negative span");
  }
#endif

  // (0) compute the suffix best score (sbs) :
  // it gives, for a suffix of length i the best suffix score
  // reached : thus is prefix_score[i] + bbs[i] < scoringthreehold
  // then it is not possible to reach the require scoringthreehold.

  vector<int> sbs(motif_span,0);
  sbs[0] = 0;

  // for all suffixes length
  for (int i = 1; i < motif_span; i++){
    int b  = motif[motif_span - i];
    int bs = -INT_INFINITY;
    // FIXME : this part can be precomputed out of the function
    for (int a = 0; a < gv_align_alphabet_size; a++ ) {
      if (MATCHES_AB(a,b)){
        bs = MAX(bs,SCORES_AB(a,b));
      }
    }
    sbs[i] = sbs[i-1] + bs;
#ifdef BUILD
    fprintf(stderr," sbs[%d]=%d\n",i,sbs[i]);
#endif
  }



  // (1) build the trie
  queue<int>  statesNbRemaining;

  // keep each state information (Fail,Rest,Level) inside this table
  vector<SeedScoreCostSet> statesSeedScoreCostSet(0);

  // create a first state [0] and put it as the final one
  int Finalstate_I = addNewState(TRUE);
  statesSeedScoreCostSet.push_back(SeedScoreCostSet(0,0,0,0));
  selfLoop(Finalstate_I);

  // create a second state [1]  : it will be the initial one
  int Initstate_I = addNewState();
  statesSeedScoreCostSet.push_back(SeedScoreCostSet(0,0,0,0));
  statesNbRemaining.push(Initstate_I);

  // bag reject set [2]
  int RejectBagstate_I   = addNewState();
  statesSeedScoreCostSet.push_back(SeedScoreCostSet(0,0,0,INT_INFINITY));
  selfLoop(RejectBagstate_I);


  while (!statesNbRemaining.empty()) {

    // current state remaining
    int stateN =  statesNbRemaining.front();
    statesNbRemaining.pop();
    int i        = SEEDSCORECOST_LEVEL(stateN);
    int final    = (i >= motif_span - 1) ? TRUE : FALSE;

    int b = motif[i];

    for (int a = 0; a < gv_align_alphabet_size; a++) {
      int cost  = SEEDSCORECOST_COST(stateN) + costs[a];
      if (MATCHES_AB(a,b)){
        if (cost > cost_threshold) {
          addNewTransition(a,stateN,RejectBagstate_I);
        } else {
          if (SEEDSCORECOST_SCORE(stateN) + SCORES_AB(a,b) + sbs[motif_span - (i+1)] >= scoringthreehold){
            int StateX = Finalstate_I; //FIXME
            if ((!nomerge &&  i < motif_span - 1) || (nomerge && i < motif_span)) {
              StateX = addNewState(final);
              statesSeedScoreCostSet.push_back(SeedScoreCostSet(0,SEEDSCORECOST_SCORE(stateN) + SCORES_AB(a,b),i+1,cost));
              if (!final) {
                statesNbRemaining.push(StateX);
              }
            }
            addNewTransition(a,stateN,StateX);
#ifdef ASSERTB
            assert(DETERMINISTIC(stateN,a));
#endif
          }
        }// costs[a] <= cost_threshold
      }//prefix penalty
    }//for a
  }//while states remaining

  sbs.clear();

  // (2) add failure links

  // init first states
  int  r  = 1;
  SEEDSCORECOST_FAIL(r) = r;
  for (int a = 0; a < gv_align_alphabet_size; a++ ) {
    if (!hasTransition(a,r)) {
      addNewTransition(a,r,SEEDSCORECOST_FAIL(r));
    } else {
      int p = _states[r]._next[a].begin()->_state;
      SEEDSCORECOST_FAIL(p) = r;
    }
#ifdef ASSERTB
    assert(DETERMINISTIC(r,a));
#endif
  }

  // init last state
  SEEDSCORECOST_FAIL(0) = 0;


  // intermediate states:  BORDER function from Aho-Corasick*/
  for (int r = 2; r < size(); r++){
    for (int a = 0; a < gv_align_alphabet_size; a++){
      if (!hasTransition(a,r)) {
#ifdef ASSERTB
        assert(DETERMINISTIC(SEEDSCORECOST_FAIL(r),a));
#endif
        //if (!(_states[r]._final)) //FIXME /*FIXME || nomerge ??*/
        int cost = SEEDSCORECOST_COST(r) + costs[a];
        if (cost > cost_threshold)
          addNewTransition(a,r,RejectBagstate_I);
        else
          addNewTransition(a,r,_states[SEEDSCORECOST_FAIL(r)]._next[a].begin()->_state);
#ifdef ASSERTB
        assert(DETERMINISTIC(r,a));
#endif

      } else {

#ifdef ASSERTB
        assert(DETERMINISTIC(r,a));
#endif
        int p = _states[r]._next[a].begin()->_state;
        int s = SEEDSCORECOST_FAIL(r);
#ifdef ASSERTB
        assert(DETERMINISTIC(s,a));
#endif
        SEEDSCORECOST_FAIL(p) = _states[s]._next[a].begin()->_state;
      }
    }
  }
  statesSeedScoreCostSet.clear();
  return 0;
}


















// Gives the probability to reach a final state after nbSteps

double automaton::PrFinal(const int nbSteps) const {
  vector<double> v0(size(),0);
  vector<double> v1(size(),0);

  // Prev
  vector< vector< vector<transition> > > PREV(size(),vector< vector<transition> >(gv_align_alphabet_size, vector<transition>(0)));
  for (int i = 0; i < size(); i++) {
    for (int a = 0; a < gv_align_alphabet_size; a++)
      for (
           vector<transition>::const_iterator
             iter  = _states[i]._next[a].begin();
           iter != _states[i]._next[a].end();
           iter ++ ) {
        PREV[iter->_state][a].push_back(transition(i,iter->_prob));
      }
  }

  // set initial state probability to 1.0
  v1[1] = 1.00;

  // compute probability at step i provided probabilities at step i-1
  for (int k = 0; k < nbSteps; k++) {
    if (k&1){
      for (int nbstate = 0; nbstate < (int)_states.size(); nbstate++) {
        double val =  0.00;
        for (int a = 0; a < gv_align_alphabet_size; a++) {
          for (vector<transition>::iterator iter = PREV[nbstate][a].begin(); iter != PREV[nbstate][a].end(); iter++) {
            val += (iter->_prob) * v0[(iter->_state)];
          }
        }
        v1[nbstate] = val;
      }

    } else {

      for (int nbstate = 0; nbstate < (int)_states.size(); nbstate++) {
        double val =  0.00;
        for (int a = 0; a < gv_align_alphabet_size; a++) {
          for (vector<transition>::iterator iter = PREV[nbstate][a].begin(); iter != PREV[nbstate][a].end(); iter++) {
            val += (iter->_prob) * v1[(iter->_state)];
          }
        }
        v0[nbstate] = val;
      }
    }
  }

  // sum final states probabilities
  double result = 0;
  if (nbSteps & 1) {
    for (int i = 0; i < (int)_states.size(); i++)
      if (_states[i]._final)
        result += v0[i];
  } else {
    for (int i = 0; i < (int)_states.size(); i++)
      if (_states[i]._final)
        result += v1[i];
  }
  v0.clear();
  v1.clear();
  return result;
}


bool  automaton::Lossless(const int nbSteps,
                          const vector<int> & costs,
                          const int cost_threshold,
                          const int Nocc) const {
  vector<int> mincost0(size(),INT_INFINITY/2);
  vector<int> finalcount0(size(),INT_INFINITY/2);

  vector<int> mincost1(size(),INT_INFINITY/2);
  vector<int> finalcount1(size(),INT_INFINITY/2);

  // Prev
  vector< vector< vector<transition> > > PREV(size(),vector< vector<transition> >(gv_align_alphabet_size, vector<transition>(0)));
  for (int i = 0; i < size(); i++) {
    for (int a = 0; a < gv_align_alphabet_size; a++)
      for (vector<transition>::const_iterator
             iter  = _states[i]._next[a].begin();
           iter != _states[i]._next[a].end();
           iter ++ ) {
        PREV[iter->_state][a].push_back(transition(i,iter->_prob));
      }
  }

  // set initial counts
  mincost1[1] = 0;
  finalcount1[1] = 0;

  for (int i= 0; i < nbSteps; i++) {
    if (i&1){
      for (int nbstate = 0; nbstate < (int)_states.size(); nbstate++) {
        int mincost    =  INT_INFINITY/2;
        int finalcount =  INT_INFINITY/2;
        for (int a = 0; a<gv_align_alphabet_size; a++) {
          for (vector<transition>::iterator iter = PREV[nbstate][a].begin(); iter != PREV[nbstate][a].end(); iter++) {
            int cost       = costs[a] + mincost0[(iter->_state)];
            if (cost <= cost_threshold) {
              mincost    = MIN(mincost,cost);
              finalcount = MIN(finalcount, _states[nbstate]._final + finalcount0[(iter->_state)]);
            }
          }
        }
        mincost1[nbstate] = mincost;
        finalcount1[nbstate] = finalcount;
      }
    } else {
      for (int nbstate = 0; nbstate < (int)_states.size(); nbstate++) {
        int mincost    =  INT_INFINITY/2;
        int finalcount =  INT_INFINITY/2;
        for (int a = 0; a<gv_align_alphabet_size; a++) {
          for (vector<transition>::iterator iter = PREV[nbstate][a].begin(); iter != PREV[nbstate][a].end(); iter++) {
            int cost       = costs[a] + mincost1[(iter->_state)];
            if (cost <= cost_threshold) {
              mincost    = MIN(mincost,cost);
              finalcount = MIN(finalcount, _states[nbstate]._final + finalcount1[(iter->_state)]);
            }
          }
        }
        mincost0[nbstate] = mincost;
        finalcount0[nbstate] = finalcount;
      }
    }
  }

  // sum final states probabilities
  bool    bNocc  = true;
  if (nbSteps & 1) {
    for (int i = 0; bNocc && (i<(int)_states.size());i++)
      if (mincost0[i] <= cost_threshold) {
        bNocc = (finalcount0[i] >= Nocc);
      }
  } else {
    for (int i = 0; bNocc && (i<(int)_states.size());i++)
      if (mincost1[i] <= cost_threshold) {
        bNocc = (finalcount1[i] >= Nocc);
      }
  }

  mincost0.clear();
  mincost1.clear();
  finalcount0.clear();
  finalcount1.clear();

#ifdef BUILD
  cerr << "Lossless : " << bNocc << endl;
#endif

  return bNocc;
}



double automaton::PrLossless(const int nbSteps,
                             const vector<int> & costs,
                             const int cost_threshold) const {

  vector< vector<double> > v0 = vector< vector<double> > (size(), vector<double>(cost_threshold+2,0.0));
  vector< vector<double> > v1 = vector< vector<double> > (size(), vector<double>(cost_threshold+2,0.0));

  // Prev
  vector< vector< vector<transition> > > PREV(size(),vector< vector<transition> >(gv_align_alphabet_size, vector<transition>(0)));
  for (int i = 0; i < size(); i++) {
    for (int a = 0; a < gv_align_alphabet_size; a++)
      for (vector<transition>::const_iterator
             iter  = _states[i]._next[a].begin();
           iter != _states[i]._next[a].end();
           iter ++ ) {
        PREV[iter->_state][a].push_back(transition(i,iter->_prob));
      }
  }

  // set initial state probability to 1.0
  v1[1][0] = 1.0;

  // compute probability at step i provided probabilities at step i-1
  for (int i= 0; i < nbSteps; i++) {
    if (i&1){
      for (int nbstate = 0; nbstate < (int)_states.size(); nbstate++) {
        for (int k = 0;k<=cost_threshold;k++)
          v1[nbstate][k] = 0.0;
        for (int kp = 0;kp<=cost_threshold+1;kp++) {
          for (int a = 0; a<gv_align_alphabet_size; a++) {
            int k  = kp + costs[a];
            for (vector<transition>::iterator iter = PREV[nbstate][a].begin(); iter != PREV[nbstate][a].end(); iter++)
              v1[nbstate][MIN(k,cost_threshold+1)] +=  (iter->_prob) * v0[(iter->_state)][kp];
          }
        }
      }
    } else {
      for (int nbstate = 0; nbstate < (int)_states.size(); nbstate++) {
        for (int k=0;k<=cost_threshold;k++)
          v0[nbstate][k] = 0.0;
        for (int kp=0;kp<=cost_threshold+1;kp++) {
          for (int a = 0; a<gv_align_alphabet_size; a++) {
            int k  = kp + costs[a];
            for (vector<transition>::iterator iter = PREV[nbstate][a].begin(); iter != PREV[nbstate][a].end(); iter++)
              v0[nbstate][MIN(k,cost_threshold+1)] +=  (iter->_prob) * v1[(iter->_state)][kp];
          }
        }
      }
    }
  }

  // sum final states probabilities
  double  PrFinalLessThanThreshold = 0.0;
  double  PrNonFinalLessThanThreshold = 0.0;
  double  PrFullLessThanThreshold = 0.0;
  if (nbSteps & 1) {
    for (int i=0;i<(int)size();i++) {
      for (int k = 0; k <= cost_threshold; k++)
        PrFullLessThanThreshold += v0[i][k];
      if(_states[i]._final)
        for (int k = 0; k <= cost_threshold; k++)
          PrFinalLessThanThreshold += v0[i][k];
      else
        for (int k = 0; k <= cost_threshold; k++)
          PrNonFinalLessThanThreshold += v0[i][k];
    }
  } else {
    for (int i=0;i<(int)size();i++) {
      for (int k = 0; k <= cost_threshold; k++)
        PrFullLessThanThreshold += v1[i][k];
      if(_states[i]._final)
        for (int k = 0; k <= cost_threshold; k++)
          PrFinalLessThanThreshold += v1[i][k];
      else
        for (int k = 0; k <= cost_threshold; k++)
          PrNonFinalLessThanThreshold += v1[i][k];
    }
  }
  v0.clear();
  v1.clear();

#ifdef BUILD
  cerr << "PrLossless : " << PrFinalLessThanThreshold << "/" << PrFullLessThanThreshold << "=" << (PrFinalLessThanThreshold / PrFullLessThanThreshold) <<endl;
#endif

  return PrFinalLessThanThreshold / PrFullLessThanThreshold;
}





void automaton::GenSeq(const int len) const {
  int q = 1;
  for (int l = 0; l < len; l++) {
    double p = (double)rand()/(double)RAND_MAX;
    int q_new = 0,a = 0;
    double psum = 0;
    for (a = 0; a < gv_align_alphabet_size; a++){
      for (vector<transition>::const_iterator iter = _states[q]._next[a].begin(); iter != _states[q]._next[a].end(); iter++) {
        double pr = iter->_prob;
        q_new = iter->_state;
        psum +=  pr;
        if (psum >= p)
          goto eol;
      }
    }
  eol:
    q = q_new;
    cout << a;
  }
  cout << endl;
}

void automaton::GenSeq(vector<int> & alignment) const {
  int q = 1;
  for (unsigned l = 0; l < alignment.size(); l++) {
    double p = (double)rand()/(double)RAND_MAX;
    int q_new = 0,a = 0;
    double psum = 0;
    for (a = 0; a < gv_align_alphabet_size; a++) {
      for (vector<transition>::const_iterator iter = _states[q]._next[a].begin(); iter != _states[q]._next[a].end(); iter++) {
        double pr = iter->_prob;
        q_new = iter->_state;
        psum +=  pr;
        if (psum >= p)
          goto eol;
      }
    }
  eol:
    q = q_new;
    alignment[l] = a;
  }
}
