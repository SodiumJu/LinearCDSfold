#ifndef _AC_MATCHER_H_
#define _AC_MATCHER_H_

#include <cstring>
#include <list>
#include <queue>
#include <unordered_map>

namespace oligodesign {

template<class T, int S> 
class ACMatcher {
public:
  ACMatcher();
  ~ACMatcher();
  void AddString(const T *, int);
  bool MakeTree();
  int Search(const T *, int);
  int Search(const T *, int, int);
  void SearchPositions(const T *text,int idx,int start,int end, std::unordered_map<int,int> &mp);
  int NumEndingInLast3(const T *, int);
  int NumEndingInLastK(const T *, int, int, int);

private:
  struct ACState
  {
    int id;
    int depth;
    int count;
    T *output;
    struct ACState *fail;
    struct ACState *transitions[S];
  };
  typedef struct ACState ACState_t;
      
  int newstate;
  ACState_t *zerostate;
  std::list<T *> cached_strings;
      
  void ACFree(ACState_t *);
};

template<class T, int S>
ACMatcher<T,S>::ACMatcher()
{
  zerostate = 0;
  newstate = 0;
}

template<class T, int S>
ACMatcher<T,S>::~ACMatcher()
{
  if (zerostate != 0) {
    for (int i = 0; i < S ;i++) {
      if (zerostate->transitions[i] != 0 &&
          zerostate->transitions[i]->id > 0) {
        ACFree(zerostate->transitions[i]);
      }
    }
  }
  // clean up cached strings
  class std::list<T *>::iterator itr;
  for( itr = cached_strings.begin(); itr != cached_strings.end(); itr++ ) {
    if (*itr != 0) {
      delete [] *itr;
    }
  }
}

template<class T, int S>
void ACMatcher<T,S>::ACFree(ACState_t *state) {
  for (int i = 0; i < S ;i++) {
    if (state->transitions[i] != 0) {
      ACFree(state->transitions[i]);
    }
  }
  delete(state);
}

template<class T, int S>
void ACMatcher<T,S>::AddString(const T *text, int n) {
  if (text == 0) return;
      
  // first cache the string locally
  cached_strings.push_back(new T[n]);
  T *cached_item = cached_strings.back();
  for(int i=0; i<n; i++) {
    cached_item[i] = text[i];
  }
      
  ACState_t *state = 0, *s = 0;
  int j = 0;

  if (zerostate == 0) {
    zerostate = new ACState_t;
    newstate = 1;
    zerostate->id = 0;
    zerostate->depth = 0;
    zerostate->output = 0;
    memset(zerostate->transitions, 0, S * sizeof(ACState_t *));
  }

  state = zerostate;

  // As long as we have transitions follow them
  while (j < n
         && (s = state->transitions[(int)*(text+j)]) != 0) {
    state = s;
    ++j;
  }

  if (j == n) {
    s->count = s->count + 1;
    return;
  }

  while (j < n) {
    // Create new state
    s = new ACState_t;
    s->id = newstate++;
    s->depth = state->depth + 1;
    memset(s->transitions, 0, S * sizeof(ACState_t *));
    // Create transition
    state->transitions[(int) * (text + j)] = s;
    state = s;
    s->output = 0;
    s->count = 1;
    ++j;
  }

  s->output = cached_strings.back();
  return;
}

template<class T, int S>
bool ACMatcher<T,S>::MakeTree() {
  if (cached_strings.empty()) return false;
  
  std::queue<ACState_t *> state_queue;
  ACState_t *state, *s, *r;
  int i;

  // Set all FAIL transition of 0 state to point to itself
  for (i = 0; i < S; i++) {
    if (zerostate->transitions[i] == 0)
      zerostate->transitions[i] = zerostate;
    // Construct fail()
    else {
      state_queue.push(zerostate->transitions[i]);
      zerostate->transitions[i]->fail = zerostate;
    }
  }

  // Set fail() for depth > 0
  while (!state_queue.empty()) {
    r = state_queue.front();
    state_queue.pop();
    for (i = 0; i < S; i++) {
      if ((s = r->transitions[i]) == 0)
        continue;
      state_queue.push(s);
      state = r->fail;
      while (state->transitions[i] == 0)
        state = state->fail;
      s->fail = state->transitions[i];
      // Join outputs missing
    }
  }
  return true;
}

template<class T, int S>
int ACMatcher<T,S>::Search(const T *text, int n) {
  if (text == 0) return 0;
  return Search(text, 0, n - 1);
}

template<class T, int S>
int ACMatcher<T,S>::Search(const T *text, int start, int end) {
  if (cached_strings.empty()) return 0;
  
  ACState_t *state = zerostate;
  int num_matches = 0;

  for (int j = start; j <= end; j++) {
    while (state->transitions[(int)*(text + j)] == 0) {
      state = state->fail;
    }
    state = state->transitions[(int)*(text + j)];
    if (state->output != 0) {
      //return state->id;
      num_matches += state->count;
      //return j;
    }
  }
  return num_matches;
  //return -1;
}

template<class T, int S>
void ACMatcher<T,S>::SearchPositions(const T *text,int idx,int start,int end, std::unordered_map<int,int> &mp) {
    if (cached_strings.empty() || text == 0)
        return ;
    ACState_t *state = zerostate;
    for (int j = start; j <= end; j++) {
        while (state->transitions[(int)*(text + j)] == 0)
            state = state->fail;
        state = state->transitions[(int)*(text + j)];
        if (state->output != 0) {
            for (int k = 0; k < state->count; k++) {
               if(mp.find(idx+j) == mp.end()){
                mp[j+idx] = state->depth;
               }
            }
        }
    }
    return;
}


template<class T, int S>
int ACMatcher<T,S>::NumEndingInLast3(const T *text, int n) {
  if (text == 0) return 0;
  return NumEndingInLastK(text, 0, n - 1, 3);
}

template<class T, int S>
int ACMatcher<T,S>::NumEndingInLastK(const T *text, int start, int end, int k) {
  if (cached_strings.empty()) return 0;
  ACState_t *state = zerostate;
  int num_matches = 0;

  for (int j = start; j <= end; j++) {
    while (state->transitions[(int)*(text + j)] == 0) {
      state = state->fail;
    }
    state = state->transitions[(int)*(text + j)];
    if (state->output != 0) {
      //return state->id;
      if (j > end - k) {
        num_matches += state->count;
      }
      //return j;
    }
  }
  return num_matches;
  //return -1;
}

} // end namespace oligodesign

#endif /* _AC_MATCHER_H_ */
