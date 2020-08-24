
/**
 * Info about this game from a run of bluffcounter.cpp:
 *
 * storing 2 doubles per iapair + 2 doubles per infoset
 * want size of index to be at least double # infosets (then 2 doubles per index size)
 *
 * for 1,1:
 *   p1 infosets 12288, p2 infosets 12288, total: 24576
 *   infoset actions pairs, p1 p2 total = 24570 24570 49140
 */


#include <iostream>
#include <string>
#include <cstring>
#include <map>
#include <cassert>
#include <cstdlib>
#include <ctime>
#include <set>

#include "bluff.h"
#include "infosetstore.h"
#include "sys/time.h"

#define LOC(b,r,c)  b[r*3 + c]

using namespace std;

// global variables
InfosetStore iss;
string filepref = "scratch/";
unsigned long long iter;
int iscWidth = 0;
double cpWidth = 10.0;
double nextCheckpoint = cpWidth;
unsigned long long nodesTouched = 0;
unsigned long long ntNextReport = 1000000;  // nodes touched base timing
unsigned long long ntMultiplier = 2;  // nodes touched base timing

// key is roll, value is # of time it shows up. Used only when determining chance outcomes
map<int,int> outcomes;

// probability of this chance move. indexed as 0-5 or 0-20
// note: CO/co stand for "chance outcome"
static int numChanceOutcomes1 = 0;
static int numChanceOutcomes2 = 0;
static double * chanceProbs1 = NULL;
static double * chanceProbs2 = NULL;
static int * chanceOutcomes1 = NULL;
static int * chanceOutcomes2 = NULL;
static int * bids = NULL;

static StopWatch stopwatch;

double getChanceProb(int player, int outcome)
{
  // outcome >= 1, so must subtract 1 from it
  int co = (player == 1 ? numChanceOutcomes1 : numChanceOutcomes2);
  assert(outcome-1 >= 0 && outcome-1 < co);
  double * cp = (player == 1 ? chanceProbs1 : chanceProbs2);
  return cp[outcome-1];
}

int numChanceOutcomes(int player)
{
  return (player == 1 ? numChanceOutcomes1 : numChanceOutcomes2);
}

void unrankco(int i, int * roll, int player)
{
  int num = 0;
  int * chanceOutcomes = (player == 1 ? chanceOutcomes1 : chanceOutcomes2);
  num = chanceOutcomes[i];

  assert(num > 0);

  int numDice = (player == 1 ? P1DICE : P2DICE);

  for (int j = numDice-1; j >= 0; j--)
  {
    roll[j] = num % 10;
    num /= 10;
  }
}


void initBids()
{
  bids = new int[BLUFFBID-1];
  int nextWildDice = 1;
  int idx = 0;
  for (int dice = 1; dice <= P1DICE + P2DICE; dice++)
  {
    for (int face = 1; face <= DIEFACES-1; face++)
    {
      bids[idx] = dice*10 + face;
      idx++;
    }

    if (dice % 2 == 1) {
      bids[idx] = nextWildDice*10 + DIEFACES;
      idx++;
      nextWildDice++;
    }
  }

  for(; nextWildDice <= (P1DICE+P2DICE); nextWildDice++) {
    bids[idx] = nextWildDice*10 + DIEFACES;
    idx++;
  }

  assert(idx == BLUFFBID-1);

}

unsigned long long getInfosetKey(GameState & gs, int player, unsigned long long bidseq)
{
  unsigned long long infosetkey = bidseq;
  infosetkey <<= iscWidth;
  if (player == 1)
  {
    infosetkey |= gs.p1roll;
    infosetkey <<= 1;
  }
  else if (player == 2)
  {
    infosetkey |= gs.p2roll;
    infosetkey <<= 1;
    infosetkey |= 1;
  }

  return infosetkey;
}

void getInfoset(GameState & gs, int player, unsigned long long bidseq, Infoset & is, unsigned long long & infosetkey, int actionshere)
{
  infosetkey = getInfosetKey(gs, player, bidseq);
  bool ret = iss.get(infosetkey, is, actionshere, 0);
  assert(ret);
}

int ceiling_log2(int val)
{
  int exp = 1, num = 2;
  do {
    if (num > val) { return exp; }
    num *= 2;
    exp++;
  }
  while (true);
}

int intpow(int x, int y)
{
  if (y == 0) return 1;
  return x * intpow(x, y-1);
}

void nextRoll(int * roll, int size)
{
  for (int i = size-1; i >= 0; i--)
  {
    // Try to increment if by 1.
    if (roll[i] < DIEFACES) {
      // if possible, do it and then set everything to the right back to 1
      roll[i]++;
      for (int j = i+1; j < size; j++)
        roll[j] = 1;

      return;
    }
  }
}

int getRollBase10(int * roll, int size)
{
  int multiplier = 1;
  int val = 0;
  for (int i = size-1; i >= 0; i--)
  {
    val += roll[i]*multiplier;
    multiplier *= 10;
  }

  return val;
}

void determineChanceOutcomes(int player)
{
  int dice = (player == 1 ? P1DICE : P2DICE);
  int rolls[dice];
  for (int r = 0; r < dice; r++) rolls[r] = 1;
  outcomes.clear();

  int permutations = intpow(DIEFACES, dice);
  int p;

  for (p = 0; p < permutations; p++) {

    // first, make a copy
    int rollcopy[dice];
    memcpy(rollcopy, rolls, dice*sizeof(int));

    // now sort
    bubsort(rollcopy, dice);

    // now convert to an integer in base 10
    int key = getRollBase10(rollcopy, dice);

    // now increment the counter for this key in the map
    outcomes[key] += 1;

    // next roll
    nextRoll(rolls, dice);
  }

  assert(p == permutations);

  int & numChanceOutcomes = (player == 1 ? numChanceOutcomes1 : numChanceOutcomes2);
  double* & chanceProbs = (player == 1 ? chanceProbs1 : chanceProbs2);
  int* & chanceOutcomes = (player == 1 ? chanceOutcomes1 : chanceOutcomes2);

  // now, transfer the map keys to the array
  numChanceOutcomes = outcomes.size();
  chanceProbs = new double[numChanceOutcomes];
  chanceOutcomes = new int[numChanceOutcomes];

  map<int,int>::iterator iter;
  int idx = 0;
  for (iter = outcomes.begin(); iter != outcomes.end(); iter++) {
    chanceOutcomes[idx] = iter->first;
    idx++;
  }

  bubsort(chanceOutcomes, numChanceOutcomes);

  for (int c = 0; c < numChanceOutcomes; c++) {
    int key = chanceOutcomes[c];
    chanceProbs[c] = static_cast<double>(outcomes[key]) / static_cast<double>(permutations);
    //cout << "player " << player << " roll " << key << " prob " << chanceProbs[c] << endl;
  }
}

void init()
{
  assert(bids == NULL);

  cout << "Initializing Bluff globals..." << endl;

  seedCurMicroSec();

  determineChanceOutcomes(1);
  determineChanceOutcomes(2);

  // iscWidth if the number of bits needed to encode the chance outcome in the integer
  int maxChanceOutcomes = (numChanceOutcomes1 > numChanceOutcomes2 ? numChanceOutcomes1 : numChanceOutcomes2);
  iscWidth = ceiling_log2(maxChanceOutcomes);

  initBids();

  cout << "Globals are: " << numChanceOutcomes1 << " " << numChanceOutcomes2 << " " << iscWidth << endl;
}


void newInfoset(Infoset & is, int actions)
{
  is.actionshere = actions;
  is.lastUpdate = 0;

  for (int i = 0; i < actions; i++)
  {
    is.cfr[i] = 0.0;
    is.totalMoveProbs[i] = 0.0;
    is.curMoveProbs[i] = 1.0 / actions;
  }
}

bool terminal(GameState & gs)
{
  return (gs.curbid == BLUFFBID);
}

// a bid is from 1 to 12, for example
void convertbid(int & dice, int & face, int bid)
{
  if (P1DICE == 1 && P2DICE == 1)
  {
    dice = (bid - 1) / DIEFACES + 1;
    face = bid % DIEFACES;
    if (face == 0) face = DIEFACES;

    assert(dice >= 1 && dice <= 2);
    assert(face >= 1 && face <= DIEFACES);
  }
  else
  {
    // stored in an array.
    int size = (P1DICE+P2DICE)*DIEFACES;
    assert((bid-1) >= 0 && (bid-1) < size);

    dice = bids[bid-1] / 10;
    face = bids[bid-1] % 10;
  }
}

void getRoll(int * roll, int chanceOutcome, int player)
{
  unrankco(chanceOutcome-1, roll, player);
}

int countMatchingDice(const GameState & gs, int player, int face)
{
  int roll[3] = {0,0,0};
  int matchingDice = 0;
  int dice = (player == 1 ? P1DICE : P2DICE);

  if (dice == 1)
  {
    if (player == 1)
      roll[1] = gs.p1roll;
    else if (player == 2)
      roll[1] = gs.p2roll;
  }
  else if (dice == 2)
  {
    if (player == 1)
      unrankco(gs.p1roll-1, roll, 1);
    else if (player == 2)
      unrankco(gs.p2roll-1, roll, 2);
  }

  for (int i = 0; i < 3; i++)
    if (roll[i] == face || roll[i] == DIEFACES)
      matchingDice++;

  return matchingDice;
}

int whowon(int bid, int bidder, int callingPlayer, int p1roll, int p2roll, int & delta)
{
  int dice = 0, face = 0;
  convertbid(dice, face, bid);

  assert(bidder != callingPlayer);

  // get the dice

  int p1rollArr[P1DICE];
  int p2rollArr[P2DICE];

  unrankco(p1roll-1, p1rollArr, 1);
  unrankco(p2roll-1, p2rollArr, 2);

  // now check the number of matches

  int matching = 0;

  for (int i = 0; i < P1DICE; i++)
    if (p1rollArr[i] == face || p1rollArr[i] == DIEFACES)
      matching++;

  for (int j = 0; j < P2DICE; j++)
    if (p2rollArr[j] == face || p2rollArr[j] == DIEFACES)
      matching++;

  delta = matching - dice;
  if (delta < 0) delta *= -1;

  if (matching >= dice)
  {
    return bidder;
  }
  else
  {
    return callingPlayer;
  }
}

int whowon(GameState & gs, int & delta)
{
  int bidder = 3 - gs.callingPlayer;
  return whowon(gs.prevbid, bidder, gs.callingPlayer, gs.p1roll, gs.p2roll, delta);
}

int whowon(GameState & gs)
{
  int bidder = 3 - gs.callingPlayer;
  int delta = 0;
  return whowon(gs.prevbid, bidder, gs.callingPlayer, gs.p1roll, gs.p2roll, delta);
}

double payoff(int winner, int player, int delta)
{
  // first thing: if it's an exact match, calling player loses 1 die
  // may as well set delta to 1 in this case
  if (delta == 0) delta = 1;

  double p1payoff = 0.0;

  if
    (P1DICE == 1 && P2DICE == 1) return (winner == player ? 1.0 : -1.0);
  else
  {
    assert(false);
  }

  return (player == 1 ? p1payoff : -p1payoff);
}

// In these functions "delta" represents the number of dice the bid is off by (not relevant for (1,1))
// Returns payoff for Liar's Dice (delta always equal to 1)
double payoff(int winner, int player)
{
  return payoff(winner, player, 1);
}

// this is the function called by all the algorithms.
// Now set to use the delta
double payoff(GameState & gs, int player)
{
  int delta = 0;
  int winner = whowon(gs, delta);
  return payoff(winner, player, delta);
}

double payoff(int bid, int bidder, int callingPlayer, int p1roll, int p2roll, int player)
{
  int delta = 0;
  int winner = whowon(bid, bidder, callingPlayer, p1roll, p2roll, delta);
  return payoff(winner, player);
}

void report(string filename, double totaltime, double bound, double conv)
{
  filename = filepref + filename;
  cout << "Reporting to " + filename + " ... " << endl;
  ofstream outf(filename.c_str(), ios::app);
  outf << iter << " " << totaltime << " " << bound << " " << conv << " " << nodesTouched << endl;
  outf.close();
}

void dumpInfosets(string prefix)
{
  string filename = filepref + prefix + "." + to_string(iter) + ".dat";
  cout << "Dumping infosets to " + filename + " ... " << endl;
  iss.dumpToDisk(filename);
}

// not even sure what I used this "meta data" for, if I ever used it....
void dumpMetaData(string prefix, double totaltime)
{
  string filename = filepref + prefix + "." + to_string(iter) + ".dat";
  cout << "Dumping metadeta to " + filename + " ... " << endl;

  ofstream outf(filename.c_str(), ios::binary);
  if (!outf.is_open()) {
    cerr << "Could not open meta data file for writing." << endl;
    return;
  }

  outf.write(reinterpret_cast<const char *>(&iter), sizeof(iter));
  outf.write(reinterpret_cast<const char *>(&nodesTouched), sizeof(nodesTouched));
  outf.write(reinterpret_cast<const char *>(&ntNextReport), sizeof(ntNextReport));
  outf.write(reinterpret_cast<const char *>(&ntMultiplier), sizeof(ntMultiplier));
  outf.write(reinterpret_cast<const char *>(&totaltime), sizeof(totaltime));

  outf.close();
}

void loadMetaData(std::string filename)
{
  ifstream inf(filename.c_str(), ios::binary);
  if (!inf.is_open()) {
    cerr << "Could not open meta data file." << endl;
    return;
  }

  double totaltime = 0;

  inf.read(reinterpret_cast<char *>(&iter), sizeof(iter));
  inf.read(reinterpret_cast<char *>(&nodesTouched), sizeof(nodesTouched));
  inf.read(reinterpret_cast<char *>(&ntNextReport), sizeof(ntNextReport));
  inf.read(reinterpret_cast<char *>(&ntMultiplier), sizeof(ntMultiplier));
  inf.read(reinterpret_cast<char *>(&totaltime), sizeof(totaltime));

  inf.close();
}

// Does a recursive tree walk setting up the information sets, creating the initial strategies
void initInfosets(GameState & gs, int player, int depth, unsigned long long bidseq)
{
  if (terminal(gs))
    return;

  // check for chance nodes
  if (gs.p1roll == 0)
  {
    for (int i = 1; i <= numChanceOutcomes1; i++)
    {
      GameState ngs = gs;
      ngs.p1roll = i;

      initInfosets(ngs, player, depth+1, bidseq);
    }

    return;
  }
  else if (gs.p2roll == 0)
  {
    for (int i = 1; i <= numChanceOutcomes2; i++)
    {
      GameState ngs = gs;
      ngs.p2roll = i;

      initInfosets(ngs, player, depth+1, bidseq);
    }

    return;
  }

  int maxBid = (gs.curbid == 0 ? BLUFFBID-1 : BLUFFBID);
  int actionshere = maxBid - gs.curbid;

  assert(actionshere > 0);
  Infoset is;
  newInfoset(is, actionshere);

  for (int i = gs.curbid+1; i <= maxBid; i++)
  {
    if (depth == 2 && i == (gs.curbid+1)) {
      cout << "InitTrees. iss stats = " << iss.getStats() << endl;
    }

    GameState ngs = gs;
    ngs.prevbid = gs.curbid;
    ngs.curbid = i;
    ngs.callingPlayer = player;
    unsigned long long newbidseq = bidseq;
    newbidseq |= (1ULL << (BLUFFBID-i));

    initInfosets(ngs, (3-player), depth+1, newbidseq);
  }

  unsigned infosetkey = 0;
  infosetkey = bidseq;
  infosetkey <<= iscWidth;
  if (player == 1)
  {
    infosetkey |= gs.p1roll;
    infosetkey <<= 1;
    iss.put(infosetkey, is, actionshere, 0);
  }
  else if (player == 2)
  {
    infosetkey |= gs.p2roll;
    infosetkey <<= 1;
    infosetkey |= 1;
    iss.put(infosetkey, is, actionshere, 0);
  }
}

void initInfosets()
{
  unsigned long long bidseq = 0;

  GameState gs;

  cout << "Initialize info set store..." << endl;
  // # doubles in total, size of index (must be at least # infosets)
  // 2 doubles per iapair + 2 per infoset =
  if (P1DICE == 1 && P2DICE == 1 && DIEFACES == 6)
    iss.init(147432, 100000);
  else
  {
    cerr << "initInfosets not defined for this PXDICE + DIEFACES" << endl;
  }

  assert(iss.getSize() > 0);

  cout << "Initializing info sets..." << endl;
  stopwatch.reset();
  initInfosets(gs, 1, 0, bidseq);

  cout << "time taken = " << stopwatch.stop() << " seconds." << endl;
  iss.stopAdding();

  cout << "Final iss stats: " << iss.getStats() << endl;
  stopwatch.reset();

  string filename = filepref + "iss.initial.dat";

  cout << "Dumping information sets to " << filename << endl;
  iss.dumpToDisk(filename);
}

